#include "Runner.h"

#include <utility>

#include "src/common/Logger.h"
#include "src/common/PrintingFunctions.h"
#include "src/eigendecompositor/EigendecompositorConstructor.h"
#include "src/entities/magnetic_susceptibility/worker/WorkerConstructor.h"
#include "src/space/optimization/OptimizedSpaceConstructor.h"

namespace runner {

Runner::Runner(model::ModelInput model) :
    Runner(
        std::move(model),
        common::physical_optimization::OptimizationList(),
        quantum::linear_algebra::FactoriesList()) {}

Runner::Runner(
    model::ModelInput model,
    common::physical_optimization::OptimizationList optimizationList) :
    Runner(
        std::move(model),
        std::move(optimizationList),
        quantum::linear_algebra::FactoriesList()) {}

Runner::Runner(
    model::ModelInput model,
    quantum::linear_algebra::FactoriesList dataStructuresFactories) :
    Runner(
        std::move(model),
        common::physical_optimization::OptimizationList(),
        std::move(dataStructuresFactories)) {}

Runner::Runner(
    model::ModelInput model,
    common::physical_optimization::OptimizationList optimizationList,
    quantum::linear_algebra::FactoriesList dataStructuresFactories) :
    consistentModelOptimizationList_(std::move(model), std::move(optimizationList)),
    dataStructuresFactories_(std::move(dataStructuresFactories)),
    space_(space::optimization::OptimizedSpaceConstructor::construct(
        consistentModelOptimizationList_,
        dataStructuresFactories)),
    eigendecompositor_(eigendecompositor::EigendecompositorConstructor::construct(
        consistentModelOptimizationList_,
        dataStructuresFactories_
        )) {}

const space::Space& runner::Runner::getSpace() const {
    return space_;
}

void Runner::BuildSpectra() {
    eigendecompositor_->BuildSpectra(
        consistentModelOptimizationList_.getOperatorsForExplicitConstruction(),
        consistentModelOptimizationList_.getDerivativeOperatorsForExplicitConstruction(),
        getSpace());
}

std::optional<std::reference_wrapper<const Matrix>>
Runner::getMatrix(common::QuantityEnum quantity_enum) {
    return getEigendecompositor()->getMatrix(quantity_enum);
}

const Spectrum& Runner::getSpectrum(common::QuantityEnum quantity_enum) {
    return getEigendecompositor()->getSpectrum(quantity_enum)->get();
}

std::optional<std::shared_ptr<const model::operators::Operator>>
Runner::getOperator(common::QuantityEnum quantity_enum) const {
    return getModel().getOperator(quantity_enum);
}

const lexicographic::IndexConverter& Runner::getIndexConverter() const {
    return getModel().getIndexConverter();
}

std::optional<std::shared_ptr<const model::operators::Operator>> Runner::getOperatorDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol) const {
    return getModel().getOperatorDerivative(quantity_enum, symbol);
}

const Spectrum& Runner::getSpectrumDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol) {
    return getEigendecompositor()->getSpectrumDerivative(quantity_enum, symbol)->get();
}

std::optional<std::reference_wrapper<const Matrix>> Runner::getMatrixDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol) {
    return getEigendecompositor()->getMatrixDerivative(quantity_enum, symbol);
}

void Runner::BuildMuSquaredWorker() {
    auto magnetic_susceptibility_worker =
        magnetic_susceptibility::worker::WorkerConstructor::construct(
            getModel(),
            getEigendecompositor(),
            getDataStructuresFactories());

    magnetic_susceptibility_controller_ = magnetic_susceptibility::MagneticSusceptibilityController(
        std::move(magnetic_susceptibility_worker));

    if (experimental_values_worker_.has_value()) {
        magnetic_susceptibility_controller_.value().initializeExperimentalValues(
            experimental_values_worker_.value());
    }
}

void Runner::initializeExperimentalValues(
    const std::shared_ptr<magnetic_susceptibility::ExperimentalValuesWorker>& experimental_values_worker) {
    if (experimental_values_worker_.has_value()) {
        throw std::invalid_argument("Experimental values have been already initialized");
    }

    experimental_values_worker_ = experimental_values_worker;

    if (magnetic_susceptibility_controller_.has_value()) {
        magnetic_susceptibility_controller_.value().initializeExperimentalValues(
            experimental_values_worker_.value());
    }
}

void Runner::initializeExperimentalValues(
    const std::vector<magnetic_susceptibility::ValueAtTemperature>& experimental_data,
    magnetic_susceptibility::ExperimentalValuesEnum experimental_quantity_type,
    double number_of_centers_ratio,
    magnetic_susceptibility::WeightingSchemeEnum weightingSchemeEnum) {

    auto experimental_values_worker =
        std::make_shared<magnetic_susceptibility::ExperimentalValuesWorker>(
            experimental_data,
            experimental_quantity_type,
            number_of_centers_ratio,
            weightingSchemeEnum);

    initializeExperimentalValues(experimental_values_worker);
}

std::map<model::symbols::SymbolName, double> Runner::calculateTotalDerivatives() {
    // TODO: it is awful. Fix it somehow, this code should be moved from Runner.

    initializeDerivatives();

    std::map<model::symbols::SymbolName, double> answer;

    for (const auto& changeable_symbol :
         getSymbolicWorker().getChangeableNames(model::symbols::J)) {
        auto derivative_vector = dataStructuresFactories_.createVector();
        for (const auto& subspectrum :
             getSpectrumDerivative(common::Energy, changeable_symbol).blocks) {
            derivative_vector->concatenate_with(subspectrum.raw_data);
        }
        auto derivative_map = std::map<
            common::QuantityEnum,
            std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>();
        derivative_map[common::Energy] = std::move(derivative_vector);
        double value = getMagneticSusceptibilityController().calculateTotalDerivative(
            model::symbols::J,
            std::move(derivative_map));
        answer[changeable_symbol] = value;
        common::Logger::debug("d(Loss function)/d{}: {}", changeable_symbol.get_name(), value);
    }

    for (const auto& changeable_symbol :
         getSymbolicWorker().getChangeableNames(model::symbols::D)) {
        auto derivative_vector = dataStructuresFactories_.createVector();
        for (const auto& subspectrum :
             getSpectrumDerivative(common::Energy, changeable_symbol).blocks) {
            derivative_vector->concatenate_with(subspectrum.raw_data);
        }
        auto derivative_map = std::map<
            common::QuantityEnum,
            std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>();
        derivative_map[common::Energy] = std::move(derivative_vector);
        double value = getMagneticSusceptibilityController().calculateTotalDerivative(
            model::symbols::D,
            std::move(derivative_map));
        answer[changeable_symbol] = value;
        common::Logger::debug("d(Loss function)/d{}: {}", changeable_symbol.get_name(), value);
    }

    for (const auto& changeable_symbol :
         getSymbolicWorker().getChangeableNames(model::symbols::g_factor)) {
        auto map = std::map<
            common::QuantityEnum,
            std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>();
        if (consistentModelOptimizationList_.getModel().is_g_sz_squared_derivatives_initialized()) {
            auto derivative_vector = dataStructuresFactories_.createVector();
            for (const auto& subspectrum :
                 getSpectrumDerivative(common::gSz_total_squared, changeable_symbol).blocks) {
                derivative_vector->concatenate_with(subspectrum.raw_data);
            }
            map[common::gSz_total_squared] = std::move(derivative_vector);
        }
        double value = getMagneticSusceptibilityController().calculateTotalDerivative(
            model::symbols::g_factor,
            std::move(map));
        answer[changeable_symbol] = value;
        common::Logger::debug("d(Loss function)/d{}: {}", changeable_symbol.get_name(), value);
    }

    // Theta calculation:
    if (!getSymbolicWorker().getChangeableNames(model::symbols::Theta).empty()) {
        model::symbols::SymbolName Theta_name =
            getSymbolicWorker().getChangeableNames(model::symbols::Theta)[0];
        auto empty_map = std::map<
            common::QuantityEnum,
            std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>();
        double value = getMagneticSusceptibilityController().calculateTotalDerivative(
            model::symbols::Theta,
            std::move(empty_map));
        answer[Theta_name] = value;
        common::Logger::debug("d(Loss function)/d{}: {}", Theta_name.get_name(), value);
    }
    common::Logger::separate(2, common::debug);

    return answer;
}

void Runner::minimizeResidualError(
    std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver> solver) {
    std::vector<model::symbols::SymbolName> changeable_names =
        getSymbolicWorker().getChangeableNames();
    std::vector<double> changeable_values;
    changeable_values.reserve(changeable_names.size());
    for (const model::symbols::SymbolName& name : changeable_names) {
        changeable_values.push_back(getSymbolicWorker().getValueOfName(name));
    }

    if (solver->doesGradientsRequired()) {
        initializeDerivatives();
    }

    // This function should calculate residual error and derivatives at a point of changeable_values
    std::function<double(const std::vector<double>&, std::vector<double>&, bool)> oneStepFunction =
        [this, capture0 = std::cref(changeable_names)](
            const std::vector<double>& changeable_values,
            std::vector<double>& gradient,
            bool isGradientRequired) {
            return stepOfRegression(capture0, changeable_values, gradient, isGradientRequired);
        };

    common::preRegressionPrint(
        consistentModelOptimizationList_.getOperatorsForExplicitConstruction(),
        consistentModelOptimizationList_.getDerivativeOperatorsForExplicitConstruction());

    solver->optimize(oneStepFunction, changeable_values);

    common::postRegressionPrint(
        changeable_names,
        changeable_values,
        getMagneticSusceptibilityController().calculateResidualError());
}

double Runner::stepOfRegression(
    const std::vector<model::symbols::SymbolName>& changeable_names,
    const std::vector<double>& changeable_values,
    std::vector<double>& gradient,
    bool isGradientRequired) {
    // At first, update actual values in SymbolicWorker:
    for (size_t i = 0; i < changeable_names.size(); ++i) {
        consistentModelOptimizationList_.setNewValueToChangeableSymbol(
            changeable_names[i],
            changeable_values[i]);
    }

    common::stepOfRegressionStartPrint(changeable_names, changeable_values);

    // (Re)build Eigendecompositor:
    BuildSpectra();
    // (Re)build MuSquaredWorker:
    BuildMuSquaredWorker();

    // Calculate residual error and write it to external variable:
    double residual_error = getMagneticSusceptibilityController().calculateResidualError();

    if (isGradientRequired) {
        // Calculate derivatives...
        std::map<model::symbols::SymbolName, double> map_gradient = calculateTotalDerivatives();
        for (size_t i = 0; i < changeable_names.size(); ++i) {
            // ...and write it to external variable:
            gradient[i] = map_gradient[changeable_names[i]];
        }
    }

    common::stepOfRegressionFinishPrint(residual_error);

    return residual_error;
}

void Runner::initializeDerivatives() {
    consistentModelOptimizationList_.InitializeDerivatives();
}

const magnetic_susceptibility::MagneticSusceptibilityController&
Runner::getMagneticSusceptibilityController() {
    if (!magnetic_susceptibility_controller_.has_value()) {
        BuildMuSquaredWorker();
    }
    return magnetic_susceptibility_controller_.value();
}

const model::symbols::SymbolicWorker& Runner::getSymbolicWorker() const {
    return getModel().getSymbolicWorker();
}

const model::Model& Runner::getModel() const {
    return consistentModelOptimizationList_.getModel();
}

const common::physical_optimization::OptimizationList& Runner::getOptimizationList() const {
    return consistentModelOptimizationList_.getOptimizationList();
}

quantum::linear_algebra::FactoriesList Runner::getDataStructuresFactories() const {
    return dataStructuresFactories_;
}

const std::unique_ptr<eigendecompositor::AbstractEigendecompositor>&
Runner::getEigendecompositor() {
    if (!eigendecompositor_->BuildSpectraWasCalled()) {
        BuildSpectra();
    }

    return eigendecompositor_;
}
}  // namespace runner
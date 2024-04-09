#ifndef SPINNER_RUNNER_H
#define SPINNER_RUNNER_H

#include <cmath>
#include <cstdint>
#include <memory>
#include <optional>
#include <utility>

#include "ConsistentModelOptimizationList.h"
#include "src/common/Quantity.h"
#include "src/eigendecompositor/AbstractEigendecompositor.h"
#include "src/entities/data_structures/FactoriesList.h"
#include "src/entities/magnetic_susceptibility/MagneticSusceptibilityController.h"
#include "src/entities/matrix/Matrix.h"
#include "src/entities/spectrum/Spectrum.h"
#include "src/nonlinear_solver/AbstractNonlinearSolver.h"
#include "src/space/Space.h"

namespace runner {
class Runner {
  public:
    Runner(
        model::ModelInput model,
        common::physical_optimization::OptimizationList optimizationList,
        quantum::linear_algebra::FactoriesList dataStructuresFactories);
    // constructor for Runner with default algebra package:
    Runner(
        model::ModelInput model,
        common::physical_optimization::OptimizationList optimizationList);
    // constructor for Runner with no optimizations:
    Runner(model::ModelInput model, quantum::linear_algebra::FactoriesList dataStructuresFactories);
    // constructor for Runner with no optimizations and default algebra package:
    explicit Runner(model::ModelInput model);

    // CHIT OPERATIONS
    void initializeExperimentalValues(
        const std::vector<magnetic_susceptibility::ValueAtTemperature>& experimental_data,
        magnetic_susceptibility::ExperimentalValuesEnum experimental_quantity_type,
        double number_of_centers_ratio,
        magnetic_susceptibility::WeightingSchemeEnum weightingSchemeEnum =
            magnetic_susceptibility::per_point);
    void initializeExperimentalValues(
        const std::shared_ptr<magnetic_susceptibility::ExperimentalValuesWorker>& experimental_values_worker);
    std::map<model::symbols::SymbolName, double> calculateTotalDerivatives();
    void minimizeResidualError(std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver>);

    const lexicographic::IndexConverter& getIndexConverter() const;
    std::optional<std::shared_ptr<const model::operators::Operator>>
        getOperator(common::QuantityEnum) const;
    const space::Space& getSpace() const;
    const Spectrum& getSpectrum(common::QuantityEnum);
    std::optional<std::reference_wrapper<const Matrix>> getMatrix(common::QuantityEnum);
    std::optional<std::shared_ptr<const model::operators::Operator>>
    getOperatorDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const;
    const Spectrum&
    getSpectrumDerivative(common::QuantityEnum, const model::symbols::SymbolName&);
    std::optional<std::reference_wrapper<const Matrix>>
    getMatrixDerivative(common::QuantityEnum, const model::symbols::SymbolName&);
    const magnetic_susceptibility::MagneticSusceptibilityController&
    getMagneticSusceptibilityController();
    const model::symbols::SymbolicWorker& getSymbolicWorker() const;

    quantum::linear_algebra::FactoriesList getDataStructuresFactories() const;

  private:
    ConsistentModelOptimizationList consistentModelOptimizationList_;
    const space::Space space_;
    quantum::linear_algebra::FactoriesList dataStructuresFactories_;
    std::unique_ptr<eigendecompositor::AbstractEigendecompositor> eigendecompositor_;

    // SPECTRUM OPERATIONS
    void BuildSpectra();
    // CHIT OPERATIONS
    void BuildMuSquaredWorker();

    void initializeDerivatives();

    const model::Model& getModel() const;
    const common::physical_optimization::OptimizationList& getOptimizationList() const;

    const std::unique_ptr<eigendecompositor::AbstractEigendecompositor>& getEigendecompositor();

    double stepOfRegression(
        const std::vector<model::symbols::SymbolName>&,
        const std::vector<double>&,
        std::vector<double>&,
        bool isGradientRequired);

    std::optional<magnetic_susceptibility::MagneticSusceptibilityController>
        magnetic_susceptibility_controller_;
    std::optional<std::shared_ptr<magnetic_susceptibility::ExperimentalValuesWorker>>
        experimental_values_worker_;

};
}  // namespace runner

#endif  //SPINNER_RUNNER_H

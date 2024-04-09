#include "Model.h"

#include <utility>

#include "src/model/operators/terms/LocalSSquaredOneCenterTerm.h"
#include "src/model/operators/terms/ScalarProductTerm.h"
#include "src/model/operators/terms/SzSzOneCenterTerm.h"
#include "src/model/operators/terms/SzSzTwoCenterTerm.h"

namespace model {
Model::Model(ModelInput modelInput) :
    numericalWorker_(modelInput.getSymbolicWorker(), modelInput.getMults().size()),
    converter_(modelInput.getMults()) {
    operators_map_[common::Energy] = std::make_shared<operators::Operator>();
    InitializeSSquared();
    if (getSymbolicWorker().isGFactorInitialized()) {
        InitializeGSzSquared();
        InitializeGSzSquaredDerivatives();
    }

    if (getSymbolicWorker().isZFSInitialized()) {
        InitializeZeroFieldSplitting();
        InitializeZeroFieldSplittingDerivative();
    }

    if (getSymbolicWorker().isIsotropicExchangeInitialized()) {
        InitializeIsotropicExchange();
        InitializeIsotropicExchangeDerivatives();
    }
}

void Model::InitializeSSquared() {
    if (operators_history_.s_squared) {
        return;
    }

    operators_map_[common::S_total_squared] =
        std::make_shared<operators::Operator>(operators::Operator::s_squared(converter_));

    operators_history_.s_squared = true;
}

void Model::InitializeGSzSquared() {
    if (operators_history_.g_sz_squared) {
        return;
    }

    operators_map_[common::gSz_total_squared] =
        std::make_shared<operators::Operator>(operators::Operator::g_sz_squared(
            converter_,
            getNumericalWorker().getGGParameters().first,
            getNumericalWorker().getGGParameters().second));

    operators_history_.g_sz_squared = true;
}

void Model::InitializeIsotropicExchange() {
    if (operators_history_.isotropic_exchange_in_hamiltonian) {
        return;
    }
    operators_map_[common::Energy]->emplace_back(
        std::make_unique<const operators::ScalarProductTerm>(
            converter_,
            getNumericalWorker().getIsotropicExchangeParameters()));
    operators_history_.isotropic_exchange_in_hamiltonian = true;
}

void Model::InitializeZeroFieldSplitting() {
    if (operators_history_.zfs_in_hamiltonian) {
        return;
    }
    auto D_parameters = getNumericalWorker().getZFSParameters().first;
    operators_map_[common::Energy]->emplace_back(
        std::make_unique<const operators::LocalSSquaredOneCenterTerm>(
            converter_,
            D_parameters,
            -1.0 / 3.0));
    operators_map_[common::Energy]->emplace_back(
        std::make_unique<const operators::SzSzOneCenterTerm>(converter_, D_parameters));
    operators_history_.zfs_in_hamiltonian = true;
}

void Model::InitializeIsotropicExchangeDerivatives() {
    if (operators_history_.isotropic_exchange_derivatives) {
        return;
    }

    for (const auto& symbol : getSymbolicWorker().getChangeableNames(symbols::SymbolTypeEnum::J)) {
        operators::Operator operator_derivative = operators::Operator();
        operator_derivative.emplace_back(
            std::make_unique<const operators::ScalarProductTerm>(
                converter_,
                getNumericalWorker().constructIsotropicExchangeDerivativeParameters(symbol)));
        derivatives_map_[{common::Energy, symbol}] =
            std::make_shared<operators::Operator>(std::move(operator_derivative));
    }

    operators_history_.isotropic_exchange_derivatives = true;
}

void Model::InitializeZeroFieldSplittingDerivative() {
    if (operators_history_.zfs_derivative) {
        return;
    }
    for (const auto& symbol : getSymbolicWorker().getChangeableNames(symbols::SymbolTypeEnum::D)) {
        operators::Operator operator_derivative = operators::Operator();
        auto derivative_parameters = getNumericalWorker().constructZFSDerivativeParameters(symbol);
        operator_derivative.emplace_back(
            std::make_unique<const operators::LocalSSquaredOneCenterTerm>(
                converter_,
                derivative_parameters,
                -1.0 / 3.0));
        operator_derivative.emplace_back(
            std::make_unique<const operators::SzSzOneCenterTerm>(
                converter_,
                derivative_parameters));
        derivatives_map_[{common::Energy, symbol}] =
            std::make_shared<operators::Operator>(std::move(operator_derivative));
    }
}

void Model::InitializeGSzSquaredDerivatives() {
    if (operators_history_.g_sz_squared_derivatives) {
        return;
    }

    for (const auto& symbol :
         getSymbolicWorker().getChangeableNames(symbols::SymbolTypeEnum::g_factor)) {
        operators::Operator operator_derivative = operators::Operator();
        auto pair_of_parameters = getNumericalWorker().constructGGDerivativeParameters(symbol);
        operator_derivative.emplace_back(
            std::make_unique<operators::SzSzOneCenterTerm>(converter_, pair_of_parameters.first));
        operator_derivative.emplace_back(
            std::make_unique<const operators::SzSzTwoCenterTerm>(
                converter_,
                pair_of_parameters.second,
                2));
        // this two from summation in Submatrix: \sum_{a=1}^N \sum_{b=a+1}^N
        derivatives_map_[{common::gSz_total_squared, symbol}] =
            std::make_shared<operators::Operator>(std::move(operator_derivative));
    }

    operators_history_.g_sz_squared_derivatives = true;
}

const lexicographic::IndexConverter& Model::getIndexConverter() const {
    return converter_;
}

std::optional<std::shared_ptr<const operators::Operator>>
Model::getOperator(common::QuantityEnum quantity_enum) const {
    if (operators_map_.contains(quantity_enum)) {
        return operators_map_.at(quantity_enum);
    }
    return std::nullopt;
}

const std::map<common::QuantityEnum, std::shared_ptr<operators::Operator>>& Model::getOperators() {
    return operators_map_;
}

std::optional<std::shared_ptr<const operators::Operator>> Model::getOperatorDerivative(
    common::QuantityEnum quantity_enum,
    const symbols::SymbolName& symbol) const {
    if (derivatives_map_.contains({quantity_enum, symbol})) {
        return derivatives_map_.at({quantity_enum, symbol});
    }
    return std::nullopt;
}

const std::
    map<std::pair<common::QuantityEnum, symbols::SymbolName>, std::shared_ptr<operators::Operator>>&
    Model::getOperatorDerivatives() {
    return derivatives_map_;
}

bool Model::is_s_squared_initialized() const {
    return operators_history_.s_squared;
}

bool Model::is_isotropic_exchange_derivatives_initialized() const {
    return operators_history_.isotropic_exchange_derivatives;
}

bool Model::is_g_sz_squared_initialized() const {
    return operators_history_.g_sz_squared;
}

bool Model::is_g_sz_squared_derivatives_initialized() const {
    return operators_history_.g_sz_squared_derivatives;
}

bool Model::is_zero_field_splitting_initialized() const {
    return operators_history_.zfs_in_hamiltonian;
}

const symbols::SymbolicWorker& Model::getSymbolicWorker() const {
    return numericalWorker_.getSymbolicWorker();
}

symbols::NumericalWorker& Model::getNumericalWorker() {
    return numericalWorker_;
}

const symbols::NumericalWorker& Model::getNumericalWorker() const {
    return numericalWorker_;
}

void Model::setNewValueToChangeableSymbol(
    const symbols::SymbolName& symbol_name,
    double new_value) {
    numericalWorker_.setNewValueToChangeableSymbol(symbol_name, new_value);
}
}  // namespace model
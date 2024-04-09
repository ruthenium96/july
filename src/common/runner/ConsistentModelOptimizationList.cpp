#include "ConsistentModelOptimizationList.h"

#include <utility>

namespace {
void checkMultiplicitiesGroupConsistence(
    const std::vector<spin_algebra::Multiplicity>& mults,
    const group::Group& group);
bool checkSymbolNamesGroupElementConsistence(
    const std::function<std::optional<model::symbols::SymbolName>(size_t, size_t)>& getter,
    group::Permutation element);
bool checkSymbolNamesGroupElementConsistence(
    const std::function<model::symbols::SymbolName(size_t)>& getter,
    group::Permutation element);
void checkAllSymbolNamesGroupConsistence(
    const model::symbols::SymbolicWorker& symbols,
    const group::Group& group);
void checkModelOptimizationListConsistence(
    const model::Model& model,
    const common::physical_optimization::OptimizationList& optimizationList);
}  // namespace

namespace runner {

const model::Model& ConsistentModelOptimizationList::getModel() const {
    return model_;
}

const common::physical_optimization::OptimizationList&
ConsistentModelOptimizationList::getOptimizationList() const {
    return optimizationList_;
}

ConsistentModelOptimizationList::ConsistentModelOptimizationList(
    model::ModelInput modelInput,
    common::physical_optimization::OptimizationList optimizationList) :
    model_(model::Model(std::move(modelInput))),
    optimizationList_(std::move(optimizationList)) {
    // this function throw an exception if ModelInput and OptimizationList are inconsistent
    checkModelOptimizationListConsistence(model_, optimizationList_);

    operators_for_explicit_construction_[common::Energy] =
        model_.getOperator(common::Energy).value();
    if (getOptimizationList().isSSquaredTransformed()) {
        return;
    }
    if (getModel().is_g_sz_squared_initialized()
        && (!getModel().getSymbolicWorker().isAllGFactorsEqual()
            || getModel().getSymbolicWorker().isZFSInitialized())) {
        operators_for_explicit_construction_[common::gSz_total_squared] =
            model_.getOperator(common::gSz_total_squared).value();
        // TODO: when there is no Sz <-> -Sz symmetry, also \sum g_aS_{az} required
    } else {
        // TODO: some tests want to calculate S^2 values. Can we fix it?
        operators_for_explicit_construction_[common::S_total_squared] =
            model_.getOperator(common::S_total_squared).value();
    }
}

void ConsistentModelOptimizationList::InitializeDerivatives() {
    if (derivatives_for_explicit_construction_.empty()) {
        for (const auto& [pair, shared_ptr] : model_.getOperatorDerivatives()) {
            derivatives_for_explicit_construction_[pair] = shared_ptr;
        }
    }
}

void ConsistentModelOptimizationList::setNewValueToChangeableSymbol(
    const model::symbols::SymbolName& symbol_name,
    double new_value) {
    model_.setNewValueToChangeableSymbol(symbol_name, new_value);
}

const std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
ConsistentModelOptimizationList::getOperatorsForExplicitConstruction() const {
    return operators_for_explicit_construction_;
}

const std::map<
    std::pair<common::QuantityEnum, model::symbols::SymbolName>,
    std::shared_ptr<const model::operators::Operator>>&
ConsistentModelOptimizationList::getDerivativeOperatorsForExplicitConstruction() const {
    return derivatives_for_explicit_construction_;
}
}  // namespace runner

namespace {
void checkModelOptimizationListConsistence(
    const model::Model& model,
    const common::physical_optimization::OptimizationList& optimizationList) {
    // Symmetrizer can be applied for all types of Hamiltonian terms...
    for (const auto& group : optimizationList.getGroupsToApply()) {
        // ...if spins invariant to group elements:
        checkMultiplicitiesGroupConsistence(model.getIndexConverter().get_mults(), group);
        // ...if Hamiltonian terms invariant to group elements:
        checkAllSymbolNamesGroupConsistence(model.getSymbolicWorker(), group);
    }
    if (optimizationList.isSSquaredTransformed()) {
        // S2-transformation can be applied only for HDvV-Hamiltionian:
        if (model.is_zero_field_splitting_initialized()) {
            throw std::invalid_argument("S2-transformation cannot be applied to ZFS-Hamiltonian");
        }
        // TODO: check if we actually can use S2-transformation even if g-factors are not equal
        // S2-transformation can be applied only if all g-factors are equal:
        if (model.getSymbolicWorker().isGFactorInitialized()
            && !model.getSymbolicWorker().isAllGFactorsEqual()) {
            throw std::invalid_argument(
                "S2-transformation cannot be applied to system with different g-factors");
        }
    }
}

void checkMultiplicitiesGroupConsistence(
    const std::vector<spin_algebra::Multiplicity>& mults,
    const group::Group& group) {
    // TODO: split code above (maybe rewrite checkSymbolNamesGroupElementConsistence with templates?)
    if (mults.size() != group.size_of_permutations()) {
        throw std::length_error(
            "The size of group elements does not equal to the number of spins.");
    }
    for (const auto& el : group.getElements()) {
        std::vector<spin_algebra::Multiplicity> permutated_mults(mults);
        for (size_t i = 0; i < group.size_of_permutations(); ++i) {
            permutated_mults[i] = mults[el[i]];
        }
        if (permutated_mults != mults) {
            throw std::invalid_argument("Group permutes centers with different multiplicities.");
        }
    }
}

bool checkSymbolNamesGroupElementConsistence(
    const std::function<std::optional<model::symbols::SymbolName>(size_t, size_t)>& getter,
    group::Permutation element) {
    // Construct here both initialSymbols and permutatedSymbols:
    std::vector<std::vector<std::optional<model::symbols::SymbolName>>> initialSymbols(
        element.size(),
        std::vector<std::optional<model::symbols::SymbolName>>(element.size(), std::nullopt));
    auto permutatedSymbols = initialSymbols;
    for (size_t i = 0; i < element.size(); ++i) {
        for (size_t j = 0; j < element.size(); ++j) {
            auto symbolName = getter(i, j);
            initialSymbols[i][j] = symbolName;
            permutatedSymbols[element[i]][element[j]] = symbolName;
        }
    }
    return initialSymbols == permutatedSymbols;
}

bool checkSymbolNamesGroupElementConsistence(
    const std::function<model::symbols::SymbolName(size_t)>& getter,
    group::Permutation element) {
    std::vector<model::symbols::SymbolName> initialSymbols(element.size());
    auto permutatedSymbols = initialSymbols;

    for (size_t i = 0; i < element.size(); ++i) {
        auto symbolName = getter(i);
        initialSymbols[i] = symbolName;
        permutatedSymbols[element[i]] = symbolName;
    }
    return initialSymbols == permutatedSymbols;
}

void checkAllSymbolNamesGroupConsistence(
    const model::symbols::SymbolicWorker& symbols,
    const group::Group& group) {
    if (symbols.isIsotropicExchangeInitialized()) {
        for (const auto& element : group.getElements()) {
            std::function<std::optional<model::symbols::SymbolName>(size_t, size_t)>
                getterFunction = [ObjectPtr = &symbols](size_t i, size_t j) {
                    return ObjectPtr->getIsotropicExchangeSymbolName(i, j);
                };
            if (!checkSymbolNamesGroupElementConsistence(getterFunction, element)) {
                throw std::invalid_argument(
                    "Isotropic exchange symbols do not match applied symmetries");
                // TODO: should we rename this exception?
            }
        }
    }
    if (symbols.isGFactorInitialized()) {
        for (const auto& element : group.getElements()) {
            std::function<model::symbols::SymbolName(size_t)> getterFunction =
                [ObjectPtr = &symbols](size_t i) { return ObjectPtr->getGFactorSymbolName(i); };
            if (!checkSymbolNamesGroupElementConsistence(getterFunction, element)) {
                throw std::invalid_argument("g factor symbols do not match applied symmetries");
                // TODO: should we rename this exception?
            }
        }
    }
}
}  // namespace
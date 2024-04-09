#include "SymbolicWorker.h"

#include <cmath>
#include <utility>

namespace model::symbols {
SymbolicWorker::SymbolicWorker(size_t number_of_spins) : number_of_spins_(number_of_spins) {}

SymbolicWorker& SymbolicWorker::assignSymbolToIsotropicExchange(
    const SymbolName& symbol_name,
    size_t center_a,
    size_t center_b) {
    if (!symbolic_isotropic_exchanges_.has_value()) {
        symbolic_isotropic_exchanges_ = std::vector<std::vector<std::optional<SymbolName>>>(
            number_of_spins_,
            std::vector<std::optional<SymbolName>>(number_of_spins_, std::nullopt));
    }
    if (center_b == center_a) {
        throw std::invalid_argument("Isotropic exchange takes place between different centers");
    }
    if (symbolic_isotropic_exchanges_.value()[center_a][center_b].has_value()) {
        throw std::invalid_argument("This parameter has been already specified");
    }

    if (getSymbolProperty(symbol_name).type_enum.has_value()
        && getSymbolProperty(symbol_name).type_enum != SymbolTypeEnum::J) {
        throw std::invalid_argument(
            symbol_name.get_name() + " has been already specified as not J parameter");
    }
    if (!getSymbolProperty(symbol_name).type_enum.has_value()) {
        modifySymbolProperty(symbol_name).type_enum = SymbolTypeEnum::J;
    }

    symbolic_isotropic_exchanges_.value()[center_a][center_b] = symbol_name;
    symbolic_isotropic_exchanges_.value()[center_b][center_a] = symbol_name;

    return *this;
}

SymbolicWorker&
SymbolicWorker::assignSymbolToGFactor(const SymbolName& symbol_name, size_t center_a) {
    if (getSymbolProperty(symbol_name).type_enum.has_value()
        && getSymbolProperty(symbol_name).type_enum != SymbolTypeEnum::g_factor) {
        throw std::invalid_argument(
            symbol_name.get_name() + " has been already specified as not g factor parameter");
    }
    if (!getSymbolProperty(symbol_name).type_enum.has_value()) {
        modifySymbolProperty(symbol_name).type_enum = SymbolTypeEnum::g_factor;
    }

    if (symbolic_g_factors_.empty()) {
        symbolic_g_factors_.resize(number_of_spins_);
    }

    symbolic_g_factors_[center_a] = symbol_name;

    return *this;
}

SymbolicWorker& SymbolicWorker::assignSymbolToTheta(const SymbolName& symbol_name) {
    if (getSymbolProperty(symbol_name).type_enum.has_value()
        && getSymbolProperty(symbol_name).type_enum != SymbolTypeEnum::Theta) {
        throw std::invalid_argument(
            symbol_name.get_name() + " has been already specified as not Theta parameter");
    }
    if (!getSymbolProperty(symbol_name).type_enum.has_value()) {
        modifySymbolProperty(symbol_name).type_enum = SymbolTypeEnum::Theta;
    }

    symbolic_Theta_ = symbol_name;

    return *this;
}

SymbolicWorker&
SymbolicWorker::assignSymbolToZFSNoAnisotropy(const SymbolName& symbol_name, size_t center_a) {
    if (getSymbolProperty(symbol_name).type_enum.has_value()
        && getSymbolProperty(symbol_name).type_enum != SymbolTypeEnum::D) {
        throw std::invalid_argument(
            symbol_name.get_name() + " has been already specified as not D parameter");
    }
    if (!getSymbolProperty(symbol_name).type_enum.has_value()) {
        modifySymbolProperty(symbol_name).type_enum = SymbolTypeEnum::D;
    }

    if (!symbolic_ZFS_.has_value()) {
        symbolic_ZFS_ = std::vector<std::optional<ZFSSymbols>>(number_of_spins_, std::nullopt);
    }

    symbolic_ZFS_.value()[center_a] = {symbol_name, std::nullopt};

    return *this;
}

SymbolName SymbolicWorker::addSymbol(
    const std::string& name_string,
    double initial_value,
    bool is_changeable,
    std::optional<SymbolTypeEnum> type_enum) {
    if (name_string.empty()) {
        throw std::invalid_argument("Name of Symbol should be non-empty");
    }
    SymbolName symbol_name(name_string);
    if (symbolsProperties_.find(symbol_name) != symbolsProperties_.end()) {
        throw std::invalid_argument(symbol_name.get_name() + " name has been already initialized");
    }
    symbolsProperties_[symbol_name] = {is_changeable, type_enum};
    symbolsValues_[symbol_name] = initial_value;
    return symbol_name;
}

std::vector<SymbolName> SymbolicWorker::getAllNames(SymbolTypeEnum type_enum) const {
    std::vector<SymbolName> answer;
    for (const auto& [symbol_name, symbol_data] : symbolsProperties_) {
        if (symbol_data.type_enum == type_enum) {
            answer.push_back(symbol_name);
        }
    }
    return answer;
}

std::vector<SymbolName> SymbolicWorker::getChangeableNames(SymbolTypeEnum type_enum) const {
    std::vector<SymbolName> answer;
    for (const auto& [symbol_name, symbol_data] : symbolsProperties_) {
        if (symbol_data.is_changeable && symbol_data.type_enum == type_enum) {
            answer.push_back(symbol_name);
        }
    }
    return answer;
}

bool SymbolicWorker::isAllGFactorsEqual() const {
    if (!isGFactorInitialized()) {
        throw std::length_error("g factor parameters have not been initialized");
    }
    for (size_t i = 1; i < symbolic_g_factors_.size(); ++i) {
        if (symbolic_g_factors_[i] != symbolic_g_factors_[0]) {
            return false;
        }
    }
    return true;
}

std::vector<SymbolName> SymbolicWorker::getChangeableNames() const {
    std::vector<SymbolName> answer;
    for (const auto& [symbol_name, symbol_data] : symbolsProperties_) {
        if (symbol_data.is_changeable) {
            answer.push_back(symbol_name);
        }
    }
    return answer;
}

double SymbolicWorker::getValueOfName(const SymbolName& symbol_name) const {
    if (symbolsValues_.find(symbol_name) == symbolsValues_.end()) {
        throw std::invalid_argument(symbol_name.get_name() + " name has not been initialized");
    }
    return symbolsValues_.at(symbol_name);
}

void SymbolicWorker::setNewValueToChangeableSymbol(
    const SymbolName& symbol_name,
    double new_value) {
    if (!getSymbolProperty(symbol_name).is_changeable) {
        throw std::invalid_argument("Cannot change value of unchangeable symbol");
    }

    symbolsValues_.at(symbol_name) = new_value;
}

bool SymbolicWorker::isIsotropicExchangeInitialized() const {
    return symbolic_isotropic_exchanges_.has_value();
}

bool SymbolicWorker::isGFactorInitialized() const {
    // TODO: check, that all symbolic g factors were initialized
    return !symbolic_g_factors_.empty();
}

bool SymbolicWorker::isThetaInitialized() const {
    return symbolic_Theta_.has_value();
}

std::optional<SymbolName> SymbolicWorker::getIsotropicExchangeSymbolName(size_t i, size_t j) const {
    if (!symbolic_isotropic_exchanges_.has_value()) {
        throw std::length_error("Isotropic exchange parameters have not been initialized");
    }
    return symbolic_isotropic_exchanges_.value()[i][j];
}

SymbolName SymbolicWorker::getGFactorSymbolName(size_t i) const {
    if (!isGFactorInitialized()) {
        throw std::length_error("G factor parameters have not been initialized");
    }
    return symbolic_g_factors_[i];
}

std::optional<ZFSSymbols> SymbolicWorker::getZFSSymbolNames(size_t i) const {
    return symbolic_ZFS_.value()[i];
}

std::optional<SymbolName> SymbolicWorker::getThetaSymbolName() const {
    return symbolic_Theta_;
}

bool SymbolicWorker::isZFSInitialized() const {
    return symbolic_ZFS_.has_value();
}

SymbolicWorker::SymbolProperty SymbolicWorker::getSymbolProperty(const SymbolName& symbol_name) const {
    if (symbolsProperties_.find(symbol_name) == symbolsProperties_.end()) {
        throw std::invalid_argument(symbol_name.get_name() + " name has not been initialized");
    }
    return symbolsProperties_.at(symbol_name);
}

SymbolicWorker::SymbolProperty& SymbolicWorker::modifySymbolProperty(const SymbolName& symbol_name) {
    if (symbolsProperties_.find(symbol_name) == symbolsProperties_.end()) {
        throw std::invalid_argument(symbol_name.get_name() + " name has not been initialized");
    }
    return symbolsProperties_.at(symbol_name);
}
}  // namespace model::symbols
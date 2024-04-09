#include "ModelInput.h"

namespace model {

ModelInput::ModelInput(const std::vector<spin_algebra::Multiplicity>& mults) :
    symbols_(mults.size()),
    mults_(mults) {}

const symbols::SymbolicWorker& ModelInput::getSymbolicWorker() const {
    return symbols_;
}

const std::vector<spin_algebra::Multiplicity>& ModelInput::getMults() const {
    return mults_;
}

ModelInput& ModelInput::assignSymbolToIsotropicExchange(
    const symbols::SymbolName& symbol_name,
    size_t center_a,
    size_t center_b) {
    symbols_.assignSymbolToIsotropicExchange(symbol_name, center_a, center_b);
    return *this;
}

ModelInput&
ModelInput::assignSymbolToGFactor(const symbols::SymbolName& symbol_name, size_t center_a) {
    symbols_.assignSymbolToGFactor(symbol_name, center_a);
    return *this;
}

ModelInput&
ModelInput::assignSymbolToZFSNoAnisotropy(const symbols::SymbolName& symbol_name, size_t center_a) {
    symbols_.assignSymbolToZFSNoAnisotropy(symbol_name, center_a);
    return *this;
}

ModelInput& ModelInput::assignSymbolToTheta(const symbols::SymbolName& symbol_name) {
    symbols_.assignSymbolToTheta(symbol_name);
    return *this;
}

symbols::SymbolName ModelInput::addSymbol(
    const std::string& name_string,
    double initial_value,
    bool is_changeable,
    std::optional<symbols::SymbolTypeEnum> type_enum) {
    return symbols_.addSymbol(name_string, initial_value, is_changeable, type_enum);
}

}  // namespace model
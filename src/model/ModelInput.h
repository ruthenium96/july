#ifndef SPINNER_MODELINPUT_H
#define SPINNER_MODELINPUT_H

#include <vector>

#include "src/common/Quantity.h"
#include "src/model/operators/Operator.h"
#include "src/model/symbols/SymbolicWorker.h"
namespace model {
// class ModelInput is responsible for user input: SymbolicWorker and Multiplicities
class ModelInput {
  public:
    explicit ModelInput(const std::vector<spin_algebra::Multiplicity>& mults);

    const symbols::SymbolicWorker& getSymbolicWorker() const;
    const std::vector<spin_algebra::Multiplicity>& getMults() const;

    ModelInput& assignSymbolToIsotropicExchange(
        const symbols::SymbolName& symbol_name,
        size_t center_a,
        size_t center_b);
    ModelInput& assignSymbolToGFactor(const symbols::SymbolName& symbol_name, size_t center_a);
    ModelInput& assignSymbolToZFSNoAnisotropy(const symbols::SymbolName& symbol_name, size_t center_a);
    ModelInput& assignSymbolToTheta(const symbols::SymbolName& symbol_name);

    symbols::SymbolName addSymbol(
        const std::string& name_string,
        double initial_value,
        bool is_changeable = true,
        std::optional<symbols::SymbolTypeEnum> type_enum = std::nullopt);

    bool operator==(const ModelInput& rhs) const = default;
    bool operator!=(const ModelInput& rhs) const = default;

  private:
    symbols::SymbolicWorker symbols_;
    std::vector<spin_algebra::Multiplicity> mults_;
};
}  // namespace model

#endif  //SPINNER_MODELINPUT_H

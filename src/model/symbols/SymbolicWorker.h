#ifndef SPINNER_SYMBOLICWORKER_H
#define SPINNER_SYMBOLICWORKER_H

#include <map>
#include <memory>
#include <optional>
#include <stdexcept>

#include "SymbolName.h"
#include "src/group/Group.h"
#include "src/model/NumericalParameters.h"

namespace model::symbols {

enum SymbolTypeEnum { J, g_factor, Theta, D };

struct ZFSSymbols {
    SymbolName D;
    std::optional<SymbolName> E;
    bool operator==(const ZFSSymbols& rhs) const = default;
    bool operator!=(const ZFSSymbols& rhs) const = default;
};

class SymbolicWorker {
  public:
    explicit SymbolicWorker(size_t number_of_spins);

    bool isAllGFactorsEqual() const;
    bool isGFactorInitialized() const;
    bool isIsotropicExchangeInitialized() const;
    bool isThetaInitialized() const;
    bool isZFSInitialized() const;

    std::optional<SymbolName> getIsotropicExchangeSymbolName(size_t i, size_t j) const;
    SymbolName getGFactorSymbolName(size_t i) const;
    std::optional<SymbolName> getThetaSymbolName() const;
    std::optional<ZFSSymbols> getZFSSymbolNames(size_t i) const;

    SymbolicWorker& assignSymbolToIsotropicExchange(
        const SymbolName& symbol_name,
        size_t center_a,
        size_t center_b);
    SymbolicWorker& assignSymbolToGFactor(const SymbolName& symbol_name, size_t center_a);
    SymbolicWorker& assignSymbolToZFSNoAnisotropy(const SymbolName& symbol_name, size_t center_a);
    SymbolicWorker& assignSymbolToTheta(const SymbolName& symbol_name);

    SymbolName addSymbol(
        const std::string& name_string,
        double initial_value,
        bool is_changeable = true,
        std::optional<SymbolTypeEnum> type_enum = std::nullopt);

    void setNewValueToChangeableSymbol(const SymbolName& symbol_name, double new_value);

    std::vector<SymbolName> getAllNames(SymbolTypeEnum type_enum) const;
    std::vector<SymbolName> getChangeableNames(SymbolTypeEnum type_enum) const;
    std::vector<SymbolName> getChangeableNames() const;
    double getValueOfName(const SymbolName& symbol_name) const;

    bool operator==(const SymbolicWorker& rhs) const = default;
    bool operator!=(const SymbolicWorker& rhs) const = default;

  private:
    size_t number_of_spins_;

    struct SymbolProperty {
        bool is_changeable;
        std::optional<SymbolTypeEnum> type_enum;
        bool operator==(const SymbolProperty& rhs) const = default;
        bool operator!=(const SymbolProperty& rhs) const = default;
    };

    std::optional<std::vector<std::vector<std::optional<SymbolName>>>>
        symbolic_isotropic_exchanges_;
    std::vector<SymbolName> symbolic_g_factors_;
    std::optional<SymbolName> symbolic_Theta_;
    std::optional<std::vector<std::optional<ZFSSymbols>>> symbolic_ZFS_;

    std::map<SymbolName, SymbolProperty> symbolsProperties_;
    std::map<SymbolName, double> symbolsValues_;

  public:
    SymbolProperty getSymbolProperty(const SymbolName& symbol_name) const;

  private:
    SymbolProperty& modifySymbolProperty(const SymbolName& symbol_name);
};

}  // namespace model::symbols
#endif  //SPINNER_SYMBOLICWORKER_H

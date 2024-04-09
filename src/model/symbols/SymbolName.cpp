#include "SymbolName.h"

namespace model::symbols {
const std::string& SymbolName::get_name() const {
    return name_;
}
SymbolName::SymbolName(std::string name) : name_(std::move(name)) {}
}  // namespace symbols
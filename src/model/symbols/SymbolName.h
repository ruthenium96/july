#ifndef SPINNER_SYMBOLNAME_H
#define SPINNER_SYMBOLNAME_H

#include <string>

namespace model::symbols {

struct SymbolName {
  public:
    explicit SymbolName(std::string);
    SymbolName() = default;
    const std::string& get_name() const;
    std::strong_ordering operator<=>(const SymbolName& rhs) const = default;

  private:
    std::string name_;
};

}  // namespace symbols

#endif  //SPINNER_SYMBOLNAME_H

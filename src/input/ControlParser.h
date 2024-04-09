#ifndef SPINNER_CONTROLPARSER_H
#define SPINNER_CONTROLPARSER_H

#include <optional>
#include <yaml-cpp/yaml.h>

#include "src/common/Logger.h"
#include "src/entities/data_structures/FactoriesList.h"

namespace input {

class ControlParser {
  public:
    explicit ControlParser(YAML::Node control_node);
    common::PrintLevel getPrintLevel() const;
    const std::optional<quantum::linear_algebra::FactoriesList>& getFactoriesList() const;
  private:
    void constructFactoriesList(YAML::Node& control_node);

    std::optional<common::PrintLevel> print_level_;
    std::optional<quantum::linear_algebra::FactoriesList> factoriesList_;
};

}  // namespace input

#endif  //SPINNER_CONTROLPARSER_H

#ifndef SPINNER_OPTIMIZATIONSPARSER_H
#define SPINNER_OPTIMIZATIONSPARSER_H

#include <optional>
#include <yaml-cpp/yaml.h>

#include "src/common/physical_optimization/OptimizationList.h"

namespace input {

class OptimizationsParser {
  public:
    explicit OptimizationsParser(YAML::Node optimizations_node);
    const std::optional<common::physical_optimization::OptimizationList>&
    getOptimizationList() const;
  private:
    std::optional<common::physical_optimization::OptimizationList> optimizations_list_;
    void customParser(YAML::Node custom_node);
    void symmetrizerParser(YAML::Node symmetrizer_node);
    void groupParser(YAML::Node group_node);
};

}  // namespace input

#endif  //SPINNER_OPTIMIZATIONSPARSER_H

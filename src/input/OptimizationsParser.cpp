#include "OptimizationsParser.h"
#include "Tools.h"

namespace input {
OptimizationsParser::OptimizationsParser(YAML::Node optimizations_node) {
    auto mode_string = extractValue<std::string>(optimizations_node, "mode");
    optimizations_list_ = common::physical_optimization::OptimizationList();
    if (mode_string == "none") {
        // do nothing
    } else if (mode_string == "custom") {
        customParser(extractValue<YAML::Node>(optimizations_node, "custom"));
    } else if (mode_string == "auto") {
        throw std::invalid_argument("optimizations::mode == auto is not implemented yet");
    } else {
        throw std::invalid_argument("Incorrect argument of optimizations::mode: " + mode_string);
    }
    throw_if_node_is_not_empty(optimizations_node);
}

const std::optional<common::physical_optimization::OptimizationList>&
OptimizationsParser::getOptimizationList() const {
    return optimizations_list_;
}

void OptimizationsParser::customParser(YAML::Node custom_node) {
    symmetrizerParser(extractValue<YAML::Node>(custom_node, "symmetrizer"));

    if (extractValue<YAML::Node>(custom_node, "tz_sorter").IsDefined()) {
        optimizations_list_->TzSort();
    }
    if (extractValue<YAML::Node>(custom_node, "positive_tz_eliminator").IsDefined()) {
        optimizations_list_->EliminatePositiveProjections();
    }
    if (extractValue<YAML::Node>(custom_node, "s2_transformer").IsDefined()) {
        optimizations_list_->SSquaredTransform();
    }

    throw_if_node_is_not_empty(custom_node);
}

void OptimizationsParser::symmetrizerParser(YAML::Node symmetrizer_node) {
    if (!symmetrizer_node.IsDefined()) {
        return;
    }
    // if symmetrizer_node is sequence, it is already correct on its level
    if (!symmetrizer_node.IsSequence()) {
        throw std::invalid_argument("Incorrect format of optimizations::symmetrizer");
    }
    for (auto group_node : symmetrizer_node) {
        groupParser(group_node);
    }
}

void OptimizationsParser::groupParser(YAML::Node group_node) {
    auto group_name = extractValue<group::Group::GroupTypeEnum>(group_node, "group_name");
    auto generators_vector = extractValue<std::vector<group::Permutation>>(group_node, "generators");

    throw_if_node_is_not_empty(group_node);

    auto group = group::Group(group_name, generators_vector);
    optimizations_list_->Symmetrize(group);
}
}  // namespace input
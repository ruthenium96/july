#ifndef SPINNER_TOOLS_H
#define SPINNER_TOOLS_H

#include <magic_enum.hpp>
#include <sstream>
#include <string>
#include <yaml-cpp/yaml.h>

namespace input {

template<class T,
         typename std::enable_if<!std::is_enum<T>{}, bool>::type = true>
T extractValue(YAML::Node& node, const std::string& key_name) {
    if (!node.IsDefined()) {
        throw std::invalid_argument("Trying to extract from an empty node");
    }
    auto value = node[key_name].as<T>();
    node.remove(key_name);
    return value;
}

template<class E,
         typename std::enable_if<std::is_enum<E>{}, bool>::type = true>
E extractValue(YAML::Node& node, const std::string& key_name) {
    if (!node.IsDefined()) {
        throw std::invalid_argument("Trying to extract from an empty node");
    }
    auto enum_value_string = node[key_name].as<std::string>();
    node.remove(key_name);

    auto mb_enum_value = magic_enum::enum_cast<E>(enum_value_string, magic_enum::case_insensitive);
    if (!mb_enum_value.has_value()) {
        std::string error_message = "Incorrect enum value: " + enum_value_string +
            "\nPossible values:";
        for (const auto& possible_enum_value : magic_enum::enum_names<E>()) {
            error_message += "\n";
            error_message += possible_enum_value;
        }
        throw std::invalid_argument(error_message);
    } else {
        return mb_enum_value.value();
    }
}

std::vector<double> range_as(YAML::Node& node);

void throw_if_node_is_not_empty(const YAML::Node& node);

}

#endif  //SPINNER_TOOLS_H

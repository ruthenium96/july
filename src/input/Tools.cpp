#include "Tools.h"

namespace input {
std::vector<double> range_as(YAML::Node& node) {
    std::vector<double> answer;
    if (node.Tag() != "!range") {
        answer = node.as<std::vector<double>>();
    } else {
        auto range = node.as<std::array<double, 3>>();
        if (range[1] - range[0] < 0) {
            throw std::invalid_argument("Invalid order of !range parameters");
        }
        if (range[2] < 0) {
            throw std::invalid_argument("Negative sign of step in !range");
        }
        for (double t = range[0]; t < range[1]; t += range[2]) {
            answer.push_back(t);
        }
    }
    return answer;
}

void throw_if_node_is_not_empty(const YAML::Node& node) {
    if (node.size() > 0) {
        std::stringstream stream;

        stream << node << std::endl;

        throw std::invalid_argument("Unknown keys were found:\n" + stream.str());
    }
}
}
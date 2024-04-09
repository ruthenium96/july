#include "ModelInputParser.h"

#include <functional>

#include "Tools.h"

template<typename T, typename U, typename W>
std::vector<W> zip_and_map(const std::vector<T>& ts,
                   const std::variant<U, std::vector<U>>& us,
                   std::function<W(const T&, const U&)> f) {
    std::vector<W> answer;
    if (std::holds_alternative<U>(us)) {
        answer.reserve(ts.size());
        for (const auto& t : ts) {
            auto w = f(t, std::get<U>(us));
            answer.emplace_back(std::move(w));
        }
        return answer;
    } else {
        auto& us_v = std::get<std::vector<U>>(us);
        if (ts.empty() || us_v.empty()) {
            throw std::invalid_argument("Cannot zip empty vectors");
        }
        if (ts.size() == 1) {
            answer.reserve(us_v.size());
            for (const auto& u : us_v) {
                auto w = f(ts[0], u);
                answer.emplace_back(std::move(w));
            }
            return answer;
        }
        if (ts.size() != us_v.size()) {
            throw std::invalid_argument("Cannot zip two vectors of different sizes");
        }
        for (size_t i = 0; i < ts.size(); ++i) {
            auto w = f(ts[i], us_v[i]);
            answer.emplace_back(std::move(w));
        }
        return answer;
    }
}

template<typename T, typename U, typename W>
std::vector<W> cartesian_and_map(const std::vector<T>& ts,
                                 const std::vector<U>& us,
                                 std::function<W(const T&, const U&)> f) {
    std::vector<W> answer;
    if (ts.empty() || us.empty()) {
        throw std::invalid_argument("Cannot cartesian-product empty vectors");
    }

    for (const auto& t : ts) {
        for (const auto& u : us) {
            auto w = f(t, u);
            answer.emplace_back(std::move(w));
        }
    }
    return answer;
}

namespace input {
ModelInputParser::ModelInputParser(YAML::Node model_input_node) {
    auto mults =
        extractValue<std::vector<spin_algebra::Multiplicity>>(model_input_node, "multiplicities");
    model_input_.emplace_back(mults);

    mode_.emplace(extractValue<ModelInputModeEnum>(model_input_node, "mode"));

    parametersParser(extractValue<YAML::Node>(model_input_node, "parameters"));

    throw_if_node_is_not_empty(model_input_node);
}

const std::vector<model::ModelInput>& ModelInputParser::getModelInputs() const {
    return model_input_;
}

void ModelInputParser::parametersParser(YAML::Node parameters_node) {
    // if parameters_node is sequence, it is already correct on its level
    if (!parameters_node.IsDefined()) {
        throw std::invalid_argument("Cannot find model::parameters block");
    }
    if (!parameters_node.IsSequence()) {
        throw std::invalid_argument("Incorrect format of model_input::parameters");
    }
    for (auto parameter_node : parameters_node) {
        parameterParser(parameter_node);
    }
}

void ModelInputParser::parameterParser(YAML::Node parameter_node) {
    auto symbol_name_string = extractValue<std::string>(parameter_node, "name");
    auto is_symbol_fixed = extractValue<bool>(parameter_node, "fixed");
    auto symbol_type = extractValue<
        model::symbols::SymbolTypeEnum
        >(parameter_node, "type");

    std::variant<double, std::vector<double>> symbol_values =
        valuesParser(extractValue<YAML::Node>(parameter_node, "value"));

    std::function<model::ModelInput(const model::ModelInput&, const double&)> function =
        [&symbol_name_string, is_symbol_fixed, symbol_type, &parameter_node]
        (const model::ModelInput& old_model_input, const double& value){
            return returnModifiedModelInput(parameter_node,
                                            old_model_input, value,
                                            symbol_name_string, is_symbol_fixed, symbol_type);
        };

    std::vector<model::ModelInput> answer;

    if (mode_.value() == Trajectory) {
        answer = zip_and_map(model_input_, symbol_values, function);
    }
    if (mode_.value() == Scan) {
        answer = cartesian_and_map(model_input_, std::get<std::vector<double>>(symbol_values), function);
    }
    if (mode_.value() == Single) {
        answer.emplace_back(function(model_input_[0], std::get<double>(symbol_values)));
    }

    model_input_ = std::move(answer);

    if (symbol_type == model::symbols::g_factor || symbol_type == model::symbols::D) {
        parameter_node.remove("centers");
    }
    if (symbol_type == model::symbols::J) {
        parameter_node.remove("pairs");
    }

    throw_if_node_is_not_empty(parameter_node);
}

std::variant<double, std::vector<double>> ModelInputParser::valuesParser(YAML::Node values_node) {
    std::variant<double, std::vector<double>> answer;

    if (mode_.value() == Single) {
        answer = values_node.as<double>();
    } else if (mode_.value() == Scan) {
        answer.emplace<std::vector<double>>();
        auto& vector = std::get<std::vector<double>>(answer);

        if (values_node.IsSequence()) {
            vector = range_as(values_node);
        } else if (values_node.IsScalar()) {
            vector.emplace_back(values_node.as<double>());
        } else {
            throw std::invalid_argument("Incorrect format of model_input::parameters::value");
        }
    } else if (mode_.value() == Trajectory) {
        if (values_node.IsSequence()) {
            answer.emplace<std::vector<double>>();
            auto& vector = std::get<std::vector<double>>(answer);
            vector = range_as(values_node);
            if (trajectory_size_.has_value() && trajectory_size_.value() != vector.size()) {
                throw std::invalid_argument("Different sizes of trajectory sequences");
            } else {
                trajectory_size_ = vector.size();
            }
        } else if (values_node.IsScalar()) {
            answer = values_node.as<double>();
        } else {
            throw std::invalid_argument("Incorrect format of model_input::parameters::value");
        }
    }
    return answer;
}

model::ModelInput ModelInputParser::returnModifiedModelInput(
    YAML::Node& parameter_node,
    const model::ModelInput& old_model_input,
    const double& value,
    const std::string& symbol_name_string,
    bool is_symbol_fixed,
    model::symbols::SymbolTypeEnum symbol_type) {
    auto new_model_input = old_model_input;
    auto symbol_name = new_model_input.addSymbol(
        symbol_name_string,
        value,
        !is_symbol_fixed,
        symbol_type);

    if (symbol_type == model::symbols::Theta) {
        new_model_input.assignSymbolToTheta(symbol_name);
    } else if (symbol_type == model::symbols::g_factor) {
        auto center_node = parameter_node["centers"];
        if (center_node.IsSequence()) {
            auto centers_vector = center_node.as<std::vector<size_t>>();
            for (auto center : centers_vector) {
                new_model_input.assignSymbolToGFactor(symbol_name, center);
            }
        } else {
            auto center_text = center_node.as<std::string>();
            if (center_text == "all") {
                for (size_t center = 0; center < new_model_input.getMults().size(); ++center) {
                    new_model_input.assignSymbolToGFactor(symbol_name, center);
                }
            } else {
                throw std::invalid_argument("Incorrect value of model::parameters::centers: " + center_text);
            }
        }
    } else if (symbol_type == model::symbols::D) {
        auto centers_vector = parameter_node["centers"].as<std::vector<size_t>>();
        for (auto center : centers_vector) {
            new_model_input.assignSymbolToZFSNoAnisotropy(symbol_name, center);
        }
    } else if (symbol_type == model::symbols::J) {
        auto pairs_vector = parameter_node["pairs"]
                                .as<std::vector<std::pair<size_t, size_t>>>();
        for (auto pair : pairs_vector) {
            new_model_input.assignSymbolToIsotropicExchange(
                symbol_name, pair.first, pair.second);
        }
    }
    return new_model_input;
}
}  // namespace input
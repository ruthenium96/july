#ifndef SPINNER_MODELINPUTPARSER_H
#define SPINNER_MODELINPUTPARSER_H

#include <yaml-cpp/yaml.h>

#include <variant>

#include "src/model/ModelInput.h"

namespace input {

class ModelInputParser {
  public:
    explicit ModelInputParser(YAML::Node model_input_node);
    const std::vector<model::ModelInput>& getModelInputs() const;

  private:
    enum ModelInputModeEnum { Single, Trajectory, Scan };

    std::optional<ModelInputModeEnum> mode_;

    std::vector<model::ModelInput> model_input_;
    void parametersParser(YAML::Node parameters_node);
    void parameterParser(YAML::Node symbol_node);
    std::variant<double, std::vector<double>> valuesParser(YAML::Node values_node);

    std::optional<size_t> trajectory_size_;

    static model::ModelInput returnModifiedModelInput(
        YAML::Node& parameter_node,
        const model::ModelInput& old_model_input,
        const double& value,
        const std::string& symbol_name_string,
        bool is_symbol_fixed,
        model::symbols::SymbolTypeEnum symbol_type);
};

}  // namespace input

#endif  //SPINNER_MODELINPUTPARSER_H

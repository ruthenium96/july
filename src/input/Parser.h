#ifndef SPINNER_PARSER_H
#define SPINNER_PARSER_H

#include "src/input/ControlParser.h"
#include "src/input/JobParser.h"
#include "src/input/ModelInputParser.h"
#include "src/input/OptimizationsParser.h"
#include "src/common/physical_optimization/OptimizationList.h"

#include <string>
#include <optional>

namespace input {

class Parser {
  public:
    explicit Parser(const std::string& filename);

    const std::vector<model::ModelInput>& getModelInputs() const;
    const std::optional<common::physical_optimization::OptimizationList>&
    getOptimizationList() const;
    const std::optional<std::vector<double>>& getTemperaturesForSimulation() const;
    const std::optional<std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver>>& getNonlinearSolver() const;
    const std::optional<std::shared_ptr<magnetic_susceptibility::ExperimentalValuesWorker>>& getExperimentalValuesWorker() const;

    const std::optional<quantum::linear_algebra::FactoriesList>& getFactoriesList() const;
  private:
    std::optional<input::ModelInputParser> model_input_parser_;
    std::optional<input::OptimizationsParser> optimizations_list_parser_;
    std::optional<input::JobParser> job_parser_;
    std::optional<input::ControlParser> control_parser_;
};

}  // namespace input

#endif  //SPINNER_PARSER_H

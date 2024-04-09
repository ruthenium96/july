#ifndef SPINNER_JOBPARSER_H
#define SPINNER_JOBPARSER_H

#include <optional>
#include <yaml-cpp/yaml.h>

#include "src/entities/magnetic_susceptibility/assistant/ExperimentalValuesWorker.h"
#include "src/nonlinear_solver/AbstractNonlinearSolver.h"

namespace input {

class JobParser {
  public:
    explicit JobParser(YAML::Node job_node);
    const std::optional<std::vector<double>>& getTemperaturesForSimulation() const;
    const std::optional<std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver>>& getNonlinearSolver() const;
    const std::optional<std::shared_ptr<magnetic_susceptibility::ExperimentalValuesWorker>>& getExperimentalValuesWorker() const;
  private:
    std::optional<std::vector<double>> temperatures_for_simulation_;
    void simulationParser(YAML::Node simulation_node);

    std::optional<std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver>> nonlinear_solver_;
    std::optional<std::shared_ptr<magnetic_susceptibility::ExperimentalValuesWorker>> experimental_values_worker_;
    void fitParser(YAML::Node fit_node);
    void experimentParser(YAML::Node experiment_node);

    static void read_data_from_file(std::vector<magnetic_susceptibility::ValueAtTemperature>& exp_data,
                             const std::string& exp_filename);
};

}  // namespace input

#endif  //SPINNER_JOBPARSER_H

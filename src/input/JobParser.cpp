#include "JobParser.h"

#include <fstream>
#include <iostream>

#include "Tools.h"
#include "src/nonlinear_solver/optimNM/optimNMAdapter.h"
#include "src/nonlinear_solver/stlbfgs/stlbfgsAdapter.h"

#ifdef _LBFGSppAdapter_BUILT
    #include "src/nonlinear_solver/LBFGSpp/LBFGSppAdapter.h"
#endif


template<>
magnetic_susceptibility::ValueAtTemperature YAML::Node::as() const {
    auto temp_and_value = as<std::array<double, 2>>();

    return {temp_and_value[0], temp_and_value[1]};
}

template<>
std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver> YAML::Node::as() const {
    auto solver_string = as<std::string>();

    if (solver_string == "optim_nm") {
        return std::make_shared<nonlinear_solver::optimNMAdapter>();
    } else if (solver_string == "stl_bfgs") {
        return std::make_shared<nonlinear_solver::stlbfgsAdapter>();
    } else if (solver_string == "lbfgs_cpp") {
#ifdef _LBFGSppAdapter_BUILT
        return std::make_shared<nonlinear_solver::LBFGSppAdapter>();
#else
        throw std::invalid_argument("LBFGSppAdapter was not built");
#endif
    } else {
        throw std::invalid_argument("Incorrect job::fit::solver: " + solver_string);
    }
}

namespace input {
const std::optional<std::vector<double>>& JobParser::getTemperaturesForSimulation() const {
    return temperatures_for_simulation_;
}

JobParser::JobParser(YAML::Node job_node) {
    auto mode_string = extractValue<std::string>(job_node, "mode");

    if (mode_string == "simulation") {
        simulationParser(extractValue<YAML::Node>(job_node, "simulation"));
    } else if (mode_string == "fit") {
        fitParser(extractValue<YAML::Node>(job_node, "fit"));
    } else {
        throw std::invalid_argument("Incorrect argument of job::mode: " + mode_string);
    }

    throw_if_node_is_not_empty(job_node);
}

void JobParser::simulationParser(YAML::Node simulation_node) {
    auto temperatures_node = extractValue<YAML::Node>(simulation_node, "temperatures");

    temperatures_for_simulation_ = range_as(temperatures_node);

    throw_if_node_is_not_empty(simulation_node);
}

void JobParser::fitParser(YAML::Node fit_node) {
    nonlinear_solver_.emplace(
        extractValue<std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver>>(
            fit_node, "solver"
            )
        );

    experimentParser(extractValue<YAML::Node>(fit_node, "experiment"));

    throw_if_node_is_not_empty(fit_node);
}

const std::optional<std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver>>&
JobParser::getNonlinearSolver() const {
    return nonlinear_solver_;
}

void JobParser::experimentParser(YAML::Node experiment_node) {
    auto ratio = extractValue<double>(experiment_node, "ratio");
    auto dimension =
        extractValue<magnetic_susceptibility::ExperimentalValuesEnum>(experiment_node,
                                                                      "dimension");
    auto weighting_scheme =
        extractValue<magnetic_susceptibility::WeightingSchemeEnum>(experiment_node,
                                                                   "weights");
    auto data_node = extractValue<YAML::Node>(experiment_node, "data");

    std::vector<magnetic_susceptibility::ValueAtTemperature> exp_data;

    if (data_node.Tag() != "!file") {
        exp_data = data_node.as<std::vector<magnetic_susceptibility::ValueAtTemperature>>();
    } else {
        auto exp_filename = data_node.as<std::string>();
        read_data_from_file(exp_data, exp_filename);
    }

    experimental_values_worker_ =
        std::make_shared<magnetic_susceptibility::ExperimentalValuesWorker>(
            exp_data,
            dimension,
            ratio,
            weighting_scheme);

    throw_if_node_is_not_empty(experiment_node);
}

const std::optional<std::shared_ptr<magnetic_susceptibility::ExperimentalValuesWorker>>&
JobParser::getExperimentalValuesWorker() const {
    return experimental_values_worker_;
}

void JobParser::read_data_from_file(
    std::vector<magnetic_susceptibility::ValueAtTemperature>& exp_data,
    const std::string& exp_filename) {
    std::ifstream exp_filestream (exp_filename);
    while (exp_filestream) {
        double temp;
        double value;
        exp_filestream >> temp >> value;
        magnetic_susceptibility::ValueAtTemperature temp_and_value = {temp, value};
        exp_data.emplace_back(temp_and_value);
    }
}
}  // namespace input
#ifndef SPINNER_EXPERIMENTALVALUESWORKER_H
#define SPINNER_EXPERIMENTALVALUESWORKER_H

#include <vector>

#include "WeightingScheme.h"

namespace magnetic_susceptibility {

constexpr double mu_squared_in_bohr_magnetons_squared_to_chiT_in_cm_cubed_kelvin_per_mol =
    0.125048612;

struct ValueAtTemperature {
    double temperature;
    double value;
};

enum ExperimentalValuesEnum {
    mu_in_bohr_magnetons,
    mu_squared_in_bohr_magnetons_squared,
    chiT_in_cm_cubed_kelvin_per_mol,
    chi_in_cm_cubed_per_mol
};

// Stores experimental data and calculates all values, depending on it.
// Namely, residual error (\sum_T (exp(T) - theor(T))^2)
// and its chain rule derivative: 2 * (exp(T) - theor(T)).
class ExperimentalValuesWorker {
  public:
    ExperimentalValuesWorker(
        std::vector<ValueAtTemperature> experimental_values,
        ExperimentalValuesEnum experimental_values_type,
        double number_of_centers_ratio,
        WeightingSchemeEnum weightingSchemeEnum = per_point);

    std::vector<double> getTemperatures() const;
    void setTheoreticalMuSquared(std::vector<ValueAtTemperature> theoretical_mu_squared);

    std::vector<ValueAtTemperature> getTheoreticalValues() const;
    std::vector<ValueAtTemperature> getExperimentalMuSquared() const;
    const std::vector<double>& getWeights() const;

    double calculateResidualError() const;
    std::vector<ValueAtTemperature> calculateDerivative() const;

  private:
    void throw_exception_if_theoretical_mu_squared_was_not_initialized() const;

    std::vector<ValueAtTemperature> experimental_mu_squared_;
    std::vector<ValueAtTemperature> theoretical_mu_squared_;
    ExperimentalValuesEnum experimental_values_type_;
    WeightingScheme weights_;
    double weighted_sum_of_squares_of_experiment_;
    // Ratio between number of centers (or molar weight) of theoretical and experimental systems.
    double number_of_centers_ratio_;
};
}  // namespace magnetic_susceptibility
#endif  //SPINNER_EXPERIMENTALVALUESWORKER_H

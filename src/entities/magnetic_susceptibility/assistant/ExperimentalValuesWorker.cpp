#include "ExperimentalValuesWorker.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "src/common/PrintingFunctions.h"

namespace magnetic_susceptibility {

ExperimentalValuesWorker::ExperimentalValuesWorker(
    std::vector<ValueAtTemperature> experimental_values,
    ExperimentalValuesEnum experimental_values_type,
    double number_of_centers_ratio,
    WeightingSchemeEnum weightingSchemeEnum) {
    if (experimental_values.empty()) {
        throw std::length_error("experimental_values cannot be empty");
    }

    experimental_values_type_ = experimental_values_type;
    number_of_centers_ratio_ = number_of_centers_ratio;

    common::initialExperimentalValuesPrint(experimental_values, experimental_values_type);

    if (experimental_values_type == mu_squared_in_bohr_magnetons_squared) {
        // DO NOTHING
    } else if (experimental_values_type == mu_in_bohr_magnetons) {
        for (ValueAtTemperature& v : experimental_values) {
            v.value = v.value * v.value;
        }
    } else if (experimental_values_type == chiT_in_cm_cubed_kelvin_per_mol) {
        for (ValueAtTemperature& v : experimental_values) {
            v.value /= mu_squared_in_bohr_magnetons_squared_to_chiT_in_cm_cubed_kelvin_per_mol;
        }
    } else if (experimental_values_type == chi_in_cm_cubed_per_mol) {
        for (ValueAtTemperature& v : experimental_values) {
            v.value = v.value * v.temperature
                / mu_squared_in_bohr_magnetons_squared_to_chiT_in_cm_cubed_kelvin_per_mol;
        }
    }

    for (ValueAtTemperature& v : experimental_values) {
        v.value *= number_of_centers_ratio_;
    }

    std::sort(experimental_values.begin(), experimental_values.end(), [](auto lhs, auto rhs) {
        return lhs.temperature < rhs.temperature;
    });

    experimental_mu_squared_ = std::move(experimental_values);

    weights_ = WeightingScheme(getTemperatures(), weightingSchemeEnum);

    weighted_sum_of_squares_of_experiment_ = 0;
    for (size_t i = 0; i < experimental_mu_squared_.size(); ++i) {
        double value = experimental_mu_squared_[i].value;
        weighted_sum_of_squares_of_experiment_ += weights_.at(i) * value * value;
    }

    common::experimentalValuesPrint(
        getExperimentalMuSquared(),
        getWeights());
}

double ExperimentalValuesWorker::calculateResidualError() const {
    throw_exception_if_theoretical_mu_squared_was_not_initialized();
    double residual_error = 0;
    for (size_t i = 0; i < theoretical_mu_squared_.size(); ++i) {
        double diff = theoretical_mu_squared_[i].value - experimental_mu_squared_[i].value;
        residual_error += weights_.at(i) * diff * diff;
    }
    residual_error /= weighted_sum_of_squares_of_experiment_;
    return residual_error;
}

std::vector<ValueAtTemperature> ExperimentalValuesWorker::calculateDerivative() const {
    throw_exception_if_theoretical_mu_squared_was_not_initialized();
    std::vector<ValueAtTemperature> derivative(experimental_mu_squared_.size());
    for (size_t i = 0; i < theoretical_mu_squared_.size(); ++i) {
        double diff = theoretical_mu_squared_[i].value - experimental_mu_squared_[i].value;
        derivative[i].temperature = theoretical_mu_squared_[i].temperature;
        derivative[i].value = 2 * diff * weights_.at(i) / weighted_sum_of_squares_of_experiment_;
    }
    return derivative;
}

std::vector<double> ExperimentalValuesWorker::getTemperatures() const {
    std::vector<double> temperatures(experimental_mu_squared_.size());
    for (size_t i = 0; i < temperatures.size(); ++i) {
        temperatures[i] = experimental_mu_squared_[i].temperature;
    }
    return temperatures;
}

void ExperimentalValuesWorker::setTheoreticalMuSquared(
    std::vector<ValueAtTemperature> theoretical_mu_squared) {
    if (theoretical_mu_squared.size() != experimental_mu_squared_.size()) {
        throw std::length_error("theoretical_mu_squared.size() != experimental_mu_squared_.size()");
    }
    theoretical_mu_squared_ = std::move(theoretical_mu_squared);
}

std::vector<ValueAtTemperature> ExperimentalValuesWorker::getTheoreticalValues() const {
    throw_exception_if_theoretical_mu_squared_was_not_initialized();

    std::vector<ValueAtTemperature> theoretical_values = theoretical_mu_squared_;

    for (ValueAtTemperature& v : theoretical_values) {
        v.value /= number_of_centers_ratio_;
    }

    if (experimental_values_type_ == mu_squared_in_bohr_magnetons_squared) {
        // DO NOTHING
    } else if (experimental_values_type_ == mu_in_bohr_magnetons) {
        for (ValueAtTemperature& v : theoretical_values) {
            v.value = sqrt(v.value);
        }
    } else if (experimental_values_type_ == chiT_in_cm_cubed_kelvin_per_mol) {
        for (ValueAtTemperature& v : theoretical_values) {
            v.value *= mu_squared_in_bohr_magnetons_squared_to_chiT_in_cm_cubed_kelvin_per_mol;
        }
    } else if (experimental_values_type_ == chi_in_cm_cubed_per_mol) {
        for (ValueAtTemperature& v : theoretical_values) {
            v.value = v.value / v.temperature
                * mu_squared_in_bohr_magnetons_squared_to_chiT_in_cm_cubed_kelvin_per_mol;
        }
    }
    return theoretical_values;
}

std::vector<ValueAtTemperature> ExperimentalValuesWorker::getExperimentalMuSquared() const {
    return experimental_mu_squared_;
}

void ExperimentalValuesWorker::throw_exception_if_theoretical_mu_squared_was_not_initialized()
    const {
    if (theoretical_mu_squared_.empty()) {
        throw std::length_error("Theoretical mu-squared was not initialized");
    }
}

const std::vector<double>& ExperimentalValuesWorker::getWeights() const {
    return weights_.getWeights();
}

}  // namespace magnetic_susceptibility
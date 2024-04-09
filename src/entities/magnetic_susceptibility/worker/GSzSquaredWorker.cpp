#include "GSzSquaredWorker.h"

namespace magnetic_susceptibility::worker {
GSzSquaredWorker::GSzSquaredWorker(
    std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>&& energy,
    std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>&& degeneracy,
    std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>&& g_sz_squared) :
    BasicWorker(std::move(energy), std::move(degeneracy)),
    g_sz_squared_(std::move(g_sz_squared)) {}

double GSzSquaredWorker::calculateTheoreticalMuSquared(double temperature) const {
    double g_sz_squared_averaged =
        3 * ensemble_averager_.ensemble_average(g_sz_squared_, temperature);
    return g_sz_squared_averaged;
}

std::vector<ValueAtTemperature> GSzSquaredWorker::calculateDerivative(
    model::symbols::SymbolTypeEnum symbol_type,
    std::map<common::QuantityEnum, std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>
        values_derivatives_map) const {
    std::vector<double> temperatures = experimental_values_worker_.value()->getTemperatures();
    std::vector<ValueAtTemperature> derivatives(temperatures.size());
    // let A = (\sum_i g_i S_{iz})^2
    if (symbol_type == model::symbols::SymbolTypeEnum::g_factor) {
        // if energy does not depend on g factors:
        // d(mu_squared)/dg = d<A>/dg = <dA/dg>
        const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>&
            g_sz_squared_derivative = values_derivatives_map[common::gSz_total_squared];
        for (size_t i = 0; i < temperatures.size(); ++i) {
            double value =
                3 * ensemble_averager_.ensemble_average(g_sz_squared_derivative, temperatures[i]);
            derivatives[i] = {temperatures[i], value};
        }
    } else {
        // d(mu_squared)/da = d(<A>)/da = (<A>*<dE/da>-<A*dE/da>)/T
        const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>& energy_derivative =
            values_derivatives_map[common::Energy];
        for (size_t i = 0; i < temperatures.size(); ++i) {
            double first_term = ensemble_averager_.ensemble_average(g_sz_squared_, temperatures[i])
                * ensemble_averager_.ensemble_average(energy_derivative, temperatures[i]);
            double second_term = ensemble_averager_.ensemble_average(
                g_sz_squared_->element_wise_multiplication(energy_derivative),
                temperatures[i]);
            double value = 3 * (first_term - second_term) / temperatures[i];
            derivatives[i] = {temperatures[i], value};
        }
    }
    return derivatives;
}

}  // namespace magnetic_susceptibility::worker
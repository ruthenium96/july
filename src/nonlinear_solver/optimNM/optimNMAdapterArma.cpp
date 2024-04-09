#include "optimNMAdapter.h"

#include <armadillo>
#include "src/common/Logger.h"

// Do not use OpenMP in solver: it leads to race condition in Runner.
#define OPTIM_DONT_USE_OPENMP
#define OPTIM_ENABLE_ARMA_WRAPPERS
#include <optim.hpp>

namespace {
std::vector<double> convertFromArmaToSTL(const arma::vec& arma_vector) {
    return arma::conv_to<std::vector<double>>::from(arma_vector);
}

arma::vec convertFromSTLToArma(const std::vector<double>& std_vector) {
    return arma::conv_to<arma::vec>::from(std_vector);
}

std::function<double(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)>
adaptSignature(const std::function<double(const std::vector<double>&, std::vector<double>&, bool)>&
                   oneStepFunction) {
    // do nothing with opt_data
    auto adaptedSignatureFunction = [oneStepFunction](
                                        const arma::vec& changeable_values_arma,
                                        arma::vec* gradient_arma,
                                        void* opt_data) {
        std::vector<double> changeable_values_stl = convertFromArmaToSTL(changeable_values_arma);
        // we do not need gradients on NM method:
        std::vector<double> gradient_stl(changeable_values_stl.size(), 0);
        // the last boolean is doesGradientsRequired
        double residual_error = oneStepFunction(changeable_values_stl, gradient_stl, false);
        // do nothing with gradient_stl
        return residual_error;
    };
    return adaptedSignatureFunction;
}
}  // namespace

namespace nonlinear_solver {

void optimNMAdapter::optimize(
    std::function<double(const std::vector<double>&, std::vector<double>&, bool)> oneStepFunction,
    std::vector<double>& changeable_values) {
    auto adaptedSignatureFunction = adaptSignature(oneStepFunction);
    auto changeable_values_arma = convertFromSTLToArma(changeable_values);

    optim::algo_settings_t algo_settings;
    algo_settings.rel_objfn_change_tol = 1E-10;
    algo_settings.rel_sol_change_tol = 1E-16;
    if (common::Logger::get_level() <= common::PrintLevel::debug) {
        algo_settings.print_level = 3;
    }

    auto spdlog_buffer = common::SpdlogStreamMsg(common::debug);
    auto coutbuf = std::cout.rdbuf(&spdlog_buffer);
    optim::nm(changeable_values_arma, adaptedSignatureFunction, nullptr, algo_settings);
    std::cout.rdbuf(coutbuf);

    changeable_values = convertFromArmaToSTL(changeable_values_arma);
}
}  // namespace nonlinear_solver
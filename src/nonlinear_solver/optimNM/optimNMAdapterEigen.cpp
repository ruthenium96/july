#include "optimNMAdapter.h"

#include <Eigen/Core>
#include "src/common/Logger.h"

// Do not use OpenMP in solver: it leads to race condition in Runner.
#define OPTIM_DONT_USE_OPENMP
#define OPTIM_ENABLE_EIGEN_WRAPPERS
#include <optim.hpp>

namespace {
std::vector<double> convertFromEigenToSTL(const Eigen::VectorXd& eigen_vector) {
    std::vector<double> std_vector;
    std_vector.resize(eigen_vector.size());
    Eigen::VectorXd::Map(&std_vector[0], eigen_vector.size()) = eigen_vector;
    return std_vector;
}
Eigen::VectorXd convertFromSTLToEigen(const std::vector<double>& std_vector) {
    Eigen::VectorXd b =
        Eigen::Map<const Eigen::VectorXd, Eigen::Unaligned>(std_vector.data(), std_vector.size());
    return b;
}

std::function<double(const Eigen::VectorXd& vals_inp, Eigen::VectorXd* grad_out, void* opt_data)>
adaptSignature(const std::function<double(const std::vector<double>&, std::vector<double>&, bool)>&
                   oneStepFunction) {
    // do nothing with opt_data
    auto adaptedSignatureFunction = [oneStepFunction](
                                        const Eigen::VectorXd& changeable_values_arma,
                                        Eigen::VectorXd* gradient_arma,
                                        void* opt_data) {
        std::vector<double> changeable_values_stl = convertFromEigenToSTL(changeable_values_arma);
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
    auto changeable_values_eigen = convertFromSTLToEigen(changeable_values);

    optim::algo_settings_t algo_settings;
    algo_settings.rel_objfn_change_tol = 1E-10;
    algo_settings.rel_sol_change_tol = 1E-16;
    if (common::Logger::get_level() <= common::PrintLevel::debug) {
        algo_settings.print_level = 3;
    }

    auto spdlog_buffer = common::SpdlogStreamMsg(common::debug);
    auto coutbuf = std::cout.rdbuf(&spdlog_buffer);
    optim::nm(changeable_values_eigen, adaptedSignatureFunction, nullptr, algo_settings);
    std::cout.rdbuf(coutbuf);

    changeable_values = convertFromEigenToSTL(changeable_values_eigen);
}
}  // namespace nonlinear_solver
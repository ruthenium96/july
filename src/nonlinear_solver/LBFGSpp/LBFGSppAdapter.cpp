#include "LBFGSppAdapter.h"

#include <LBFGS.h>

#include <Eigen/Core>

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
std::function<double(const Eigen::VectorXd&, Eigen::VectorXd&)>
adaptSignature(const std::function<double(const std::vector<double>&, std::vector<double>&, bool)>&
                   oneStepFunction) {
    auto adaptedSignatureFunction = [oneStepFunction](
                                        const Eigen::VectorXd& changeable_values_eigen,
                                        Eigen::VectorXd& gradient_eigen) {
        std::vector<double> changeable_values_stl = convertFromEigenToSTL(changeable_values_eigen);
        std::vector<double> gradient_stl = convertFromEigenToSTL(gradient_eigen);
        // the last boolean is doesGradientsRequired
        double residual_error = oneStepFunction(changeable_values_stl, gradient_stl, true);
        gradient_eigen = convertFromSTLToEigen(gradient_stl);
        return residual_error;
    };
    return adaptedSignatureFunction;
}
}  // namespace

namespace nonlinear_solver {
void LBFGSppAdapter::optimize(
    std::function<double(const std::vector<double>&, std::vector<double>&, bool)> oneStepFunction,
    std::vector<double>& changeable_values) {
    auto adaptedSignatureFunction = adaptSignature(oneStepFunction);
    auto changeable_values_eigen = convertFromSTLToEigen(changeable_values);

    LBFGSpp::LBFGSParam<double> param;
    param.epsilon = 1e-14;
    param.epsilon_rel = 1e-14;
    param.max_iterations = 1000;

    LBFGSpp::LBFGSSolver<double> solver(param);
    double fx;
    solver.minimize(adaptedSignatureFunction, changeable_values_eigen, fx);

    changeable_values = convertFromEigenToSTL(changeable_values_eigen);
}
}  // namespace nonlinear_solver
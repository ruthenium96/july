#include "stlbfgsAdapter.h"

#include <stlbfgs.h>

namespace {
std::function<void(const std::vector<double>&, double&, std::vector<double>&)>
adaptSignature(const std::function<double(const std::vector<double>&, std::vector<double>&, bool)>&
                   oneStepFunction) {
    auto adaptedSignatureFunction = [oneStepFunction](
                                        const std::vector<double>& changeable_values,
                                        double& residual_error,
                                        std::vector<double>& gradient) {
        // the last boolean is doesGradientsRequired
        residual_error = oneStepFunction(changeable_values, gradient, true);
    };
    return adaptedSignatureFunction;
}
}  // namespace

namespace nonlinear_solver {
void stlbfgsAdapter::optimize(
    std::function<double(const std::vector<double>&, std::vector<double>&, bool)> oneStepFunction,
    std::vector<double>& changeable_values) {
    auto adaptedSignatureFunction = adaptSignature(oneStepFunction);
    STLBFGS::Optimizer optimizer {adaptedSignatureFunction};
    optimizer.verbose = false;
    optimizer.ftol = 1e-16;
    optimizer.gtol = 1e-18;
    // Run calculation from initial guess. STLBFGS updates changeable_values every iteration.
    optimizer.run(changeable_values);
}
}  // namespace nonlinear_solver
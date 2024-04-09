#ifndef SPINNER_UNITABSTRACTSOLVER_TESTS_H
#define SPINNER_UNITABSTRACTSOLVER_TESTS_H

#include <random>

#include "gtest/gtest.h"
#include "src/nonlinear_solver/AbstractNonlinearSolver.h"
#include "tests/tools/concreteSolverConstructors/createConcreteSolver.h"

#define RESIDUAL_ERROR_EPSILON 1e-3

template<class T>
class find_local_minima: public testing::Test {
  protected:
    find_local_minima() : solver_(createConcreteSolver<T>()) {}
    std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver> const solver_;
};

TYPED_TEST_SUITE_P(find_local_minima);

TYPED_TEST_P(find_local_minima, threeD_paraboloid) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> value_dist(-10, +10);

    for (size_t _ = 0; _ < 20; ++_) {
        const double a = value_dist(rng);
        const double b = value_dist(rng);

        std::function<double(const std::vector<double>&, std::vector<double>&, bool)>
            oneStepFunction = [a,
                               b](const std::vector<double>& changeable_values,
                                  std::vector<double>& gradient,
                                  bool isGradientRequired) {
                // z = (x - a)^2 + (y - b)^2
                if (isGradientRequired) {
                    gradient.resize(2);
                    gradient[0] = 2 * (changeable_values[0] - a);
                    gradient[1] = 2 * (changeable_values[1] - b);
                }
                return (changeable_values[0] - a) * (changeable_values[0] - a)
                    + (changeable_values[1] - b) * (changeable_values[1] - b);
            };

        std::vector<double> values = {0, 0};

        this->solver_->optimize(oneStepFunction, values);
        EXPECT_NEAR(values[0], a, RESIDUAL_ERROR_EPSILON);
        EXPECT_NEAR(values[1], b, RESIDUAL_ERROR_EPSILON);
    }
}

TYPED_TEST_P(find_local_minima, fourD_paraboloid) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> value_dist(-10, +10);

    for (size_t _ = 0; _ < 20; ++_) {
        const double a = value_dist(rng);
        const double b = value_dist(rng);
        const double c = value_dist(rng);

        std::function<double(const std::vector<double>&, std::vector<double>&, bool)>
            oneStepFunction = [a, b, c](
                                  const std::vector<double>& changeable_values,
                                  std::vector<double>& gradient,
                                  bool isGradientRequired) {
                // w = (x - a)^2 + (y - b)^2 + (z - c)^2
                if (isGradientRequired) {
                    gradient.resize(3);
                    gradient[0] = 2 * (changeable_values[0] - a);
                    gradient[1] = 2 * (changeable_values[1] - b);
                    gradient[2] = 2 * (changeable_values[2] - c);
                }
                return (changeable_values[0] - a) * (changeable_values[0] - a)
                    + (changeable_values[1] - b) * (changeable_values[1] - b)
                    + (changeable_values[2] - c) * (changeable_values[2] - c);
            };

        std::vector<double> values = {0, 0, 0};

        this->solver_->optimize(oneStepFunction, values);
        EXPECT_NEAR(values[0], a, RESIDUAL_ERROR_EPSILON);
        EXPECT_NEAR(values[1], b, RESIDUAL_ERROR_EPSILON);
        EXPECT_NEAR(values[2], c, RESIDUAL_ERROR_EPSILON);
    }
}

inline double smooth_heaviside(double x) {
    if (std::abs(x) < 1e-6) {
        return x / 1e-6;
    } else if (x > 0) {
        return 1.0;
    } else {
        return -1.0;
    }
}

TYPED_TEST_P(find_local_minima, twoD_absolute) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> value_dist(-5, +5);

    for (size_t _ = 0; _ < 20; ++_) {
        const double a = value_dist(rng);

        std::function<double(const std::vector<double>&, std::vector<double>&, bool)>
            oneStepFunction = [a](const std::vector<double>& changeable_values,
                                  std::vector<double>& gradient,
                                  bool isGradientRequired) {
                // y = |x - a|
                if (isGradientRequired) {
                    gradient.resize(1);
                    gradient[0] = smooth_heaviside(changeable_values[0] - a);
                }
                return std::abs(changeable_values[0] - a);
            };

        std::vector<double> values = {0};

        this->solver_->optimize(oneStepFunction, values);
        EXPECT_NEAR(values[0], a, RESIDUAL_ERROR_EPSILON);
    }
}

TYPED_TEST_P(find_local_minima, threeD_trigonometric) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> value_dist(-0.5, +0.5);

    for (size_t _ = 0; _ < 20; ++_) {
        const double a = value_dist(rng);
        const double b = value_dist(rng);

        std::function<double(const std::vector<double>&, std::vector<double>&, bool)>
            oneStepFunction = [a,
                               b](const std::vector<double>& changeable_values,
                                  std::vector<double>& gradient,
                                  bool isGradientRequired) {
                // z = sin((x - a)^2 + (y - b)^2)
                double sum_of_squares = (changeable_values[0] - a) * (changeable_values[0] - a)
                    + (changeable_values[1] - b) * (changeable_values[1] - b);
                if (isGradientRequired) {
                    gradient.resize(2);
                    gradient[0] = 2 * (changeable_values[0] - a) * cos(sum_of_squares);
                    gradient[1] = 2 * (changeable_values[1] - b) * cos(sum_of_squares);
                }
                return sin(sum_of_squares);
            };

        std::vector<double> values = {0, 0};

        this->solver_->optimize(oneStepFunction, values);
        EXPECT_NEAR(values[0], a, RESIDUAL_ERROR_EPSILON);
        EXPECT_NEAR(values[1], b, RESIDUAL_ERROR_EPSILON);
    }
}

REGISTER_TYPED_TEST_SUITE_P(
    find_local_minima,
    threeD_paraboloid,
    fourD_paraboloid,
    twoD_absolute,
    threeD_trigonometric);


#endif  //SPINNER_UNITABSTRACTSOLVER_TESTS_H

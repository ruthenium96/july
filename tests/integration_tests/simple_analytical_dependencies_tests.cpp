#include <random>

#include "gtest/gtest.h"
#include "src/common/runner/Runner.h"

double sum_of_s_squared(const lexicographic::IndexConverter& indexConverter) {
    double sum_of_s_squared = 0;
    for (auto spin : indexConverter.get_spins()) {
        sum_of_s_squared += spin * (spin + 1);
    }
    return sum_of_s_squared;
}

// mu_2 = g^2 \sum_a s_a (s_a + 1)

TEST(simple_analytical_dependencies, nothing) {
    std::vector<std::vector<spin_algebra::Multiplicity>> mults_cases =
        {{2}, {3}, {4}, {5}, {2, 2}, {2, 2, 2, 2}, {2, 3, 4, 5}};

    for (const auto& mults : mults_cases) {
        const double g_exact = 2.0;

        std::vector<magnetic_susceptibility::ValueAtTemperature> values;

        {
            model::ModelInput model(mults);
            double g_value = g_exact;
            auto g = model.addSymbol("g", g_value);
            for (size_t i = 0; i < mults.size(); ++i) {
                model.assignSymbolToGFactor(g, i);
            }

            runner::Runner runner(model);

            double sum_of_s_squared_ = sum_of_s_squared(runner.getIndexConverter());

            for (size_t temperature = 1; temperature < 301; ++temperature) {
                double calculated_value =
                    runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(
                        temperature);
                double exact_value = (g_exact * g_exact * sum_of_s_squared_);
                // TODO: epsilon
                EXPECT_NEAR(calculated_value, exact_value, 1e-9);
            }
        }
    }
}

// mu_2 = T/(T-Theta) g^2 \sum_a s_a (s_a + 1)

TEST(simple_analytical_dependencies, Theta) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> Theta_dist(-100, -1.0);

    std::vector<std::vector<spin_algebra::Multiplicity>> mults_cases =
        {{2}, {3}, {4}, {5}, {2, 2}, {2, 2, 2, 2}, {2, 3, 4, 5}};

    for (const auto& mults : mults_cases) {
        for (size_t _ = 0; _ < 10; ++_) {
            const double Theta_exact = Theta_dist(rng);
            const double g_exact = 2.0;

            std::vector<magnetic_susceptibility::ValueAtTemperature> values;

            {
                model::ModelInput model(mults);
                auto Theta = model.addSymbol("Theta", Theta_exact);
                model.assignSymbolToTheta(Theta);
                double g_value = g_exact;
                auto g = model.addSymbol("g", g_value);
                for (size_t i = 0; i < mults.size(); ++i) {
                    model.assignSymbolToGFactor(g, i);
                }

                runner::Runner runner(model);

                double sum_of_s_squared_ = sum_of_s_squared(runner.getIndexConverter());

                for (size_t temperature = 1; temperature < 301; ++temperature) {
                    double calculated_value =
                        runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(
                            temperature);
                    double exact_value = (g_exact * g_exact * sum_of_s_squared_) * temperature
                        / (temperature - Theta_exact);
                    // TODO: epsilon
                    EXPECT_NEAR(calculated_value, exact_value, 1e-9);
                }
            }
        }
    }
}

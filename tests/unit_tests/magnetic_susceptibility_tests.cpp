#include <random>
#include <utility>

#include "gtest/gtest.h"
#include "src/common/runner/Runner.h"

TEST(magnetic_susceptibility, throw_experimental_values_worker_set_values_length_error) {
    std::vector<magnetic_susceptibility::ValueAtTemperature> exp_values = {{1, 1}, {2, 2}, {3, 3}};

    magnetic_susceptibility::ExperimentalValuesWorker experimental_values_worker(
        exp_values,
        magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
        1);

    std::vector<magnetic_susceptibility::ValueAtTemperature> theor_values = exp_values;
    theor_values.pop_back();

    EXPECT_THROW(
        experimental_values_worker.setTheoreticalMuSquared(theor_values),
        std::length_error);
}

TEST(magnetic_susceptibility, throw_experimental_values_worker_get_values_empty) {
    std::vector<magnetic_susceptibility::ValueAtTemperature> exp_values = {{1, 1}, {2, 2}, {3, 3}};

    magnetic_susceptibility::ExperimentalValuesWorker experimental_values_worker(
        exp_values,
        magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
        1);

    EXPECT_THROW(experimental_values_worker.getTheoreticalValues(), std::length_error);
}

TEST(magnetic_susceptibility, throw_experimental_values_worker_residual_error_empty) {
    std::vector<magnetic_susceptibility::ValueAtTemperature> exp_values = {{1, 1}, {2, 2}, {3, 3}};

    magnetic_susceptibility::ExperimentalValuesWorker experimental_values_worker(
        exp_values,
        magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
        1);

    EXPECT_THROW(experimental_values_worker.calculateResidualError(), std::length_error);
}

TEST(magnetic_susceptibility, throw_experimental_values_worker_derivative_empty) {
    std::vector<magnetic_susceptibility::ValueAtTemperature> exp_values = {{1, 1}, {2, 2}, {3, 3}};

    magnetic_susceptibility::ExperimentalValuesWorker experimental_values_worker(
        exp_values,
        magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
        1);

    EXPECT_THROW(experimental_values_worker.calculateDerivative(), std::length_error);
}

TEST(magnetic_susceptibility, throw_experimental_values_worker_empty_experimental_values) {
    std::vector<magnetic_susceptibility::ValueAtTemperature> exp_values = {};

    EXPECT_THROW(
        magnetic_susceptibility::ExperimentalValuesWorker experimental_values_worker(
            exp_values,
            magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
            1),
        std::length_error);
}

TEST(magnetic_susceptibility, do_not_throw_experimental_before_theoretical) {
    std::vector<magnetic_susceptibility::ValueAtTemperature> exp_values = {{1, 20}};
    std::vector<spin_algebra::Multiplicity> mults = {2, 2};
    model::ModelInput model(mults);
    double J_value = 10;
    auto J = model.addSymbol("J", J_value);
    model.assignSymbolToIsotropicExchange(J, 0, 1);
    double g_value = 2.0;
    auto g = model.addSymbol("g", g_value);
    for (size_t i = 0; i < mults.size(); ++i) {
        model.assignSymbolToGFactor(g, i);
    }

    runner::Runner runner(model);

    runner.initializeExperimentalValues(
        exp_values,
        magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
        1);

    EXPECT_NO_THROW(runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(100));
}

TEST(magnetic_susceptibility, do_not_throw_theoretical_before_experimental) {
    std::vector<magnetic_susceptibility::ValueAtTemperature> exp_values = {{1, 20}};
    std::vector<spin_algebra::Multiplicity> mults = {2, 2};
    model::ModelInput model(mults);
    double J_value = 10;
    auto J = model.addSymbol("J", J_value);
    model.assignSymbolToIsotropicExchange(J, 0, 1);
    double g_value = 2.0;
    auto g = model.addSymbol("g", g_value);
    for (size_t i = 0; i < mults.size(); ++i) {
        model.assignSymbolToGFactor(g, i);
    }

    runner::Runner runner(model);

    runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(100);
    EXPECT_NO_THROW(runner.initializeExperimentalValues(
        exp_values,
        magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
        1));
}

TEST(magnetic_susceptibility, sort_experimental_temperatues) {
    std::vector<magnetic_susceptibility::ValueAtTemperature> exp_values =
        {{10, 20}, {4, 20}, {5, 20}, {11, 20}, {8, 20}, {9, 20}};
    auto worker = magnetic_susceptibility::ExperimentalValuesWorker(
        exp_values,
        magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
        1);
    auto temperatures = worker.getTemperatures();
    EXPECT_TRUE(std::is_sorted(temperatures.begin(), temperatures.end()));
}

TEST(magnetic_susceptibility, equvalence_of_weightong_schemes_on_equidistant_data) {
    std::vector<magnetic_susceptibility::ValueAtTemperature> exp_values;

    for (size_t i = 1; i < 300; i += 3) {
        magnetic_susceptibility::ValueAtTemperature value_at_temperature = {
            static_cast<double>(i),
            1.0};
        exp_values.push_back(value_at_temperature);
    }

    magnetic_susceptibility::ExperimentalValuesWorker worker_per_point(
        exp_values,
        magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
        1.0,
        magnetic_susceptibility::per_point);
    magnetic_susceptibility::ExperimentalValuesWorker worker_per_interval(
        exp_values,
        magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
        1.0,
        magnetic_susceptibility::per_interval);

    const auto& weights_per_point = worker_per_point.getWeights();
    const auto& weights_per_interval = worker_per_interval.getWeights();

    EXPECT_EQ(weights_per_interval.size(), weights_per_point.size());
    for (size_t i = 0; i < weights_per_interval.size(); ++i) {
        EXPECT_NEAR(weights_per_interval.at(i), weights_per_point.at(i), 1e-9);
    }
}

TEST(magnetic_susceptibility, value_mu_squared_reversibility) {
    std::vector<magnetic_susceptibility::ExperimentalValuesEnum> values_enum = {
        magnetic_susceptibility::mu_in_bohr_magnetons,
        magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
        magnetic_susceptibility::chiT_in_cm_cubed_kelvin_per_mol,
        magnetic_susceptibility::chi_in_cm_cubed_per_mol};

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(1, 100);

    std::vector<magnetic_susceptibility::ValueAtTemperature> exp_values;

    for (size_t i = 1; i < 300; ++i) {
        magnetic_susceptibility::ValueAtTemperature value_at_temperature = {
            static_cast<double>(i),
            static_cast<double>(dist(rng))};
        exp_values.push_back(value_at_temperature);
    }

    for (const auto& value_enum : values_enum) {
        for (size_t ratio = 1; ratio < 10; ++ratio) {
            magnetic_susceptibility::ExperimentalValuesWorker worker(exp_values, value_enum, ratio);
            worker.setTheoreticalMuSquared(worker.getExperimentalMuSquared());
            auto double_transformed = worker.getTheoreticalValues();
            for (size_t i = 0; i < double_transformed.size(); ++i) {
                EXPECT_DOUBLE_EQ(double_transformed[i].value, exp_values[i].value);
            }
        }
    }
}

TEST(magnetic_susceptibility, unique_g_different_g_difference) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> g_dist(1.8, 3.0);

    for (size_t _ = 0; _ < 20; ++_) {
        const double g_factor = g_dist(rng);
        std::vector<spin_algebra::Multiplicity> mults = {3, 3, 3, 3};

        std::vector<magnetic_susceptibility::ValueAtTemperature> values_unique;

        {
            model::ModelInput model(mults);
            auto g = model.addSymbol("g", g_factor);
            for (size_t i = 0; i < mults.size(); ++i) {
                model.assignSymbolToGFactor(g, i);
            }

            runner::Runner runner(model);

            for (size_t i = 1; i < 301; ++i) {
                magnetic_susceptibility::ValueAtTemperature value_at_temperature = {
                    static_cast<double>(i),
                    runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(i)};
                values_unique.push_back(value_at_temperature);
            }
        }

        std::vector<magnetic_susceptibility::ValueAtTemperature> values_different;
        {
            model::ModelInput model(mults);
            auto g_one = model.addSymbol("g1", g_factor);
            auto g_two = model.addSymbol("g2", g_factor);
            for (size_t i = 0; i < mults.size() / 2; ++i) {
                model.assignSymbolToGFactor(g_one, i);
            }
            for (size_t i = mults.size() / 2; i < mults.size(); ++i) {
                model.assignSymbolToGFactor(g_two, i);
            }

            runner::Runner runner(model);

            for (size_t i = 1; i < 301; ++i) {
                magnetic_susceptibility::ValueAtTemperature value_at_temperature = {
                    static_cast<double>(i),
                    runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(i)};
                values_different.push_back(value_at_temperature);
            }
        }

        EXPECT_EQ(values_unique.size(), values_different.size());
        for (size_t i = 0; i < values_unique.size(); ++i) {
            const auto& [temperature_unique, value_unique] = values_unique[i];
            const auto& [temperature_different, value_different] = values_different[i];
            EXPECT_EQ(temperature_unique, temperature_different);
            EXPECT_NEAR(value_unique, value_different, 1e-9);
        }
    }
}

TEST(magnetic_susceptibility, unique_g_different_g_difference_J) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> J_dist(-100, -10);
    std::uniform_real_distribution<double> g_dist(1.8, 3.0);

    for (size_t _ = 0; _ < 10; ++_) {
        const double J_exact = J_dist(rng);
        const double g_factor = g_dist(rng);

        std::vector<spin_algebra::Multiplicity> mults = {3, 3, 3, 3};

        std::vector<magnetic_susceptibility::ValueAtTemperature> values_unique;

        {
            model::ModelInput model(mults);
            auto J = model.addSymbol("J", J_exact);
            model.assignSymbolToIsotropicExchange(J, 0, 1)
                .assignSymbolToIsotropicExchange(J, 1, 2)
                .assignSymbolToIsotropicExchange(J, 2, 3)
                .assignSymbolToIsotropicExchange(J, 3, 0);

            auto g = model.addSymbol("g", g_factor);
            for (size_t i = 0; i < mults.size(); ++i) {
                model.assignSymbolToGFactor(g, i);
            }

            runner::Runner runner(model);

            for (size_t i = 1; i < 301; ++i) {
                magnetic_susceptibility::ValueAtTemperature value_at_temperature = {
                    static_cast<double>(i),
                    runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(i)};
                values_unique.push_back(value_at_temperature);
            }
        }

        std::vector<magnetic_susceptibility::ValueAtTemperature> values_different;
        {
            model::ModelInput model(mults);
            auto J = model.addSymbol("J", J_exact);
            model.assignSymbolToIsotropicExchange(J, 0, 1)
                .assignSymbolToIsotropicExchange(J, 1, 2)
                .assignSymbolToIsotropicExchange(J, 2, 3)
                .assignSymbolToIsotropicExchange(J, 3, 0);
            auto g_one = model.addSymbol("g1", g_factor);
            auto g_two = model.addSymbol("g2", g_factor);
            for (size_t i = 0; i < mults.size() / 2; ++i) {
                model.assignSymbolToGFactor(g_one, i);
            }
            for (size_t i = mults.size() / 2; i < mults.size(); ++i) {
                model.assignSymbolToGFactor(g_two, i);
            }

            runner::Runner runner(model);

            for (size_t i = 1; i < 301; ++i) {
                magnetic_susceptibility::ValueAtTemperature value_at_temperature = {
                    static_cast<double>(i),
                    runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(i)};
                values_different.push_back(value_at_temperature);
            }
        }

        EXPECT_EQ(values_unique.size(), values_different.size());
        for (size_t i = 0; i < values_unique.size(); ++i) {
            const auto& [temperature_unique, value_unique] = values_unique[i];
            const auto& [temperature_different, value_different] = values_different[i];
            EXPECT_EQ(temperature_unique, temperature_different);
            EXPECT_NEAR(value_unique, value_different, 1e-9);
        }
    }
}

model::ModelInput constructFourCenterModel_g_J(
    const std::vector<spin_algebra::Multiplicity>& mults,
    double J_value,
    double g_value) {
    model::ModelInput model(mults);
    auto J = model.addSymbol("J", J_value);
    model.assignSymbolToIsotropicExchange(J, 0, 1)
        .assignSymbolToIsotropicExchange(J, 1, 2)
        .assignSymbolToIsotropicExchange(J, 2, 3)
        .assignSymbolToIsotropicExchange(J, 3, 0);

    auto g = model.addSymbol("g", g_value);
    for (size_t i = 0; i < mults.size(); ++i) {
        model.assignSymbolToGFactor(g, i);
    }
    return model;
}

double calculateResidualError(
    model::ModelInput model,
    const std::vector<magnetic_susceptibility::ValueAtTemperature>& values) {
    runner::Runner runner(std::move(model));

    runner.initializeExperimentalValues(
        values,
        magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
        1);

    return runner.getMagneticSusceptibilityController().calculateResidualError();
}

TEST(magnetic_susceptibility, analytical_derivative_vs_finite_differences_J_g) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> J_dist(-100, -10);
    std::uniform_real_distribution<double> g_dist(1.8, 3.0);

    for (size_t _ = 0; _ < 20; ++_) {
        const double J_exact = J_dist(rng);
        const double g_factor = g_dist(rng);

        std::vector<spin_algebra::Multiplicity> mults = {3, 3, 3, 3};

        std::vector<magnetic_susceptibility::ValueAtTemperature> values;
        {
            auto model = constructFourCenterModel_g_J(mults, J_exact, g_factor);
            runner::Runner runner(model);

            for (size_t i = 1; i < 301; ++i) {
                magnetic_susceptibility::ValueAtTemperature value_at_temperature = {
                    static_cast<double>(i),
                    runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(i)};
                values.push_back(value_at_temperature);
            }
        }

        double J_value = -50;
        double g_value = 2.0;

        double analytical_dR2_wrt_dJ;
        double analytical_dR2_wrt_dg;

        double delta_J = 1e-5;
        double delta_g = 1e-8;
        double finite_difference_dR2_wrt_dJ;
        double finite_difference_dR2_wrt_dg;

        {
            auto model = constructFourCenterModel_g_J(mults, J_value, g_value);

            runner::Runner runner(model);

            runner.initializeExperimentalValues(
                values,
                magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
                1);

            auto J = runner.getSymbolicWorker().getChangeableNames(model::symbols::J)[0];
            auto g = runner.getSymbolicWorker().getChangeableNames(model::symbols::g_factor)[0];

            auto derivative_map = runner.calculateTotalDerivatives();
            analytical_dR2_wrt_dJ = derivative_map[J];
            analytical_dR2_wrt_dg = derivative_map[g];
        }

        {
            auto model_initial = constructFourCenterModel_g_J(mults, J_value, g_value);
            double initial_R2 = calculateResidualError(model_initial, values);

            auto model_delta_J = constructFourCenterModel_g_J(mults, J_value + delta_J, g_value);
            double delta_J_R2 = calculateResidualError(model_delta_J, values);
            finite_difference_dR2_wrt_dJ = (delta_J_R2 - initial_R2) / delta_J;

            auto model_delta_g = constructFourCenterModel_g_J(mults, J_value, g_value + delta_g);
            double delta_g_R2 = calculateResidualError(model_delta_g, values);
            finite_difference_dR2_wrt_dg = (delta_g_R2 - initial_R2) / delta_g;
        }

        double dR2_wrt_dJ_range = std::abs(analytical_dR2_wrt_dJ / 1000);
        double dR2_wrt_dg_range = std::abs(analytical_dR2_wrt_dg / 1000);

        EXPECT_NEAR(analytical_dR2_wrt_dJ, finite_difference_dR2_wrt_dJ, dR2_wrt_dJ_range)
            << "J_exact = " << J_exact << "\n"
            << "g_exact = " << g_factor << "\n";
        EXPECT_NEAR(analytical_dR2_wrt_dg, finite_difference_dR2_wrt_dg, dR2_wrt_dg_range)
            << "J_exact = " << J_exact << "\n"
            << "g_exact = " << g_factor << "\n";
    }
}

model::ModelInput constructRingModel_J_g_D(
    const std::vector<spin_algebra::Multiplicity>& mults,
    double J_value,
    double g_value,
    double D_value) {
    model::ModelInput model(mults);
    auto J = model.addSymbol("J", J_value);
    size_t size = mults.size();
    for (size_t i = 0; i < size; ++i) {
        model.assignSymbolToIsotropicExchange(J, i, (i + 1) % size);
    }

    auto D = model.addSymbol("D", D_value);
    model.assignSymbolToZFSNoAnisotropy(D, 0);

    auto g = model.addSymbol("g", g_value, false);
    for (size_t i = 0; i < mults.size(); ++i) {
        model.assignSymbolToGFactor(g, i);
    }
    return model;
}

TEST(magnetic_susceptibility, analytical_derivative_vs_finite_differences_J_g_D) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> J_dist(-100, -10);
    std::uniform_real_distribution<double> D_dist(-10, -50);

    for (size_t _ = 0; _ < 20; ++_) {
        const double J_exact = J_dist(rng);
        const double D_exact = D_dist(rng);
        const double g_factor = 2.0;

        std::vector<spin_algebra::Multiplicity> mults = {3, 3, 3, 3};

        std::vector<magnetic_susceptibility::ValueAtTemperature> values;
        {
            auto model = constructRingModel_J_g_D(mults, J_exact, g_factor, D_exact);
            runner::Runner runner(model);

            for (size_t i = 1; i < 301; ++i) {
                magnetic_susceptibility::ValueAtTemperature value_at_temperature = {
                    static_cast<double>(i),
                    runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(i)};
                values.push_back(value_at_temperature);
            }
        }

        double J_value = -50;
        double D_value = 0.1;

        double analytical_dR2_wrt_dJ;
        double analytical_dR2_wrt_dD;

        double delta_J = 1e-5;
        double delta_D = 1e-5;
        double finite_difference_dR2_wrt_dJ;
        double finite_difference_dR2_wrt_dD;

        {
            auto model = constructRingModel_J_g_D(mults, J_value, g_factor, D_value);

            runner::Runner runner(model);

            runner.initializeExperimentalValues(
                values,
                magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
                1);

            auto J = runner.getSymbolicWorker().getChangeableNames(model::symbols::J)[0];
            auto D = runner.getSymbolicWorker().getChangeableNames(model::symbols::D)[0];

            auto derivative_map = runner.calculateTotalDerivatives();
            analytical_dR2_wrt_dJ = derivative_map[J];
            analytical_dR2_wrt_dD = derivative_map[D];
        }

        {
            auto model_initial = constructRingModel_J_g_D(mults, J_value, g_factor, D_value);
            double initial_R2 = calculateResidualError(model_initial, values);

            auto model_delta_J =
                constructRingModel_J_g_D(mults, J_value + delta_J, g_factor, D_value);
            double delta_J_R2 = calculateResidualError(model_delta_J, values);
            finite_difference_dR2_wrt_dJ = (delta_J_R2 - initial_R2) / delta_J;

            auto model_delta_D =
                constructRingModel_J_g_D(mults, J_value, g_factor, D_value + delta_D);
            double delta_D_R2 = calculateResidualError(model_delta_D, values);
            finite_difference_dR2_wrt_dD = (delta_D_R2 - initial_R2) / delta_D;
        }

        double dR2_wrt_dJ_range = std::abs(analytical_dR2_wrt_dJ / 1000);
        double dR2_wrt_dD_range = std::abs(analytical_dR2_wrt_dD / 100);

        EXPECT_NEAR(analytical_dR2_wrt_dJ, finite_difference_dR2_wrt_dJ, dR2_wrt_dJ_range);
        EXPECT_NEAR(analytical_dR2_wrt_dD, finite_difference_dR2_wrt_dD, dR2_wrt_dD_range);
    }
}
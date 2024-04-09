#include "gtest/gtest.h"
#include "src/common/runner/Runner.h"

struct EnergyAndSSquared {
    std::partial_ordering operator<=>(const EnergyAndSSquared&) const = default;
    double energy;
    double s_squared;
};

std::vector<EnergyAndSSquared> construct_final_vector(runner::Runner& runner) {
    std::vector<EnergyAndSSquared> vector;

    auto factory = runner.getDataStructuresFactories();
    auto energy_vector = factory.createVector();
    auto s_squared_vector = factory.createVector();
    auto degeneracy_vector = factory.createVector();

    for (const auto& subspectrum : runner.getSpectrum(common::Energy).blocks) {
        energy_vector->concatenate_with(subspectrum.raw_data);
        degeneracy_vector->add_identical_values(
            subspectrum.raw_data->size(),
            subspectrum.properties.degeneracy);
    }
    energy_vector->subtract_minimum();

    for (const auto& subspectrum : runner.getSpectrum(common::S_total_squared).blocks) {
        s_squared_vector->concatenate_with(subspectrum.raw_data);
    }

    for (size_t i = 0; i < degeneracy_vector->size(); ++i) {
        EnergyAndSSquared energy_and_s_squared = {energy_vector->at(i), s_squared_vector->at(i)};
        for (size_t j = 0; j < (size_t)degeneracy_vector->at(i); ++j) {
            vector.push_back(energy_and_s_squared);
        }
    }
    std::stable_sort(vector.begin(), vector.end());
    return vector;
}

void expect_final_vectors_equivalence(runner::Runner& simple, runner::Runner& second) {

    auto first_vector = construct_final_vector(simple);
    auto second_vector = construct_final_vector(second);

    double s_squared_sum_first = 0;
    double s_squared_sum_second = 0;
    double last_energy = INFINITY;

    for (size_t i = 0; i < first_vector.size(); ++i) {
        // TODO: epsilon
        EXPECT_NEAR(first_vector[i].energy, second_vector[i].energy, 1e-9);
        // TODO: now s2-squared values can be incorrect,
        // EXPECT_NEAR(first_vector[i].s_squared, second_vector[i].s_squared, 1e-9);

        // TODO: but sum of s_squared for degenerate eigenvectors is correct
        if (std::abs(last_energy - first_vector[i].energy) > 1e-5) {
            EXPECT_NEAR(s_squared_sum_first, s_squared_sum_second, 1e-3);
            s_squared_sum_first = 0;
            s_squared_sum_second = 0;
        }
        s_squared_sum_first += first_vector[i].s_squared;
        s_squared_sum_second += second_vector[i].s_squared;
        last_energy = first_vector[i].energy;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void initialize_five_center_mirror_symmetry_exchange_chain(
    model::ModelInput& model,
    double external,
    double internal) {
    auto Jexternal = model.addSymbol("J1", external);
    model.assignSymbolToIsotropicExchange(Jexternal, 0, 1)
        .assignSymbolToIsotropicExchange(Jexternal, 3, 4);

    auto Jinternal = model.addSymbol("J2", internal);
    model.assignSymbolToIsotropicExchange(Jinternal, 1, 2)
        .assignSymbolToIsotropicExchange(Jinternal, 2, 3);
}

TEST(spectrum_final_equivalence, five_center_mirror_symmetry_chain) {
    std::vector<std::vector<spin_algebra::Multiplicity>> multss = {
        {2, 2, 2, 2, 2},
        {2, 2, 3, 2, 2},
        {2, 2, 4, 2, 2},
        {3, 3, 2, 3, 3},
        {3, 3, 3, 3, 3},
        {3, 3, 4, 3, 3}};
    std::vector<std::pair<double, double>> js =
        {{10, 15}, {-10, 15}, {10, -15}, {-10, -15}, {-10, -10}, {10, 10}};

    for (const auto& mults : multss) {
        for (auto [Jfirst, Jsecond] : js) {
            group::Group group(group::Group::S2, {{4, 3, 2, 1, 0}});

            model::ModelInput model(mults);
            initialize_five_center_mirror_symmetry_exchange_chain(model, Jfirst, Jsecond);

            runner::Runner runner_simple(model);

            // TZ_SORTER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort();
                runner::Runner runner_tz_sorted(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner_tz_sorted);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().EliminatePositiveProjections();
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.Symmetrize(group);
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().Symmetrize(group);
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().EliminatePositiveProjections().Symmetrize(group);

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + S2_TRANSFORMER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().EliminatePositiveProjections().SSquaredTransform();

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + SYMMETRIZER + S2_TRANSFORMER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort()
                    .EliminatePositiveProjections()
                    .Symmetrize(group)
                    .SSquaredTransform();

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void initialize_four_centers_exchange_rectangle(
    model::ModelInput& model,
    double first,
    double second) {
    auto Jfirst = model.addSymbol("J1", first);
    model.assignSymbolToIsotropicExchange(Jfirst, 0, 1)
        .assignSymbolToIsotropicExchange(Jfirst, 2, 3);

    auto Jsecond = model.addSymbol("J2", second);
    model.assignSymbolToIsotropicExchange(Jsecond, 1, 2)
        .assignSymbolToIsotropicExchange(Jsecond, 3, 0);
}

TEST(spectrum_final_equivalence, rectangle) {
    std::vector<std::vector<spin_algebra::Multiplicity>> multss = {
        {2, 2, 2, 2},
        {3, 3, 3, 3},
        {4, 4, 4, 4}};
    std::vector<std::pair<double, double>> js =
        {{10, 15}, {-10, 15}, {10, -15}, {-10, -15}, {-10, -10}, {10, 10}};

    for (const auto& mults : multss) {
        for (auto [Jfirst, Jsecond] : js) {
            group::Group first_direction(group::Group::S2, {{1, 0, 3, 2}});
            group::Group second_direction(group::Group::S2, {{3, 2, 1, 0}});

            model::ModelInput model(mults);
            initialize_four_centers_exchange_rectangle(model, Jfirst, Jsecond);

            runner::Runner runner_simple(model);

            // TZ_SORTER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort();
                runner::Runner runner_tz_sorted(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner_tz_sorted);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().EliminatePositiveProjections();
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.Symmetrize(first_direction);
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.Symmetrize(second_direction);
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().Symmetrize(first_direction);
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().Symmetrize(second_direction);
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // SYMMETRIZER + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.Symmetrize(first_direction).Symmetrize(second_direction);
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + SYMMETRIZER + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().Symmetrize(first_direction).Symmetrize(second_direction);

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + SYMMETRIZER + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort()
                    .EliminatePositiveProjections()
                    .Symmetrize(first_direction)
                    .Symmetrize(second_direction);

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + S2_TRANSFORMER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().EliminatePositiveProjections().SSquaredTransform();

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + SYMMETRIZER + S2_TRANSFORMER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort()
                    .EliminatePositiveProjections()
                    .Symmetrize(first_direction)
                    .SSquaredTransform();

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort()
                    .EliminatePositiveProjections()
                    .Symmetrize(second_direction)
                    .SSquaredTransform();

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + SYMMETRIZER + SYMMETRIZER + S2_TRANSFORMER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort()
                    .EliminatePositiveProjections()
                    .Symmetrize(first_direction)
                    .Symmetrize(second_direction)
                    .SSquaredTransform();

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void initialize_three_centers_exchange_triangle(model::ModelInput& model, double first) {
    auto Jfirst = model.addSymbol("J1", first);
    model.assignSymbolToIsotropicExchange(Jfirst, 0, 1)
        .assignSymbolToIsotropicExchange(Jfirst, 1, 2)
        .assignSymbolToIsotropicExchange(Jfirst, 2, 0);
}

TEST(spectrum_final_equivalence, triangle) {
    std::vector<std::vector<spin_algebra::Multiplicity>> multss = {{2, 2, 2}, {3, 3, 3}, {4, 4, 4}};
    std::vector<double> js = {10, 17.17, 33};

    for (const auto& mults : multss) {
        for (auto Jfirst : js) {
            group::Group triangle(group::Group::S3, {{2, 0, 1}, {0, 2, 1}});
            model::ModelInput model(mults);
            initialize_three_centers_exchange_triangle(model, Jfirst);

            runner::Runner runner_simple(model);

            // TZ_SORTER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort();
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().EliminatePositiveProjections();
                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.Symmetrize(triangle);

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            {
                group::Group triangle_diff(group::Group::S3, {{1, 2, 0}, {2, 1, 0}});
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.Symmetrize(triangle_diff);

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().Symmetrize(triangle);

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + SYMMETRIZER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().EliminatePositiveProjections().Symmetrize(triangle);

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
            // SYMMETRIZER + NON_ABELIAN_SIMPLIFIER
            {}  // TZ_SORTER + SYMMETRIZER + NON_ABELIAN_SIMPLIFIER
            {}  // SYMMETRIZER + TZ_SORTER + NON_ABELIAN_SIMPLIFIER
            {}  // SYMMETRIZER + NON_ABELIAN_SIMPLIFIER + TZ_SORTER
            // TZ_SORTER + POSITIVE_PROJECTIONS_ELIMINATOR + S2_TRANSFORMER
            {
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort().EliminatePositiveProjections().SSquaredTransform();

                runner::Runner runner(model, optimizationList);
                expect_final_vectors_equivalence(runner_simple, runner);
            }
        }
    }
}
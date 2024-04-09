#include "gtest/gtest.h"
#include "src/common/runner/ConsistentModelOptimizationList.h"

TEST(consistentModelOptimizationList_tests, throw_groups_with_different_sizes_of_permutations) {
    common::physical_optimization::OptimizationList optimizationList;
    optimizationList.Symmetrize(group::Group::S2, {{0, 2, 1}});
    EXPECT_THROW(optimizationList.Symmetrize(group::Group::S2, {{1, 0}}), std::invalid_argument);
}

TEST(consistentModelOptimizationList_tests, throw_noncommutative_groups) {
    common::physical_optimization::OptimizationList optimizationList;
    optimizationList.Symmetrize(group::Group::S2, {{0, 2, 1}});
    EXPECT_THROW(optimizationList.Symmetrize(group::Group::S2, {{1, 0, 2}}), std::invalid_argument);
}

TEST(consistentModelOptimizationList_tests, throw_wrong_size_of_pemutation) {
    std::vector<spin_algebra::Multiplicity> mults = {4, 4, 4};
    model::ModelInput model(mults);
    common::physical_optimization::OptimizationList optimizationList;
    optimizationList.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});

    EXPECT_THROW(
        runner::ConsistentModelOptimizationList(model, optimizationList),
        std::length_error);
}

TEST(consistentModelOptimizationList_tests, throw_permutes_different_multiplicities) {
    std::vector<spin_algebra::Multiplicity> mults = {4, 4, 4, 3};
    model::ModelInput model(mults);
    common::physical_optimization::OptimizationList optimizationList;
    optimizationList.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});

    EXPECT_THROW(
        runner::ConsistentModelOptimizationList(model, optimizationList),
        std::invalid_argument);
}

TEST(consistentModelOptimizationList_tests, throw_2222_isotropic_inconsistent_symmetry) {
    std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2};
    model::ModelInput model(mults);
    double J_value = 10;
    auto tripledJ = model.addSymbol("3J", 3 * J_value);
    auto J = model.addSymbol("J", J_value);
    auto doubledJ = model.addSymbol("2J", 2 * J_value);
    model
        .assignSymbolToIsotropicExchange(tripledJ, 0, 1)
        .assignSymbolToIsotropicExchange(J, 1, 2)
        .assignSymbolToIsotropicExchange(doubledJ, 2, 3)
        .assignSymbolToIsotropicExchange(J, 3, 0);

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList.TzSort()
        .Symmetrize(group::Group::S2, {{1, 0, 3, 2}})
        .Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

    EXPECT_THROW(
        runner::ConsistentModelOptimizationList(model, optimizationList),
        std::invalid_argument);
}

TEST(consistentModelOptimizationList_tests, throw_2222_isotropic_accidental_symmetry) {
    std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2};
    model::ModelInput model(mults);
    double J = 10;
    auto firstJ = model.addSymbol("J1", J);
    auto secondJ = model.addSymbol("J2", J);
    model.assignSymbolToIsotropicExchange(firstJ, 0, 1)
        .assignSymbolToIsotropicExchange(firstJ, 1, 2)
        .assignSymbolToIsotropicExchange(firstJ, 2, 3)
        .assignSymbolToIsotropicExchange(secondJ, 3, 0);

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList.TzSort()
        .Symmetrize(group::Group::S2, {{1, 0, 3, 2}})
        .Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

    EXPECT_THROW(
        runner::ConsistentModelOptimizationList(model, optimizationList),
        std::invalid_argument);
}

TEST(consistentModelOptimizationList_tests, throw_2222_gfactor_inconsistent_symmetry) {
    std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2};
    model::ModelInput model(mults);
    double J_value = 10;
    auto J = model.addSymbol("J", J_value);
    auto firstg = model.addSymbol("g1", 2.0);
    auto secondg = model.addSymbol("g2", 3.0);
    model.assignSymbolToIsotropicExchange(J, 0, 1)
        .assignSymbolToIsotropicExchange(J, 1, 2)
        .assignSymbolToIsotropicExchange(J, 2, 3)
        .assignSymbolToIsotropicExchange(J, 3, 0)
        .assignSymbolToGFactor(firstg, 0)
        .assignSymbolToGFactor(firstg, 1)
        .assignSymbolToGFactor(firstg, 2)
        .assignSymbolToGFactor(secondg, 3);

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList.TzSort()
        .Symmetrize(group::Group::S2, {{1, 0, 3, 2}})
        .Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

    EXPECT_THROW(
        runner::ConsistentModelOptimizationList(model, optimizationList),
        std::invalid_argument);
}

TEST(consistentModelOptimizationList_tests, throw_2222_gfactor_accidental_symmetry) {
    std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2};
    model::ModelInput model(mults);
    double J_value = 10;
    double g_value = 2.0;
    auto J = model.addSymbol("J", J_value);
    auto firstg = model.addSymbol("g1", g_value);
    auto secondg = model.addSymbol("g2", g_value);
    model.assignSymbolToIsotropicExchange(J, 0, 1)
        .assignSymbolToIsotropicExchange(J, 1, 2)
        .assignSymbolToIsotropicExchange(J, 2, 3)
        .assignSymbolToIsotropicExchange(J, 3, 0)
        .assignSymbolToGFactor(firstg, 0)
        .assignSymbolToGFactor(firstg, 1)
        .assignSymbolToGFactor(firstg, 2)
        .assignSymbolToGFactor(secondg, 3);

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList.TzSort()
        .Symmetrize(group::Group::S2, {{1, 0, 3, 2}})
        .Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

    EXPECT_THROW(
        runner::ConsistentModelOptimizationList(model, optimizationList),
        std::invalid_argument);
}

TEST(consistentModelOptimizationList_tests, throw_SSquaredTransformation_of_ZFS) {
    std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2};
    model::ModelInput model(mults);
    double J_value = 10;
    double g_value = 2.0;
    double D_value = 5.0;
    auto J = model.addSymbol("J", J_value);
    auto g = model.addSymbol("g1", g_value);
    auto D = model.addSymbol("D", D_value);
    model.assignSymbolToIsotropicExchange(J, 0, 1)
        .assignSymbolToIsotropicExchange(J, 1, 2)
        .assignSymbolToIsotropicExchange(J, 2, 3)
        .assignSymbolToIsotropicExchange(J, 3, 0)
        .assignSymbolToGFactor(g, 0)
        .assignSymbolToGFactor(g, 1)
        .assignSymbolToGFactor(g, 2)
        .assignSymbolToGFactor(g, 3)
        .assignSymbolToZFSNoAnisotropy(D, 0);

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList.TzSort().EliminatePositiveProjections().SSquaredTransform();

    EXPECT_THROW(
        runner::ConsistentModelOptimizationList(model, optimizationList),
        std::invalid_argument);
}

TEST(consistentModelOptimizationList_tests, throw_SSquaredTransformation_of_different_gs) {
    std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2};
    model::ModelInput model(mults);
    double J_value = 10;
    double g_value = 2.0;
    double D_value = 5.0;
    auto J = model.addSymbol("J", J_value);
    auto g_one = model.addSymbol("g1", g_value);
    auto g_two = model.addSymbol("g2", g_value);
    model.assignSymbolToIsotropicExchange(J, 0, 1)
        .assignSymbolToIsotropicExchange(J, 1, 2)
        .assignSymbolToIsotropicExchange(J, 2, 3)
        .assignSymbolToIsotropicExchange(J, 3, 0)
        .assignSymbolToGFactor(g_one, 0)
        .assignSymbolToGFactor(g_one, 1)
        .assignSymbolToGFactor(g_two, 2)
        .assignSymbolToGFactor(g_two, 3);

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList.TzSort().EliminatePositiveProjections().SSquaredTransform();

    EXPECT_THROW(
        runner::ConsistentModelOptimizationList(model, optimizationList),
        std::invalid_argument);
}
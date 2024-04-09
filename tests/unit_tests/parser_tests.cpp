#include "gtest/gtest.h"

#include "src/input/JobParser.h"
#include "src/input/ModelInputParser.h"
#include "src/input/OptimizationsParser.h"
#include "src/input/Tools.h"

#include <yaml-cpp/yaml.h>

TEST(parser_tests, extractValuesEnum) {
    enum TestEnum {
        one, TWO, Three, fOUR
    };
    std::vector<std::string> correctValues = {
        "key: one",
        "key: two",
        "key: three",
        "key: four",
    };
    for (const auto& s : correctValues) {
        auto node = YAML::Load(s);
        EXPECT_NO_THROW(input::extractValue<TestEnum>(node, "key"));
        EXPECT_NO_THROW(input::throw_if_node_is_not_empty(node));
    }
    std::string wrongValue = "key: five";
    auto wrong_node = YAML::Load(wrongValue);
    EXPECT_ANY_THROW(input::extractValue<TestEnum>(wrong_node, "key"));
}

TEST(parser_tests, range_as_no_tag) {
    std::vector<std::string> data_to_parse = {
        "[1.0, 2.0, 2.4, 1.5, 3.0]",
        "[1.0, 2]",
        "[100, 123.4, 111]",
    };
    std::vector<std::vector<double>> answers = {
        {1.0, 2.0, 2.4, 1.5, 3.0},
        {1.0, 2},
        {100, 123.4, 111}
    };

    for (size_t i = 0; i < data_to_parse.size(); ++i) {
        auto node = YAML::Load(data_to_parse[i]);
        auto parsed_node = input::range_as(node);
        const auto& answer = answers[i];
        EXPECT_EQ(parsed_node, answer);
    }
}

TEST(parser_tests, range_as_range_tag) {
    std::vector<std::string> data_to_parse = {
        "!range [1.0, 20.0, 0.1]",
        "!range [2.0, 100.1, 10.0]",
        "!range [100, 123.4, 111]",
    };
    std::vector<std::array<double, 3>> answers = {
        {1.0, 20.0, 0.1},
        {2.0, 100.1, 10.0},
        {100, 123.4, 111},
    };

    for (size_t i = 0; i < data_to_parse.size(); ++i) {
        auto node = YAML::Load(data_to_parse[i]);
        auto parsed_node = input::range_as(node);
        const auto& answer = answers[i];

        size_t number_of_element_to_generate = ceil((answer[1] - answer[0]) / answer[2]);
        EXPECT_EQ(number_of_element_to_generate, parsed_node.size());

        double acc = answer[0];
        for (double j : parsed_node) {
            EXPECT_NEAR(j, acc, 1e-9);
            EXPECT_TRUE(acc <= answer[1]);

            acc += answer[2];
        }
    }
}

TEST(parser_tests, optimizations_parser_modes) {
    std::string none_string = R""""(
mode: none
)"""";
    EXPECT_NO_THROW(input::OptimizationsParser(YAML::Load(none_string)));

    std::string empty_custom_string = R""""(
mode: custom
custom:
)"""";
    EXPECT_NO_THROW(input::OptimizationsParser(YAML::Load(empty_custom_string)));
}

TEST(parser_tests, optimizations_parser_throw) {
    {
        std::string string = R""""(
)"""";
        EXPECT_ANY_THROW(input::OptimizationsParser(YAML::Load(string)));
    }
    {
        std::string string = R""""(
mode: _incorrect_value_
)"""";
        EXPECT_ANY_THROW(input::OptimizationsParser(YAML::Load(string)));
    }
    {
        std::string string = R""""(
mode: custom
)"""";
        EXPECT_ANY_THROW(input::OptimizationsParser(YAML::Load(string)));
    }
    {
        std::string string = R""""(
mode: none
custom:
)"""";
        EXPECT_ANY_THROW(input::OptimizationsParser(YAML::Load(string)));
    }
}

TEST(parser_tests, optimizations_parser_custom) {
    {
        std::string string = R""""(
mode: custom
custom:
)"""";
        auto parser = input::OptimizationsParser(YAML::Load(string));

        EXPECT_NO_THROW(parser.getOptimizationList().value());

        const auto& optimizationList = parser.getOptimizationList().value();

        EXPECT_FALSE(optimizationList.isTzSorted());
        EXPECT_FALSE(optimizationList.isPositiveProjectionsEliminated());
        EXPECT_FALSE(optimizationList.isSSquaredTransformed());
        EXPECT_TRUE(optimizationList.getGroupsToApply().empty());
    }
    {
        std::string string = R""""(
mode: custom
custom:
  tz_sorter:
  positive_tz_eliminator:
  s2_transformer:
)"""";
        auto parser = input::OptimizationsParser(YAML::Load(string));

        EXPECT_NO_THROW(parser.getOptimizationList().value());

        const auto& optimizationList = parser.getOptimizationList().value();

        EXPECT_TRUE(optimizationList.isTzSorted());
        EXPECT_TRUE(optimizationList.isPositiveProjectionsEliminated());
        EXPECT_TRUE(optimizationList.isSSquaredTransformed());
        EXPECT_TRUE(optimizationList.getGroupsToApply().empty());
    }
    {
        std::string string = R""""(
mode: custom
custom:
  symmetrizer:
    - group_name: S2
      generators: [[3, 4, 5, 0, 1, 2]]
    - group_name: S3
      generators: [[1, 2, 0, 4, 5, 3], [0, 2, 1, 3, 5, 4]]
)"""";
        auto parser = input::OptimizationsParser(YAML::Load(string));

        EXPECT_NO_THROW(parser.getOptimizationList().value());

        const auto& optimizationList = parser.getOptimizationList().value();

        EXPECT_FALSE(optimizationList.isTzSorted());
        EXPECT_FALSE(optimizationList.isPositiveProjectionsEliminated());
        EXPECT_FALSE(optimizationList.isSSquaredTransformed());
        EXPECT_EQ(optimizationList.getGroupsToApply().size(), 2);
    }
}

TEST(parser_tests, job_parser_modes) {
    {
        std::string string = R""""(
mode: simulation
simulation:
  temperatures: [1.0, 20.0, 200.0]
)"""";
        EXPECT_NO_THROW(input::JobParser(YAML::Load(string)));
    }
    {
        std::string string = R""""(
mode: fit
fit:
  solver: optim_nm
  experiment:
    data: [[1.0, 1.0], [2.0, 3.0], [6.0, 7.0], [11.1, 22.2]]
    dimension: chiT_in_cm_cubed_kelvin_per_mol
    ratio: 1.0
    weights: per_interval
)"""";
        EXPECT_NO_THROW(input::JobParser(YAML::Load(string)));
    }
}

TEST(parser_tests, job_parser_simulation) {
    {
        std::string string = R""""(
mode: simulation
simulation:
  temperatures: [1.0, 20.0, 200.0]
)"""";
        auto parser = input::JobParser(YAML::Load(string));
        EXPECT_NO_THROW(parser.getTemperaturesForSimulation().value());
        EXPECT_FALSE(parser.getNonlinearSolver().has_value());
        EXPECT_FALSE(parser.getExperimentalValuesWorker().has_value());
    }
}

TEST(parser_tests, job_parser_fit) {
    {
        std::string string = R""""(
mode: fit
fit:
  solver: optim_nm
  experiment:
    data: [[1.0, 1.0], [2.0, 3.0], [6.0, 7.0], [11.1, 22.2]]
    dimension: chiT_in_cm_cubed_kelvin_per_mol
    ratio: 1.0
    weights: per_interval
)"""";
        auto parser = input::JobParser(YAML::Load(string));
        EXPECT_FALSE(parser.getTemperaturesForSimulation().has_value());
        EXPECT_NO_THROW(parser.getNonlinearSolver().value());
        EXPECT_NO_THROW(parser.getExperimentalValuesWorker().value());
    }
}

TEST(parser_tests, job_parser_throw) {
    {
        std::string string = R""""(
)"""";
        EXPECT_ANY_THROW(input::JobParser(YAML::Load(string)));
    }
    {
        std::string string = R""""(
mode: _incorrect_value_
)"""";
        EXPECT_ANY_THROW(input::JobParser(YAML::Load(string)));
    }
    {
        std::string string = R""""(
mode: simulation
)"""";
        EXPECT_ANY_THROW(input::JobParser(YAML::Load(string)));
    }
    {
        std::string string = R""""(
mode: simulation
simulation:
)"""";
        EXPECT_ANY_THROW(input::JobParser(YAML::Load(string)));
    }
    {
        std::string string = R""""(
mode: fit
)"""";
        EXPECT_ANY_THROW(input::JobParser(YAML::Load(string)));
    }
    {
        std::string string = R""""(
mode: fit
fit:
)"""";
        EXPECT_ANY_THROW(input::JobParser(YAML::Load(string)));
    }
}

TEST(parser_tests, model_input_parser_simple) {
    {
        std::string string = R""""(
multiplicities: [2, 2, 2, 2]
mode: single
parameters:
  - name: J1
    type: J
    value: -10.
    fixed: true
    pairs: [[0, 1], [2, 3]]
  - name: J2
    type: J
    value: +10.0
    fixed: false
    pairs: [[3, 0], [1, 2]]
  - name: g
    type: g_factor
    value: 2.0
    fixed: false
    centers: [0, 1, 2, 3]
  - name: theta
    type: theta
    value: -3
    fixed: false
)"""";
        auto parser = input::ModelInputParser(YAML::Load(string));

        auto model = model::ModelInput({2, 2, 2, 2});
        auto J1 = model.addSymbol("J1", -10, false);
        auto J2 = model.addSymbol("J2", +10, true);
        auto theta = model.addSymbol("theta", -3.0, true);
        auto g = model.addSymbol("g", 2.0, true, model::symbols::g_factor);
        model.assignSymbolToIsotropicExchange(J1, 0, 1)
            .assignSymbolToIsotropicExchange(J1, 2, 3)
            .assignSymbolToIsotropicExchange(J2, 3, 0)
            .assignSymbolToIsotropicExchange(J2, 1, 2)
            .assignSymbolToGFactor(g, 0)
            .assignSymbolToGFactor(g, 1)
            .assignSymbolToGFactor(g, 2)
            .assignSymbolToGFactor(g, 3)
            .assignSymbolToTheta(theta);

        EXPECT_NO_THROW(parser.getModelInputs());
        EXPECT_EQ(parser.getModelInputs().size(), 1);
        EXPECT_EQ(parser.getModelInputs().at(0), model);
    }
}

TEST(parser_tests, model_input_parser_simple_g_factor_all) {
    {
        std::string string = R""""(
multiplicities: [2, 2, 2, 2]
mode: single
parameters:
  - name: g
    type: g_factor
    value: 2.0
    fixed: false
    centers: all
)"""";
        auto parser = input::ModelInputParser(YAML::Load(string));

        EXPECT_NO_THROW(parser.getModelInputs());
        EXPECT_EQ(parser.getModelInputs().size(), 1);
        EXPECT_TRUE(parser.getModelInputs().at(0).getSymbolicWorker().isAllGFactorsEqual());
    }
}

TEST(parser_tests, model_input_parser_scan) {
    {
        std::string string = R""""(
multiplicities: [2, 2, 2, 2]
mode: scan
parameters:
  - name: J1
    type: J
    value: [-10., -9, -8, -7]
    fixed: true
    pairs: [[0, 1], [2, 3]]
  - name: J2
    type: J
    value: [+10.0, 9, 8]
    fixed: false
    pairs: [[3, 0], [1, 2]]
  - name: g
    type: g_factor
    value: [2.0, 2.1]
    fixed: false
    centers: [0, 1, 2, 3]
  - name: theta
    type: theta
    value: [-3]
    fixed: false
)"""";
        auto parser = input::ModelInputParser(YAML::Load(string));

        EXPECT_NO_THROW(parser.getModelInputs());
        EXPECT_EQ(parser.getModelInputs().size(), 24);
    }
}

TEST(parser_tests, model_input_parser_trajectory) {
    {
        std::string string = R""""(
multiplicities: [2, 2, 2, 2]
mode: trajectory
parameters:
  - name: J1
    type: J
    value: [-10., -9, -8, -7]
    fixed: true
    pairs: [[0, 1], [2, 3]]
  - name: J2
    type: J
    value: [+10.0, 9, 8, 7.0]
    fixed: false
    pairs: [[3, 0], [1, 2]]
  - name: g
    type: g_factor
    value: 2.0
    fixed: false
    centers: [0, 1, 2, 3]
  - name: theta
    type: theta
    value: -3
    fixed: false
)"""";
        auto parser = input::ModelInputParser(YAML::Load(string));

        EXPECT_NO_THROW(parser.getModelInputs());
        EXPECT_EQ(parser.getModelInputs().size(), 4);
    }
    {
        std::string string = R""""(
multiplicities: [2, 2, 2, 2]
mode: trajectory
parameters:
  - name: J1
    type: J
    value: [-10.]
    fixed: true
    pairs: [[0, 1], [2, 3]]
  - name: J2
    type: J
    value: [+10.0]
    fixed: false
    pairs: [[3, 0], [1, 2]]
  - name: g
    type: g_factor
    value: 2.0
    fixed: false
    centers: [0, 1, 2, 3]
  - name: theta
    type: theta
    value: -3
    fixed: false
)"""";
        auto parser = input::ModelInputParser(YAML::Load(string));

        EXPECT_NO_THROW(parser.getModelInputs());
        EXPECT_EQ(parser.getModelInputs().size(), 1);
    }
    {
        std::string string = R""""(
multiplicities: [2, 2, 2, 2]
mode: trajectory
parameters:
  - name: J1
    type: J
    value: -10.
    fixed: true
    pairs: [[0, 1], [2, 3]]
  - name: J2
    type: J
    value: +10.0
    fixed: false
    pairs: [[3, 0], [1, 2]]
  - name: g
    type: g_factor
    value: 2.0
    fixed: false
    centers: [0, 1, 2, 3]
  - name: theta
    type: theta
    value: -3
    fixed: false
)"""";
        auto parser = input::ModelInputParser(YAML::Load(string));

        EXPECT_NO_THROW(parser.getModelInputs());
        EXPECT_EQ(parser.getModelInputs().size(), 1);
    }
}

TEST(parser_tests, model_input_parser_trajectory_throw) {
    {
        std::string string = R""""(
multiplicities: [2, 2, 2, 2]
mode: trajectory
parameters:
  - name: J1
    type: J
    value: [-10., -9, -8, -7]
    fixed: true
    pairs: [[0, 1], [2, 3]]
  - name: J2
    type: J
    value: [+10.0, 9]
    fixed: false
    pairs: [[3, 0], [1, 2]]
)"""";
        EXPECT_ANY_THROW(input::ModelInputParser(YAML::Load(string)));
    }
    {
        std::string string = R""""(
multiplicities: [2, 2, 2, 2]
mode: trajectory
parameters:
  - name: J1
    type: J
    value: [-10., -9, -8, -7]
    fixed: true
    pairs: [[0, 1], [2, 3]]
  - name: J2
    type: J
    value: [+10.0]
    fixed: false
    pairs: [[3, 0], [1, 2]]
)"""";
        EXPECT_ANY_THROW(input::ModelInputParser(YAML::Load(string)));
    }
}

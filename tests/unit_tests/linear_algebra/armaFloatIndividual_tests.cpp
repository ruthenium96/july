#include "abstractIndividual_tests.h"
#include "src/entities/data_structures/arma/ArmaFactories.h"

template<>
std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory>
createConcreteFactory<std::pair<quantum::linear_algebra::ArmaDenseTransformAndDiagonalizeFactory, float>>() {
    auto factory =
        std::make_shared<quantum::linear_algebra::ArmaDenseTransformAndDiagonalizeFactory>();
    factory->setPrecision(quantum::linear_algebra::SINGLE);
    return factory;
};


typedef testing::Types<std::pair<quantum::linear_algebra::ArmaDenseTransformAndDiagonalizeFactory, float>> ArmaFloat;
INSTANTIATE_TYPED_TEST_SUITE_P(
    ArmaFloatIndividualTests,
    AbstractDenseTransformAndDiagonalizeFactoryIndividualTest,
    ArmaFloat);

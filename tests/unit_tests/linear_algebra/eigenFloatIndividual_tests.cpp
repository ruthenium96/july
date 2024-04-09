#include "abstractIndividual_tests.h"
#include "src/entities/data_structures/eigen/EigenFactories.h"

template<>
std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory>
createConcreteFactory<std::pair<quantum::linear_algebra::EigenDenseTransformAndDiagonalizeFactory, float>>() {
    return std::make_shared<quantum::linear_algebra::EigenDenseTransformAndDiagonalizeFactory>();
}

typedef testing::Types<std::pair<quantum::linear_algebra::EigenDenseTransformAndDiagonalizeFactory, float>> EigenFloat;
INSTANTIATE_TYPED_TEST_SUITE_P(
    EigenFloatIndividualTests,
    AbstractDenseTransformAndDiagonalizeFactoryIndividualTest,
    EigenFloat);
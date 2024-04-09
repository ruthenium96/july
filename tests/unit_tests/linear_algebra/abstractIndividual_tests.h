#ifndef SPINNER_ABSTRACT_INDIVIDUAL_TESTS_H
#define SPINNER_ABSTRACT_INDIVIDUAL_TESTS_H

#include "gtest/gtest.h"
#include "src/entities/data_structures/AbstractFactories.h"
#include "tests/tools/GenerateSameMatrix.h"

template<class T>
std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory>
createConcreteFactory();

template<class T>
class AbstractDenseTransformAndDiagonalizeFactoryIndividualTest: public testing::Test {
  protected:
    AbstractDenseTransformAndDiagonalizeFactoryIndividualTest() :
        factory_(createConcreteFactory<T>()) {}
    std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory> const
        factory_;
};

TYPED_TEST_SUITE_P(AbstractDenseTransformAndDiagonalizeFactoryIndividualTest);

// TODO: implement tests for abstract factory, abstract matrix and abstract vector
TYPED_TEST_P(AbstractDenseTransformAndDiagonalizeFactoryIndividualTest, NonNullptrObjects) {
    EXPECT_FALSE(this->factory_->createDenseDiagonalizableMatrix(0) == nullptr);
}

TYPED_TEST_P(
    AbstractDenseTransformAndDiagonalizeFactoryIndividualTest,
    unitary_transformation_and_unitary_transformation_and_return_main_diagonal_equivalence) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    for (size_t size = 2; size < 100; ++size) {
        auto denseDiagonalizableMatrix =
            generateDenseDiagonalizableMatrix(size, this->factory_, dist, rng);
        auto denseUnitaryMatrix = generateDenseUnitaryMatrix(size, this->factory_, dist, rng);

        auto denseDiagonalizableMatrixTransformedMatrix =
            denseUnitaryMatrix->unitaryTransform(denseDiagonalizableMatrix);
        auto denseDiagonalizableMatrixTransformedMainDiagonal =
            denseUnitaryMatrix->unitaryTransformAndReturnMainDiagonal(denseDiagonalizableMatrix);

        // check equality:
        for (size_t i = 0; i < size; ++i) {
            double epsilon =
                std::abs(denseDiagonalizableMatrixTransformedMainDiagonal->at(i) * 1e-6);
            EXPECT_NEAR(
                denseDiagonalizableMatrixTransformedMainDiagonal->at(i),
                denseDiagonalizableMatrixTransformedMatrix->at(i, i),
                epsilon);
        }
    }
}

REGISTER_TYPED_TEST_SUITE_P(
    AbstractDenseTransformAndDiagonalizeFactoryIndividualTest,
    NonNullptrObjects,
    unitary_transformation_and_unitary_transformation_and_return_main_diagonal_equivalence);

#endif  //SPINNER_ABSTRACT_INDIVIDUAL_TESTS_H

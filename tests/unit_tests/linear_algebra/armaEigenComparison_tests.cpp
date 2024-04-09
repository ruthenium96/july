#include <random>

#include "gtest/gtest.h"
#include "tests/tools/AllSymmetricMatrixFactories.h"
#include "tests/tools/GenerateSameMatrix.h"

inline int sign_f(double a, double b)
{
    return -2 * (std::signbit(a) ^ std::signbit(b)) + 1;
}

TEST(linearAlgebraFactories, throw_combination_of_different_objects) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    auto denseUnitaryMatrices = generateDenseUnitaryMatrices(
        10,
        constructAllDenseTransformAndDiagonalizeFactories(),
        dist,
        rng);
    auto denseDiagonalizableMatrices = generateDenseDiagonalizableMatrices(
        10,
        constructAllDenseTransformAndDiagonalizeFactories(),
        dist,
        rng);

    for (int i = 0; i < denseDiagonalizableMatrices.size(); ++i) {
        auto& denseDiagonalizableMatrix_i = denseDiagonalizableMatrices[i];
        for (int j = 0; j < denseUnitaryMatrices.size(); ++j) {
            auto& denseUnitatyMatrix_j = denseUnitaryMatrices[j];
            if (i != j) {
                EXPECT_ANY_THROW(denseUnitatyMatrix_j->
                                 unitaryTransformAndReturnMainDiagonal(denseDiagonalizableMatrix_i));
            } else {
                EXPECT_NO_THROW(denseUnitatyMatrix_j->
                                unitaryTransformAndReturnMainDiagonal(denseDiagonalizableMatrix_i));
            }
        }
    }

    std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>> denseVectors;
    for (int i = 0; i < denseUnitaryMatrices.size(); ++i) {
        denseVectors.emplace_back(denseDiagonalizableMatrices[i]->diagonalizeValues());
    }

    for (int i = 0; i < denseVectors.size(); ++i) {
        auto& vector_i = denseVectors[i];
        for (int j = 0; j < denseVectors.size(); ++j) {
            auto& vector_j = denseVectors[j];
            if (i != j) {
                EXPECT_ANY_THROW(vector_i->concatenate_with(vector_j));
                EXPECT_ANY_THROW(vector_i->dot(vector_j));
                EXPECT_ANY_THROW(vector_i->element_wise_multiplication(vector_j));
            } else {
                EXPECT_NO_THROW(vector_i->concatenate_with(vector_j));
                EXPECT_NO_THROW(vector_i->dot(vector_j));
                EXPECT_NO_THROW(vector_i->element_wise_multiplication(vector_j));
            }
        }
    }
}

// Check if different linear algebra packages make the same eigen decomposition
TEST(linearAlgebraFactories, eigendecomposition) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    for (size_t size = 2; size < 100; ++size) {
        // construct identical dense diagonalizable matrix:
        auto matrices = generateDenseDiagonalizableMatrices(
            size,
            constructAllDenseTransformAndDiagonalizeFactories(),
            dist,
            rng);
        auto& armaMatrix = matrices[0];
        auto& eigenMatrix = matrices[2];
        // only-values-eigendecomposition:
        {
            // decomposition
            std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>> denseVectors;
            for (const auto& matrix : matrices) {
                denseVectors.emplace_back(matrix->diagonalizeValues());
            }
            // check equality:
            for (size_t i = 0; i < denseVectors.size(); ++i) {
                const auto& denseVector_i = denseVectors[i];
                for (size_t j = 0; j < denseVectors.size(); ++j) {
                    if (i == j) {
                        continue;
                    }
                    const auto& denseVector_j = denseVectors[j];
                    for (size_t k = 0; k < size; ++k) {
                        double epsilon = std::abs(denseVector_i->at(k) * 5e-3);
                        EXPECT_NEAR(denseVector_i->at(k), denseVector_j->at(k), epsilon);
                    }
                }
            }
        }
        // eigendecomposition-with-eigenvectors
        {
            // decomposition
            std::vector<quantum::linear_algebra::EigenCouple> denseEigenCouples;
            for (const auto& matrix : matrices) {
                denseEigenCouples.emplace_back(matrix->diagonalizeValuesVectors());
            }
            // check equality:
            for (size_t i = 0; i < denseEigenCouples.size(); ++i) {
                auto& eigencouple_i = denseEigenCouples[i];
                for (size_t j = 0; j < denseEigenCouples.size(); ++j) {
                    if (i == j) {
                        continue;
                    }
                    auto& eigencouple_j = denseEigenCouples[j];
                    for (size_t k = 0; k < size; ++k) {
                        double epsilon_values = std::abs(eigencouple_i.eigenvalues->at(k) * 5e-3);
                        EXPECT_NEAR(
                            eigencouple_i.eigenvalues->at(k),
                            eigencouple_j.eigenvalues->at(k),
                            epsilon_values);
                        int sign;
                        for (size_t l = 0; l < size; ++l) {
                            if (l == 0) {
                                sign = sign_f(eigencouple_i.eigenvectors->at(k, l),
                                              eigencouple_j.eigenvectors->at(k, l));
                            }
                            // TODO: epsilon
                            EXPECT_NEAR(
                                eigencouple_i.eigenvectors->at(k, l),
                                sign * eigencouple_j.eigenvectors->at(k, l),
                                5e-4);
                        }
                    }
                }
            }
        }
    }
}

// Check if different linear algebra packages make the same unitary transformation
TEST(linearAlgebraFactories, unitary_transformation) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    for (size_t size = 2; size < 100; ++size) {
        auto denseDiagonalizableMatrices = generateDenseDiagonalizableMatrices(
            size,
            constructAllDenseTransformAndDiagonalizeFactories(),
            dist,
            rng);
        auto denseUnitaryMatrices = generateDenseUnitaryMatrices(
            size,
            constructAllDenseTransformAndDiagonalizeFactories(),
            dist,
            rng);

        // unitary transformation:
        std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDiagonalizableMatrix>>
            denseTransformedMatrices;
        for (size_t i = 0; i < denseDiagonalizableMatrices.size(); ++i) {
            const auto& denseDiagonalizableMatrix = denseDiagonalizableMatrices[i];
            const auto& denseUnitaryMatrix = denseUnitaryMatrices[i];

            denseTransformedMatrices.emplace_back(
                denseUnitaryMatrix->unitaryTransform(denseDiagonalizableMatrix));
        }

        // check equality:
        for (size_t i = 0; i < denseTransformedMatrices.size(); ++i) {
            const auto& transformedMatrix_i = denseTransformedMatrices[i];
            for (size_t j = 0; j < denseTransformedMatrices.size(); ++j) {
                if (i == j) {
                    continue;
                }
                const auto& transformedMatrix_j = denseTransformedMatrices[j];
                for (size_t k = 0; k < size; ++k) {
                    for (size_t l = 0; l < size; ++l) {
                        double epsilon = std::max(
                            std::abs(transformedMatrix_i->at(k, l) * 1e-5),
                            3e-2);
                        EXPECT_NEAR(
                            transformedMatrix_i->at(k, l),
                            transformedMatrix_j->at(k, l),
                            epsilon);
                    }
                }
            }
        }
    }
}

// Check if different linear algebra packages make the same unitary transformation
TEST(linearAlgebraFactories, unitary_transformation_and_return_main_diagonal) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(-1000, +1000);

    for (size_t size = 2; size < 100; ++size) {
        auto denseDiagonalizableMatrices = generateDenseDiagonalizableMatrices(
            size,
            constructAllDenseTransformAndDiagonalizeFactories(),
            dist,
            rng);
        auto denseUnitaryMatrices = generateDenseUnitaryMatrices(
            size,
            constructAllDenseTransformAndDiagonalizeFactories(),
            dist,
            rng);

        // unitary transformation and returning main diagonal:
        std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>
            denseTransformedMainDiagonal;
        for (size_t i = 0; i < denseDiagonalizableMatrices.size(); ++i) {
            const auto& denseDiagonalizableMatrix = denseDiagonalizableMatrices[i];
            const auto& denseUnitaryMatrix = denseUnitaryMatrices[i];

            denseTransformedMainDiagonal.emplace_back(
                denseUnitaryMatrix->unitaryTransformAndReturnMainDiagonal(denseDiagonalizableMatrix));
        }

        // check equality:
        for (size_t i = 0; i < denseTransformedMainDiagonal.size(); ++i) {
            const auto& transformedMainDiagonal_i = denseTransformedMainDiagonal[i];
            for (size_t j = 0; j < denseTransformedMainDiagonal.size(); ++j) {
                if (i == j) {
                    continue;
                }
                const auto& transformedMainDiagonal_j = denseTransformedMainDiagonal[j];
                for (size_t k = 0; k < size; ++k) {
                    double epsilon = std::max(
                        std::abs(transformedMainDiagonal_i->at(k) * 1e-4),
                        5e-3);
                    EXPECT_NEAR(
                        transformedMainDiagonal_i->at(k),
                        transformedMainDiagonal_j->at(k),
                        epsilon);
                }
            }
        }
    }
}
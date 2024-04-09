#include "ArmaFactories.h"

#include "ArmaDenseDiagonalizableMatrix.h"
#include "ArmaDenseSemiunitaryMatrix.h"
#include "ArmaDenseVector.h"
#include "ArmaSparseDiagonalizableMatrix.h"
#include "ArmaSparseSemiunitaryMatrix.h"

namespace quantum::linear_algebra {

std::unique_ptr<AbstractDiagonalizableMatrix>
ArmaDenseTransformAndDiagonalizeFactory::createDenseDiagonalizableMatrix(uint32_t size) {
    if (getPrecision() == Precision::SINGLE) {
        auto matrix = std::make_unique<ArmaDenseDiagonalizableMatrix<float>>();
        matrix->resize(size);
        return matrix;
    } else {
        auto matrix = std::make_unique<ArmaDenseDiagonalizableMatrix<double>>();
        matrix->resize(size);
        return matrix;
    }
}

std::unique_ptr<AbstractDiagonalizableMatrix>
ArmaDenseTransformAndDiagonalizeFactory::createSparseDiagonalizableMatrix(uint32_t size) {
    if (getPrecision() == Precision::SINGLE) {
        auto matrix = std::make_unique<ArmaSparseDiagonalizableMatrix<float>>();
        matrix->resize(size);
        return matrix;
    } else {
        auto matrix = std::make_unique<ArmaSparseDiagonalizableMatrix<double>>();
        matrix->resize(size);
        return matrix;
    }
}

std::unique_ptr<AbstractDenseSemiunitaryMatrix>
ArmaDenseTransformAndDiagonalizeFactory::createDenseSemiunitaryMatrix(
    uint32_t cols,
    uint32_t rows) {
    if (getPrecision() == Precision::SINGLE) {
        auto matrix = std::make_unique<ArmaDenseSemiunitaryMatrix<float>>();
        matrix->resize(rows, cols);
        return matrix;
    } else {
        auto matrix = std::make_unique<ArmaDenseSemiunitaryMatrix<double>>();
        matrix->resize(rows, cols);
        return matrix;
    }
}

std::unique_ptr<AbstractDenseVector> ArmaDenseTransformAndDiagonalizeFactory::createVector() {
    if (getPrecision() == Precision::SINGLE) {
        auto vector = std::make_unique<ArmaDenseVector<float>>();
        return vector;
    } else {
        auto vector = std::make_unique<ArmaDenseVector<double>>();
        return vector;
    }
}

std::unique_ptr<AbstractSparseSemiunitaryMatrix>
ArmaSparseTransformFactory::createSparseSemiunitaryMatrix(uint32_t cols, uint32_t rows) {
    auto matrix = std::make_unique<ArmaSparseSemiunitaryMatrix>();
    matrix->resize(cols, rows);
    return matrix;
}

std::unique_ptr<AbstractSymmetricMatrix>
ArmaSparseTransformFactory::createSparseSymmetricMatrix(uint32_t size) {
    auto matrix = std::make_unique<ArmaSparseSymmetricMatrix<double>>();
    matrix->resize(size);
    return matrix;
}
}  // namespace quantum::linear_algebra
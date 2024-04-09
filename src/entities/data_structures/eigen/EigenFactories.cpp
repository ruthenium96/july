#include "EigenFactories.h"

#include "EigenDenseDiagonalizableMatrix.h"
#include "EigenDenseSemiunitaryMatrix.h"
#include "EigenDenseVector.h"
#include "EigenSparseDiagonalizableMatrix.h"

namespace quantum::linear_algebra {

std::unique_ptr<AbstractDiagonalizableMatrix>
EigenDenseTransformAndDiagonalizeFactory::createDenseDiagonalizableMatrix(uint32_t size) {
    if (getPrecision() == Precision::SINGLE) {
        auto matrix = std::make_unique<EigenDenseDiagonalizableMatrix<float>>();
        matrix->resize(size);
        return matrix;
    } else {
        auto matrix = std::make_unique<EigenDenseDiagonalizableMatrix<double>>();
        matrix->resize(size);
        return matrix;
    }
}

std::unique_ptr<AbstractDiagonalizableMatrix>
EigenDenseTransformAndDiagonalizeFactory::createSparseDiagonalizableMatrix(uint32_t size) {
    if (getPrecision() == Precision::SINGLE) {
        auto matrix = std::make_unique<EigenSparseDiagonalizableMatrix<float>>();
        matrix->resize(size);
        return matrix;
    } else {
        auto matrix = std::make_unique<EigenSparseDiagonalizableMatrix<double>>();
        matrix->resize(size);
        return matrix;
    }
}

std::unique_ptr<AbstractDenseSemiunitaryMatrix>
EigenDenseTransformAndDiagonalizeFactory::createDenseSemiunitaryMatrix(
    uint32_t cols,
    uint32_t rows) {
    if (getPrecision() == Precision::SINGLE) {
        auto matrix = std::make_unique<EigenDenseSemiunitaryMatrix<float>>();
        matrix->resize(rows, cols);
        return matrix;
    } else {
        auto matrix = std::make_unique<EigenDenseSemiunitaryMatrix<double>>();
        matrix->resize(rows, cols);
        return matrix;
    }
}

std::unique_ptr<AbstractDenseVector> EigenDenseTransformAndDiagonalizeFactory::createVector() {
    if (getPrecision() == Precision::SINGLE) {
        auto vector = std::make_unique<EigenDenseVector<float>>();
        return vector;
    } else {
        auto vector = std::make_unique<EigenDenseVector<double>>();
        return vector;
    }
}

}  // namespace quantum::linear_algebra
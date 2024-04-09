#include "ArmaLogic.h"

#include "ArmaDenseDiagonalizableMatrix.h"
#include "ArmaDenseSemiunitaryMatrix.h"
#include "ArmaDenseVector.h"
#include "ArmaSparseDiagonalizableMatrix.h"

namespace quantum::linear_algebra {
template <typename T>
std::unique_ptr<AbstractDenseVector> ArmaLogic<T>::diagonalizeValues(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& diagonalizableMatrix) const {
    return diagonalizeValues(*diagonalizableMatrix);
}

template <typename T>
std::unique_ptr<AbstractDenseVector>
ArmaLogic<T>::diagonalizeValues(const AbstractDiagonalizableMatrix& diagonalizableMatrix) const {
    auto eigenvalues_ = std::make_unique<ArmaDenseVector<T>>();

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseDiagonalizableMatrix<T>*>(&diagonalizableMatrix)) {
        arma::eig_sym(
            eigenvalues_->modifyDenseVector(),
            maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix());
    } else if (
        auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseDiagonalizableMatrix<T>*>(&diagonalizableMatrix)) {
        auto sparseToDenseCopy = arma::Mat<T>(maybeSparseSymmetricMatrix->getSparseSymmetricMatrix());
        arma::eig_sym(eigenvalues_->modifyDenseVector(), sparseToDenseCopy);
    } else {
        throw std::bad_cast();
    }

    return eigenvalues_;
}

template <typename T>
EigenCouple ArmaLogic<T>::diagonalizeValuesVectors(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix) const {
    return diagonalizeValuesVectors(*symmetricMatrix);
}

template <typename T>
EigenCouple ArmaLogic<T>::diagonalizeValuesVectors(
    const AbstractDiagonalizableMatrix& diagonalizableMatrix) const {
    auto eigenvalues_ = std::make_unique<ArmaDenseVector<T>>();
    auto eigenvectors_ = std::make_unique<ArmaDenseSemiunitaryMatrix<T>>();

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseDiagonalizableMatrix<T>*>(&diagonalizableMatrix)) {
        arma::eig_sym(
            eigenvalues_->modifyDenseVector(),
            eigenvectors_->modifyDenseSemiunitaryMatrix(),
            maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix());
    } else if (
        auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseDiagonalizableMatrix<T>*>(&diagonalizableMatrix)) {
        auto sparseToDenseCopy = arma::Mat<T>(maybeSparseSymmetricMatrix->getSparseSymmetricMatrix());
        arma::eig_sym(
            eigenvalues_->modifyDenseVector(),
            eigenvectors_->modifyDenseSemiunitaryMatrix(),
            sparseToDenseCopy);
    } else {
        throw std::bad_cast();
    }

    EigenCouple answer;
    answer.eigenvalues = std::move(eigenvalues_);
    answer.eigenvectors = std::move(eigenvectors_);
    return answer;
}

template <typename T>
std::unique_ptr<AbstractDenseVector> ArmaLogic<T>::unitaryTransformAndReturnMainDiagonal(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
    const std::unique_ptr<AbstractDenseSemiunitaryMatrix>& denseSemiunitaryMatrix) const {
    return unitaryTransformAndReturnMainDiagonal(symmetricMatrix, *denseSemiunitaryMatrix);
}

template <typename T>
std::unique_ptr<AbstractDenseVector> ArmaLogic<T>::unitaryTransformAndReturnMainDiagonal(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
    const AbstractDenseSemiunitaryMatrix& denseSemiunitaryMatrix) const {
    auto main_diagonal = std::make_unique<ArmaDenseVector<T>>();

    auto maybeDenseSemiunitaryMatrix =
        dynamic_cast<const ArmaDenseSemiunitaryMatrix<T>*>(&denseSemiunitaryMatrix);
    if (maybeDenseSemiunitaryMatrix == nullptr) {
        throw std::bad_cast();
    }

    if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        arma::Mat<T> firstMultiplicationResult =
            maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().t()
            * maybeSparseSymmetricMatrix->getSparseSymmetricMatrix();
        main_diagonal->modifyDenseVector() = arma::diagvec(
            firstMultiplicationResult * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix());
        return std::move(main_diagonal);
    }

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        arma::Mat<T> firstMultiplicationResult =
            maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().t()
            * maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix();
        main_diagonal->modifyDenseVector() = arma::diagvec(
            firstMultiplicationResult * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix());
        return std::move(main_diagonal);
    }

    throw std::bad_cast();
}

template <typename T>
std::unique_ptr<AbstractDiagonalizableMatrix> ArmaLogic<T>::unitaryTransform(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
    const AbstractDenseSemiunitaryMatrix& denseSemiunitaryMatrix) const {
    auto result = std::make_unique<ArmaDenseDiagonalizableMatrix<T>>();

    auto maybeDenseSemiunitaryMatrix =
        dynamic_cast<const ArmaDenseSemiunitaryMatrix<T>*>(&denseSemiunitaryMatrix);
    if (maybeDenseSemiunitaryMatrix == nullptr) {
        throw std::bad_cast();
    }

    if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        result->modifyDenseDiagonalizableMatrix() =
            maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().t()
            * maybeSparseSymmetricMatrix->getSparseSymmetricMatrix()
            * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix();
        return std::move(result);
    }

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        result->modifyDenseDiagonalizableMatrix() =
            maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().t()
            * maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix()
            * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix();
        return std::move(result);
    }

    throw std::bad_cast();
}

template class ArmaLogic<double>;
template class ArmaLogic<float>;
}  // namespace quantum::linear_algebra
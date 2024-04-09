#include "EigenLogic.h"

#include "EigenDenseDiagonalizableMatrix.h"
#include "EigenDenseSemiunitaryMatrix.h"
#include "EigenDenseVector.h"
#include "EigenSparseDiagonalizableMatrix.h"

namespace quantum::linear_algebra {

template <typename T>
std::unique_ptr<AbstractDenseVector> EigenLogic<T>::unitaryTransformAndReturnMainDiagonal(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
    const AbstractDenseSemiunitaryMatrix& denseSemiunitaryMatrix) const {
    auto main_diagonal = std::make_unique<EigenDenseVector<T>>();

    auto maybeDenseSemiunitaryMatrix =
        dynamic_cast<const EigenDenseSemiunitaryMatrix<T>*>(&denseSemiunitaryMatrix);
    if (maybeDenseSemiunitaryMatrix == nullptr) {
        throw std::bad_cast();
    }

    if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const EigenSparseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        Eigen::Matrix<T, -1, -1> firstMultiplicationResult =
            maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().transpose()
            * maybeSparseSymmetricMatrix->getSparseDiagonalizableMatrix().transpose();
        main_diagonal->resize(maybeSparseSymmetricMatrix->size());
#pragma omp parallel for shared( \
        main_diagonal, \
            firstMultiplicationResult, \
            maybeDenseSemiunitaryMatrix) default(shared)
        for (size_t i = 0; i < main_diagonal->size(); ++i) {
            auto value = firstMultiplicationResult.row(i)
                * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().col(i);
            main_diagonal->modifyDenseVector().coeffRef(i) = value;
        }
        return std::move(main_diagonal);
    }

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const EigenDenseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        Eigen::Matrix<T, -1, -1> firstMultiplicationResult =
            maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().transpose()
            * maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix();
        main_diagonal->resize(maybeDenseSymmetricMatrix->size());
#pragma omp parallel for shared( \
        main_diagonal, \
            firstMultiplicationResult, \
            maybeDenseSemiunitaryMatrix) default(none)
        for (size_t i = 0; i < main_diagonal->size(); ++i) {
            auto value = firstMultiplicationResult.row(i)
                * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().col(i);
            main_diagonal->modifyDenseVector().coeffRef(i) = value;
        }
        return std::move(main_diagonal);
    }

    throw std::bad_cast();
}

template <typename T>
std::unique_ptr<AbstractDiagonalizableMatrix> EigenLogic<T>::unitaryTransform(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
    const AbstractDenseSemiunitaryMatrix& denseSemiunitaryMatrix) const {
    auto result = std::make_unique<EigenDenseDiagonalizableMatrix<T>>();

    auto maybeDenseSemiunitaryMatrix =
        dynamic_cast<const EigenDenseSemiunitaryMatrix<T>*>(&denseSemiunitaryMatrix);
    if (maybeDenseSemiunitaryMatrix == nullptr) {
        throw std::bad_cast();
    }

    if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const EigenSparseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        // TODO: check if all transpositions are efficient
        result->modifyDenseDiagonalizableMatrix() =
            maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().transpose()
            * maybeSparseSymmetricMatrix->getSparseDiagonalizableMatrix().transpose()
            * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix();
        return std::move(result);
    }

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const EigenDenseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        result->modifyDenseDiagonalizableMatrix() =
            maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().transpose()
            * maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix()
            * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix();
        return std::move(result);
    }

    throw std::bad_cast();
}

template <typename T>
std::unique_ptr<AbstractDenseVector>
EigenLogic<T>::diagonalizeValues(const AbstractDiagonalizableMatrix& symmetricMatrix) const {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, -1, -1>> es;

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const EigenDenseDiagonalizableMatrix<T>*>(&symmetricMatrix)) {
        es.compute(
            maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix(),
            Eigen::EigenvaluesOnly);
    } else if (
        auto maybeSparseSymmetricMatrix =
            dynamic_cast<const EigenSparseDiagonalizableMatrix<T>*>(&symmetricMatrix)) {
        Eigen::Matrix<T, -1, -1> sparseToDenseCopy =
            maybeSparseSymmetricMatrix->getSparseDiagonalizableMatrix();
        es.compute(sparseToDenseCopy, Eigen::EigenvaluesOnly);
    } else {
        throw std::bad_cast();
    }

    auto eigenvalues_ = std::make_unique<EigenDenseVector<T>>();
    eigenvalues_->modifyDenseVector() = es.eigenvalues();

    return eigenvalues_;
}

template <typename T>
EigenCouple
EigenLogic<T>::diagonalizeValuesVectors(const AbstractDiagonalizableMatrix& symmetricMatrix) const {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, -1, -1>> es;

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const EigenDenseDiagonalizableMatrix<T>*>(&symmetricMatrix)) {
        es.compute(
            maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix(),
            Eigen::ComputeEigenvectors);
    } else if (
        auto maybeSparseSymmetricMatrix =
            dynamic_cast<const EigenSparseDiagonalizableMatrix<T>*>(&symmetricMatrix)) {
        Eigen::Matrix<T, -1, -1> sparseToDenseCopy =
            maybeSparseSymmetricMatrix->getSparseDiagonalizableMatrix();
        es.compute(sparseToDenseCopy, Eigen::ComputeEigenvectors);
    } else {
        throw std::bad_cast();
    }

    auto eigenvalues_ = std::make_unique<EigenDenseVector<T>>();
    auto eigenvectors_ = std::make_unique<EigenDenseSemiunitaryMatrix<T>>();
    eigenvalues_->modifyDenseVector() = es.eigenvalues();
    eigenvectors_->modifyDenseSemiunitaryMatrix() = es.eigenvectors();

    EigenCouple answer;
    answer.eigenvalues = std::move(eigenvalues_);
    answer.eigenvectors = std::move(eigenvectors_);
    return answer;
}

template class EigenLogic<double>;
template class EigenLogic<float>;
}  // namespace quantum::linear_algebra
#include "ArmaSparseDiagonalizableMatrix.h"

#include "ArmaLogic.h"

namespace quantum::linear_algebra {

template <typename T>
void ArmaSparseDiagonalizableMatrix<T>::add_to_position(double value, uint32_t i, uint32_t j) {
    sparseDiagonalizableMatrix_.at(i, j) += value;
    if (i != j) {
        sparseDiagonalizableMatrix_(j, i) += value;
    }
}

template <typename T>
void ArmaSparseDiagonalizableMatrix<T>::resize(uint32_t matrix_in_space_basis_size_i) {
    sparseDiagonalizableMatrix_.resize(matrix_in_space_basis_size_i, matrix_in_space_basis_size_i);
}

template <typename T>
EigenCouple ArmaSparseDiagonalizableMatrix<T>::diagonalizeValuesVectors() const {
    ArmaLogic<T> logic;

    return logic.diagonalizeValuesVectors(*this);
}

template <typename T>
std::unique_ptr<AbstractDenseVector> ArmaSparseDiagonalizableMatrix<T>::diagonalizeValues() const {
    ArmaLogic<T> logic;

    return logic.diagonalizeValues(*this);
}

template <typename T>
std::unique_ptr<AbstractDiagonalizableMatrix>
ArmaSparseDiagonalizableMatrix<T>::multiply_by(double multiplier) const {
    auto answer = std::make_unique<ArmaSparseDiagonalizableMatrix>();
    answer->resize(size());
    answer->sparseDiagonalizableMatrix_ = sparseDiagonalizableMatrix_ * multiplier;
    return answer;
}

template <typename T>
uint32_t ArmaSparseDiagonalizableMatrix<T>::size() const {
    return sparseDiagonalizableMatrix_.n_rows;
}

template <typename T>
double ArmaSparseDiagonalizableMatrix<T>::at(uint32_t i, uint32_t j) const {
    return sparseDiagonalizableMatrix_.at(i, j);
}

template <typename T>
void ArmaSparseDiagonalizableMatrix<T>::print(std::ostream& os) const {
    os << sparseDiagonalizableMatrix_ << std::endl;
}

template <typename T>
const arma::SpMat<T>& ArmaSparseDiagonalizableMatrix<T>::getSparseSymmetricMatrix() const {
    return sparseDiagonalizableMatrix_;
}

template <typename T>
arma::SpMat<T>& ArmaSparseDiagonalizableMatrix<T>::modifySparseSymmetricMatrix() {
    return sparseDiagonalizableMatrix_;
}

template class ArmaSparseDiagonalizableMatrix<double>;
template class ArmaSparseDiagonalizableMatrix<float>;
}  // namespace quantum::linear_algebra
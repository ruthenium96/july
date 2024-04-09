#include "EigenSparseDiagonalizableMatrix.h"

#include "EigenLogic.h"

namespace quantum::linear_algebra {

template <typename T>
void EigenSparseDiagonalizableMatrix<T>::add_to_position(double value, uint32_t i, uint32_t j) {
    sparseDiagonalizableMatrix_.coeffRef(i, j) += value;
    if (i != j) {
        sparseDiagonalizableMatrix_.coeffRef(j, i) += value;
    }
}

template <typename T>
EigenCouple EigenSparseDiagonalizableMatrix<T>::diagonalizeValuesVectors() const {
    EigenLogic<T> eigenLogic;

    return eigenLogic.diagonalizeValuesVectors(*this);
}

template <typename T>
std::unique_ptr<AbstractDenseVector> EigenSparseDiagonalizableMatrix<T>::diagonalizeValues() const {
    EigenLogic<T> eigenLogic;

    return eigenLogic.diagonalizeValues(*this);
}

template <typename T>
std::unique_ptr<AbstractDiagonalizableMatrix>
EigenSparseDiagonalizableMatrix<T>::multiply_by(double multiplier) const {
    auto answer = std::make_unique<EigenSparseDiagonalizableMatrix>();
    answer->resize(size());

    answer->sparseDiagonalizableMatrix_ = sparseDiagonalizableMatrix_ * multiplier;

    return answer;
}

template <typename T>
uint32_t EigenSparseDiagonalizableMatrix<T>::size() const {
    return sparseDiagonalizableMatrix_.innerSize();
}

template <typename T>
double EigenSparseDiagonalizableMatrix<T>::at(uint32_t i, uint32_t j) const {
    return sparseDiagonalizableMatrix_.coeff(i, j);
}

template <typename T>
void EigenSparseDiagonalizableMatrix<T>::resize(size_t size) {
    sparseDiagonalizableMatrix_.resize(size, size);
}

template <typename T>
const Eigen::SparseMatrix<T>&
EigenSparseDiagonalizableMatrix<T>::getSparseDiagonalizableMatrix() const {
    return sparseDiagonalizableMatrix_;
}

template <typename T>
Eigen::SparseMatrix<T>& EigenSparseDiagonalizableMatrix<T>::modifySparseDiagonalizableMatrix() {
    return sparseDiagonalizableMatrix_;
}

template <typename T>
void EigenSparseDiagonalizableMatrix<T>::print(std::ostream& os) const {
    os << sparseDiagonalizableMatrix_ << std::endl;
}

template class EigenSparseDiagonalizableMatrix<double>;
template class EigenSparseDiagonalizableMatrix<float>;
}  // namespace quantum::linear_algebra
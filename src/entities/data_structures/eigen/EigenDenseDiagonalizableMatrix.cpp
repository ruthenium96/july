#include "EigenDenseDiagonalizableMatrix.h"

#include "EigenLogic.h"

namespace quantum::linear_algebra {

template <typename T>
void EigenDenseDiagonalizableMatrix<T>::add_to_position(double value, uint32_t i, uint32_t j) {
    denseDiagonalizableMatrix_(i, j) += value;
    if (i != j) {
        denseDiagonalizableMatrix_(j, i) += value;
    }
}

template <typename T>
EigenCouple EigenDenseDiagonalizableMatrix<T>::diagonalizeValuesVectors() const {
    EigenLogic<T> eigenLogic;

    return eigenLogic.diagonalizeValuesVectors(*this);
}

template <typename T>
std::unique_ptr<AbstractDenseVector> EigenDenseDiagonalizableMatrix<T>::diagonalizeValues() const {
    EigenLogic<T> eigenLogic;

    return eigenLogic.diagonalizeValues(*this);
}

template <typename T>
std::unique_ptr<AbstractDiagonalizableMatrix>
EigenDenseDiagonalizableMatrix<T>::multiply_by(double multiplier) const {
    auto answer = std::make_unique<EigenDenseDiagonalizableMatrix>();
    answer->resize(size());

    answer->denseDiagonalizableMatrix_ = denseDiagonalizableMatrix_ * multiplier;

    return answer;
}

template <typename T>
uint32_t EigenDenseDiagonalizableMatrix<T>::size() const {
    return denseDiagonalizableMatrix_.cols();
}

template <typename T>
double EigenDenseDiagonalizableMatrix<T>::at(uint32_t i, uint32_t j) const {
    return denseDiagonalizableMatrix_(i, j);
}

template <typename T>
void EigenDenseDiagonalizableMatrix<T>::print(std::ostream& os) const {
    os << denseDiagonalizableMatrix_ << std::endl;
}

template <typename T>
void EigenDenseDiagonalizableMatrix<T>::resize(uint32_t size) {
    denseDiagonalizableMatrix_.resize(size, size);
    denseDiagonalizableMatrix_.fill(0);
}

template <typename T>
const Eigen::Matrix<T, -1, -1>& EigenDenseDiagonalizableMatrix<T>::getDenseDiagonalizableMatrix() const {
    return denseDiagonalizableMatrix_;
}

template <typename T>
Eigen::Matrix<T, -1, -1>& EigenDenseDiagonalizableMatrix<T>::modifyDenseDiagonalizableMatrix() {
    return denseDiagonalizableMatrix_;
}

template class EigenDenseDiagonalizableMatrix<double>;
template class EigenDenseDiagonalizableMatrix<float>;
}  // namespace quantum::linear_algebra
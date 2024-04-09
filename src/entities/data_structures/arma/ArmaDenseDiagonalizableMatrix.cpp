#include "ArmaDenseDiagonalizableMatrix.h"

#include "ArmaLogic.h"

namespace quantum::linear_algebra {

template <typename T>
void ArmaDenseDiagonalizableMatrix<T>::add_to_position(double value, uint32_t i, uint32_t j) {
    denseDiagonalizableMatrix_(i, j) += value;
    if (i != j) {
        denseDiagonalizableMatrix_(j, i) += value;
    }
}

template <typename T>
EigenCouple ArmaDenseDiagonalizableMatrix<T>::diagonalizeValuesVectors() const {
    ArmaLogic<T> logic;

    return logic.diagonalizeValuesVectors(*this);
}

template <typename T>
std::unique_ptr<AbstractDenseVector> ArmaDenseDiagonalizableMatrix<T>::diagonalizeValues() const {
    ArmaLogic<T> logic;

    return logic.diagonalizeValues(*this);
}

template <typename T>
std::unique_ptr<AbstractDiagonalizableMatrix>
ArmaDenseDiagonalizableMatrix<T>::multiply_by(double multiplier) const {
    auto answer = std::make_unique<ArmaDenseDiagonalizableMatrix>();
    answer->resize(size());
    answer->denseDiagonalizableMatrix_ = denseDiagonalizableMatrix_ * multiplier;
    return answer;
}

template <typename T>
uint32_t ArmaDenseDiagonalizableMatrix<T>::size() const {
    return denseDiagonalizableMatrix_.n_rows;
}

template <typename T>
double ArmaDenseDiagonalizableMatrix<T>::at(uint32_t i, uint32_t j) const {
    return denseDiagonalizableMatrix_.at(i, j);
}

template <typename T>
void ArmaDenseDiagonalizableMatrix<T>::print(std::ostream& os) const {
    os << denseDiagonalizableMatrix_ << std::endl;
}

template <typename T>
void ArmaDenseDiagonalizableMatrix<T>::resize(uint32_t size) {
    denseDiagonalizableMatrix_.resize(size, size);
}

template <typename T>
const arma::Mat<T>& ArmaDenseDiagonalizableMatrix<T>::getDenseDiagonalizableMatrix() const {
    return denseDiagonalizableMatrix_;
}

template <typename T>
arma::Mat<T>& ArmaDenseDiagonalizableMatrix<T>::modifyDenseDiagonalizableMatrix() {
    return denseDiagonalizableMatrix_;
}

template class ArmaDenseDiagonalizableMatrix<double>;
template class ArmaDenseDiagonalizableMatrix<float>;
}  // namespace quantum::linear_algebra
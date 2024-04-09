#include "ArmaDenseSemiunitaryMatrix.h"

#include "ArmaLogic.h"

namespace quantum::linear_algebra {

template <typename T>
uint32_t ArmaDenseSemiunitaryMatrix<T>::size_rows() const {
    return denseSemiunitaryMatrix_.n_rows;
}

template <typename T>
uint32_t ArmaDenseSemiunitaryMatrix<T>::size_cols() const {
    return denseSemiunitaryMatrix_.n_cols;
}

template <typename T>
double ArmaDenseSemiunitaryMatrix<T>::at(uint32_t i, uint32_t j) const {
    return denseSemiunitaryMatrix_.at(j, i);
}

template <typename T>
void ArmaDenseSemiunitaryMatrix<T>::print(std::ostream& os) const {
    os << denseSemiunitaryMatrix_ << std::endl;
}

template <typename T>
void ArmaDenseSemiunitaryMatrix<T>::resize(size_t size_rows, size_t size_cols) {
    denseSemiunitaryMatrix_.resize(size_rows, size_cols);
    denseSemiunitaryMatrix_.fill(arma::fill::zeros);
}

template <typename T>
std::unique_ptr<AbstractDenseVector>
ArmaDenseSemiunitaryMatrix<T>::unitaryTransformAndReturnMainDiagonal(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& matrix_to_transform) const {
    ArmaLogic<T> logic;

    return logic.unitaryTransformAndReturnMainDiagonal(matrix_to_transform, *this);
}

template <typename T>
std::unique_ptr<AbstractDiagonalizableMatrix> ArmaDenseSemiunitaryMatrix<T>::unitaryTransform(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& matrix_to_transform) const {
    ArmaLogic<T> logic;

    return logic.unitaryTransform(matrix_to_transform, *this);
}

template <typename T>
const arma::Mat<T>& ArmaDenseSemiunitaryMatrix<T>::getDenseSemiunitaryMatrix() const {
    return denseSemiunitaryMatrix_;
}

template <typename T>
arma::Mat<T>& ArmaDenseSemiunitaryMatrix<T>::modifyDenseSemiunitaryMatrix() {
    return denseSemiunitaryMatrix_;
}

template <typename T>
void ArmaDenseSemiunitaryMatrix<T>::add_to_position(double value, uint32_t i, uint32_t j) {
    denseSemiunitaryMatrix_.at(j, i) += value;
}

template <typename T>
void ArmaDenseSemiunitaryMatrix<T>::normalize() {
    denseSemiunitaryMatrix_ = arma::normalise(denseSemiunitaryMatrix_, 2, 0);
}

template class ArmaDenseSemiunitaryMatrix<double>;
template class ArmaDenseSemiunitaryMatrix<float>;
}  // namespace quantum::linear_algebra
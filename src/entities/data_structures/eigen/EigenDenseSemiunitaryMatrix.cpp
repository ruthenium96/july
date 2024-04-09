#include "EigenDenseSemiunitaryMatrix.h"

#include "EigenDenseVector.h"
#include "EigenLogic.h"

namespace quantum::linear_algebra {

template <typename T>
uint32_t EigenDenseSemiunitaryMatrix<T>::size_rows() const {
    return denseSemiunitaryMatrix_.rows();
}

template <typename T>
uint32_t EigenDenseSemiunitaryMatrix<T>::size_cols() const {
    return denseSemiunitaryMatrix_.cols();
}

template <typename T>
void EigenDenseSemiunitaryMatrix<T>::resize(size_t size_rows, size_t size_cols) {
    denseSemiunitaryMatrix_.resize(size_rows, size_cols);
    denseSemiunitaryMatrix_.fill(0);
}

template <typename T>
double EigenDenseSemiunitaryMatrix<T>::at(uint32_t i, uint32_t j) const {
    return denseSemiunitaryMatrix_(j, i);
}

template <typename T>
void EigenDenseSemiunitaryMatrix<T>::print(std::ostream& os) const {
    os << denseSemiunitaryMatrix_ << std::endl;
}

template <typename T>
std::unique_ptr<AbstractDenseVector>
EigenDenseSemiunitaryMatrix<T>::unitaryTransformAndReturnMainDiagonal(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& matrix_to_transform) const {
    EigenLogic<T> eigenLogic;

    return eigenLogic.unitaryTransformAndReturnMainDiagonal(matrix_to_transform, *this);
}

template <typename T>
std::unique_ptr<AbstractDiagonalizableMatrix> EigenDenseSemiunitaryMatrix<T>::unitaryTransform(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& matrix_to_transform) const {
    EigenLogic<T> eigenLogic;

    return eigenLogic.unitaryTransform(matrix_to_transform, *this);
}

template <typename T>
const Eigen::Matrix<T, -1, -1>& EigenDenseSemiunitaryMatrix<T>::getDenseSemiunitaryMatrix() const {
    return denseSemiunitaryMatrix_;
}

template <typename T>
Eigen::Matrix<T, -1, -1>& EigenDenseSemiunitaryMatrix<T>::modifyDenseSemiunitaryMatrix() {
    return denseSemiunitaryMatrix_;
}

template <typename T>
void EigenDenseSemiunitaryMatrix<T>::add_to_position(double value, uint32_t i, uint32_t j) {
    denseSemiunitaryMatrix_(j, i) += value;
}

template <typename T>
void EigenDenseSemiunitaryMatrix<T>::normalize() {
    for (int i = 0; i < denseSemiunitaryMatrix_.cols(); i++) {
        denseSemiunitaryMatrix_.col(i).normalize();
    }
}

template class EigenDenseSemiunitaryMatrix<double>;
template class EigenDenseSemiunitaryMatrix<float>;
}  // namespace quantum::linear_algebra
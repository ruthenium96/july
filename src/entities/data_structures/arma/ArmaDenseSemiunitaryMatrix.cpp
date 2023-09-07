#include "ArmaDenseSemiunitaryMatrix.h"

#include "ArmaLogic.h"

namespace quantum::linear_algebra {
uint32_t ArmaDenseSemiunitaryMatrix::size_rows() const {
    return denseSemiunitaryMatrix_.n_rows;
}

uint32_t ArmaDenseSemiunitaryMatrix::size_cols() const {
    return denseSemiunitaryMatrix_.n_cols;
}

double ArmaDenseSemiunitaryMatrix::at(uint32_t i, uint32_t j) const {
    return denseSemiunitaryMatrix_.at(j, i);
}

void ArmaDenseSemiunitaryMatrix::print(std::ostream& os) const {
    os << denseSemiunitaryMatrix_ << std::endl;
}

void ArmaDenseSemiunitaryMatrix::resize(size_t size_rows, size_t size_cols) {
    denseSemiunitaryMatrix_.resize(size_rows, size_cols);
    denseSemiunitaryMatrix_.fill(arma::fill::zeros);
}

std::unique_ptr<AbstractDenseVector>
ArmaDenseSemiunitaryMatrix::unitaryTransformAndReturnMainDiagonal(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& matrix_to_transform) const {
    ArmaLogic logic;

    return logic.unitaryTransformAndReturnMainDiagonal(matrix_to_transform, *this);
}

std::unique_ptr<AbstractDiagonalizableMatrix> ArmaDenseSemiunitaryMatrix::unitaryTransform(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& matrix_to_transform) const {
    ArmaLogic logic;

    return logic.unitaryTransform(matrix_to_transform, *this);
}

const arma::dmat& ArmaDenseSemiunitaryMatrix::getDenseSemiunitaryMatrix() const {
    return denseSemiunitaryMatrix_;
}

arma::dmat& ArmaDenseSemiunitaryMatrix::modifyDenseSemiunitaryMatrix() {
    return denseSemiunitaryMatrix_;
}

void ArmaDenseSemiunitaryMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    denseSemiunitaryMatrix_.at(j, i) += value;
}

void ArmaDenseSemiunitaryMatrix::normalize() {
    denseSemiunitaryMatrix_ = arma::normalise(denseSemiunitaryMatrix_, 2, 0);
}
}  // namespace quantum::linear_algebra
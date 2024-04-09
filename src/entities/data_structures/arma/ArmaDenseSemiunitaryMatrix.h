#ifndef SPINNER_ARMADENSESEMIUNITARYMATRIX_H
#define SPINNER_ARMADENSESEMIUNITARYMATRIX_H

#include <armadillo>

#include "src/entities/data_structures/AbstractDenseSemiunitaryMatrix.h"

namespace quantum::linear_algebra {
template <typename T>
class ArmaDenseSemiunitaryMatrix: public AbstractDenseSemiunitaryMatrix {
  public:
    uint32_t size_rows() const override;
    uint32_t size_cols() const override;
    double at(uint32_t i, uint32_t j) const override;
    void add_to_position(double value, uint32_t i, uint32_t j) override;

    void print(std::ostream& os) const override;
    std::unique_ptr<AbstractDenseVector> unitaryTransformAndReturnMainDiagonal(
        const std::unique_ptr<AbstractDiagonalizableMatrix>& matrix_to_transform) const override;
    std::unique_ptr<AbstractDiagonalizableMatrix> unitaryTransform(
        const std::unique_ptr<AbstractDiagonalizableMatrix>& matrix_to_transform) const override;
    void resize(size_t size_rows, size_t size_cols);
    const arma::Mat<T>& getDenseSemiunitaryMatrix() const;
    arma::Mat<T>& modifyDenseSemiunitaryMatrix();
    void normalize() override;

  private:
    arma::Mat<T> denseSemiunitaryMatrix_;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_ARMADENSESEMIUNITARYMATRIX_H

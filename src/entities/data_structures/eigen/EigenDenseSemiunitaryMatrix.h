#ifndef SPINNER_EIGENDENSESEMIUNITARYMATRIX_H
#define SPINNER_EIGENDENSESEMIUNITARYMATRIX_H

#include <Eigen/Core>

#include "src/entities/data_structures/AbstractDenseSemiunitaryMatrix.h"

namespace quantum::linear_algebra {
template <typename T>
class EigenDenseSemiunitaryMatrix: public AbstractDenseSemiunitaryMatrix {
  public:
    uint32_t size_rows() const override;
    uint32_t size_cols() const override;
    void resize(size_t size_rows, size_t size_cols);
    double at(uint32_t i, uint32_t j) const override;
    void add_to_position(double value, uint32_t i, uint32_t j) override;

    void print(std::ostream& os) const override;
    std::unique_ptr<AbstractDenseVector> unitaryTransformAndReturnMainDiagonal(
        const std::unique_ptr<AbstractDiagonalizableMatrix>& matrix_to_transform) const override;
    std::unique_ptr<AbstractDiagonalizableMatrix> unitaryTransform(
        const std::unique_ptr<AbstractDiagonalizableMatrix>& matrix_to_transform) const override;

    const Eigen::Matrix<T, -1, -1>& getDenseSemiunitaryMatrix() const;
    Eigen::Matrix<T, -1, -1>& modifyDenseSemiunitaryMatrix();
    void normalize() override;

  private:
    Eigen::Matrix<T, -1, -1> denseSemiunitaryMatrix_;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_EIGENDENSESEMIUNITARYMATRIX_H

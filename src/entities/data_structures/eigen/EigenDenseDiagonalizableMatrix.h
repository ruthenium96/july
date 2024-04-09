#ifndef SPINNER_EIGENDENSEDIAGONALIZABLEMATRIX_H
#define SPINNER_EIGENDENSEDIAGONALIZABLEMATRIX_H

#include <Eigen/Core>

#include "src/entities/data_structures/AbstractDiagonalizableMatrix.h"

namespace quantum::linear_algebra {
template <typename T>
class EigenDenseDiagonalizableMatrix: public AbstractDiagonalizableMatrix {
  public:
    void add_to_position(double value, uint32_t i, uint32_t j) override;
    EigenCouple diagonalizeValuesVectors() const override;
    std::unique_ptr<AbstractDenseVector> diagonalizeValues() const override;
    std::unique_ptr<AbstractDiagonalizableMatrix> multiply_by(double multiplier) const override;
    uint32_t size() const override;
    double at(uint32_t i, uint32_t j) const override;
    void print(std::ostream& os) const override;
    void resize(uint32_t size);

    const Eigen::Matrix<T, -1, -1>& getDenseDiagonalizableMatrix() const;
    Eigen::Matrix<T, -1, -1>& modifyDenseDiagonalizableMatrix();

  private:
    Eigen::Matrix<T, -1, -1> denseDiagonalizableMatrix_;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_EIGENDENSEDIAGONALIZABLEMATRIX_H

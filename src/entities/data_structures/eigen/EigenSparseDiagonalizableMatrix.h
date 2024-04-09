#ifndef SPINNER_EIGENSPARSEDIAGONALIZABLEMATRIX_H
#define SPINNER_EIGENSPARSEDIAGONALIZABLEMATRIX_H

#include "Eigen/Sparse"
#include "src/entities/data_structures/AbstractDiagonalizableMatrix.h"

namespace quantum::linear_algebra {
template <typename T>
class EigenSparseDiagonalizableMatrix: public AbstractDiagonalizableMatrix {
  public:
    void add_to_position(double value, uint32_t i, uint32_t j) override;
    EigenCouple diagonalizeValuesVectors() const override;
    std::unique_ptr<AbstractDenseVector> diagonalizeValues() const override;
    std::unique_ptr<AbstractDiagonalizableMatrix> multiply_by(double multiplier) const override;
    uint32_t size() const override;
    double at(uint32_t i, uint32_t j) const override;
    void print(std::ostream& os) const override;
    ~EigenSparseDiagonalizableMatrix() override = default;
    void resize(size_t size);

    const Eigen::SparseMatrix<T>& getSparseDiagonalizableMatrix() const;
    Eigen::SparseMatrix<T>& modifySparseDiagonalizableMatrix();

  private:
    Eigen::SparseMatrix<T> sparseDiagonalizableMatrix_;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_EIGENSPARSEDIAGONALIZABLEMATRIX_H

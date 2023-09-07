#ifndef SPINNER_ARMADENSEDIAGONALIZABLEMATRIX_H
#define SPINNER_ARMADENSEDIAGONALIZABLEMATRIX_H

#include <armadillo>

#include "src/entities/data_structures/AbstractDiagonalizableMatrix.h"

namespace quantum::linear_algebra {
class ArmaDenseDiagonalizableMatrix: public AbstractDiagonalizableMatrix {
  public:
    void add_to_position(double value, uint32_t i, uint32_t j) override;
    EigenCouple diagonalizeValuesVectors() const override;
    std::unique_ptr<AbstractDenseVector> diagonalizeValues() const override;
    uint32_t size() const override;
    double at(uint32_t i, uint32_t j) const override;
    void print(std::ostream& os) const override;
    void resize(uint32_t size);
    const arma::dmat& getDenseDiagonalizableMatrix() const;
    arma::dmat& modifyDenseDiagonalizableMatrix();

  private:
    arma::dmat denseDiagonalizableMatrix_;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_ARMADENSEDIAGONALIZABLEMATRIX_H

#ifndef SPINNER_ARMALOGIC_H
#define SPINNER_ARMALOGIC_H

#include <memory>

#include "src/entities/data_structures/AbstractDenseSemiunitaryMatrix.h"
#include "src/entities/data_structures/AbstractDenseVector.h"
#include "src/entities/data_structures/AbstractSparseSemiunitaryMatrix.h"
#include "src/entities/data_structures/AbstractSymmetricMatrix.h"

namespace quantum::linear_algebra {
template <typename T>
class ArmaLogic {
  public:
    std::unique_ptr<AbstractDenseVector>
    diagonalizeValues(const AbstractDiagonalizableMatrix& diagonalizableMatrix) const;
    std::unique_ptr<AbstractDenseVector> diagonalizeValues(
        const std::unique_ptr<AbstractDiagonalizableMatrix>& diagonalizableMatrix) const;

    EigenCouple
    diagonalizeValuesVectors(const AbstractDiagonalizableMatrix& diagonalizableMatrix) const;
    EigenCouple diagonalizeValuesVectors(
        const std::unique_ptr<AbstractDiagonalizableMatrix>& diagonalizableMatrix) const;

    std::unique_ptr<AbstractDenseVector> unitaryTransformAndReturnMainDiagonal(
        const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
        const AbstractDenseSemiunitaryMatrix& denseSemiunitaryMatrix) const;
    std::unique_ptr<AbstractDenseVector> unitaryTransformAndReturnMainDiagonal(
        const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
        const std::unique_ptr<AbstractDenseSemiunitaryMatrix>& denseSemiunitaryMatrix) const;
    std::unique_ptr<AbstractDiagonalizableMatrix> unitaryTransform(
        const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
        const AbstractDenseSemiunitaryMatrix& denseSemiunitaryMatrix) const;
};
}  // namespace quantum::linear_algebra

#endif  //SPINNER_ARMALOGIC_H

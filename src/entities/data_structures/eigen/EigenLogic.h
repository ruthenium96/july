#ifndef SPINNER_EIGENLOGIC_H
#define SPINNER_EIGENLOGIC_H

#include <memory>

#include "src/entities/data_structures/AbstractDenseSemiunitaryMatrix.h"
#include "src/entities/data_structures/AbstractDenseVector.h"
#include "src/entities/data_structures/AbstractSparseSemiunitaryMatrix.h"
#include "src/entities/data_structures/AbstractSymmetricMatrix.h"

namespace quantum::linear_algebra {
template <typename T>
class EigenLogic {
  public:
    std::unique_ptr<AbstractDenseVector>
    diagonalizeValues(const AbstractDiagonalizableMatrix& symmetricMatrix) const;

    EigenCouple diagonalizeValuesVectors(const AbstractDiagonalizableMatrix& symmetricMatrix) const;

    std::unique_ptr<AbstractDenseVector> unitaryTransformAndReturnMainDiagonal(
        const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
        const AbstractDenseSemiunitaryMatrix& denseSemiunitaryMatrix) const;

    std::unique_ptr<AbstractDiagonalizableMatrix> unitaryTransform(
        const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
        const AbstractDenseSemiunitaryMatrix& denseSemiunitaryMatrix) const;
};

}  // namespace quantum::linear_algebra

#endif  //SPINNER_EIGENLOGIC_H

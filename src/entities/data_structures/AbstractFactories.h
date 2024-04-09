#ifndef SPINNER_ABSTRACTFACTORIES_H
#define SPINNER_ABSTRACTFACTORIES_H

#include "AbstractDenseSemiunitaryMatrix.h"
#include "AbstractDiagonalizableMatrix.h"
#include "AbstractSparseSemiunitaryMatrix.h"
#include "AbstractSymmetricMatrix.h"

namespace quantum::linear_algebra {

enum Precision { SINGLE, DOUBLE };

class AbstractDenseTransformAndDiagonalizeFactory {
  public:
    static std::shared_ptr<AbstractDenseTransformAndDiagonalizeFactory> defaultFactory();
    virtual std::unique_ptr<AbstractDiagonalizableMatrix>
    createDenseDiagonalizableMatrix(uint32_t size) = 0;
    virtual std::unique_ptr<AbstractDiagonalizableMatrix>
    createSparseDiagonalizableMatrix(uint32_t size) = 0;
    virtual std::unique_ptr<AbstractDenseSemiunitaryMatrix>
    createDenseSemiunitaryMatrix(uint32_t cols, uint32_t rows) = 0;
    virtual std::unique_ptr<AbstractDenseVector> createVector() = 0;

    ~AbstractDenseTransformAndDiagonalizeFactory() = default;
  private:
    Precision precision_ = Precision::DOUBLE;
  public:
    void setPrecision(Precision precision) {
        precision_ = precision;
    }
    Precision getPrecision() const {
        return precision_;
    }
};

class AbstractSparseTransformFactory {
  public:
    static std::shared_ptr<AbstractSparseTransformFactory> defaultSparseFactory();
    virtual std::unique_ptr<AbstractSparseSemiunitaryMatrix>
    createSparseSemiunitaryMatrix(uint32_t cols, uint32_t rows) = 0;
    virtual std::unique_ptr<AbstractSymmetricMatrix> createSparseSymmetricMatrix(uint32_t size) = 0;

    ~AbstractSparseTransformFactory() = default;
};

}  // namespace quantum::linear_algebra

#endif  //SPINNER_ABSTRACTFACTORIES_H

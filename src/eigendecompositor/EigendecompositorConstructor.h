#ifndef SPINNER_EIGENDECOMPOSITORCONSTRUCTOR_H
#define SPINNER_EIGENDECOMPOSITORCONSTRUCTOR_H

#include "AbstractEigendecompositor.h"
#include "src/common/runner/ConsistentModelOptimizationList.h"

namespace eigendecompositor {

class EigendecompositorConstructor {
  public:
    static std::unique_ptr<AbstractEigendecompositor> construct(
        const runner::ConsistentModelOptimizationList& consistentModelOptimizationList,
        const quantum::linear_algebra::FactoriesList& factories);

};

}  // namespace eigendecompositor

#endif  //SPINNER_EIGENDECOMPOSITORCONSTRUCTOR_H

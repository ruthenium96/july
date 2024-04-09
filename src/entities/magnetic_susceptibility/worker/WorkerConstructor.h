#ifndef SPINNER_WORKERCONSTRUCTOR_H
#define SPINNER_WORKERCONSTRUCTOR_H

#include "AbstractWorker.h"
#include "src/eigendecompositor/AbstractEigendecompositor.h"
#include "src/model/Model.h"

namespace magnetic_susceptibility::worker {

class WorkerConstructor {
  public:
    static std::unique_ptr<AbstractWorker> construct(
        const model::Model& model,
        const std::unique_ptr<eigendecompositor::AbstractEigendecompositor>& eigendecompositor,
        const quantum::linear_algebra::FactoriesList& factories);
};

}  // namespace magnetic_susceptibility::worker

#endif  //SPINNER_WORKERCONSTRUCTOR_H

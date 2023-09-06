#ifndef SPINNER_CONSISTENTMODELOPTIMIZATIONLIST_H
#define SPINNER_CONSISTENTMODELOPTIMIZATIONLIST_H

#include "src/common/physical_optimization/OptimizationList.h"
#include "src/model/Model.h"

namespace runner {
// This class keeps consistent pair (Model, OptimizationList).
// If inconsistent pair was passed to constructor, it throws a suitable exception.
class ConsistentModelOptimizationList {
  public:
    ConsistentModelOptimizationList(
        model::ModelInput,
        common::physical_optimization::OptimizationList);

    void InitializeDerivatives();
    void
    setNewValueToChangeableSymbol(const model::symbols::SymbolName& symbol_name, double new_value);

    const model::Model& getModel() const;
    const common::physical_optimization::OptimizationList& getOptimizationList() const;

  private:
    model::Model model_;
    common::physical_optimization::OptimizationList optimizationList_;
};
}  // namespace runner

#endif  //SPINNER_CONSISTENTMODELOPTIMIZATIONLIST_H

#include "OptimizedSpaceConstructor.h"

#include "src/common/Logger.h"
#include "src/space/optimization/NonAbelianSimplifier.h"
#include "src/space/optimization/PositiveProjectionsEliminator.h"
#include "src/space/optimization/S2Transformer.h"
#include "src/space/optimization/Symmetrizer.h"
#include "src/space/optimization/TzSorter.h"
#include "src/spin_algebra/GroupAdapter.h"
#include "src/spin_algebra/OrderOfSummation.h"

namespace {
std::vector<size_t> sizes_of_blocks(const space::Space& space) {
    std::vector<size_t> answer(space.getBlocks().size());
    for (int i = 0; i < space.getBlocks().size(); ++i) {
        answer[i] = space.getBlocks()[i].decomposition->size_cols();
    }
    return answer;
}
}

namespace space::optimization {

// TODO: if Space/Subspace will be transformed into classes,
//  it should be static constructor of Space.
Space OptimizedSpaceConstructor::construct(
    const runner::ConsistentModelOptimizationList& consistentModelOptimizationList,
    const quantum::linear_algebra::FactoriesList& factories) {
    const lexicographic::IndexConverter& indexConverter =
        consistentModelOptimizationList.getModel().getIndexConverter();
    const common::physical_optimization::OptimizationList& optimizationList =
        consistentModelOptimizationList.getOptimizationList();

    Space space = Space(
        consistentModelOptimizationList.getModel().getIndexConverter().get_total_space_size(),
        factories);

    bool spaceIsNormalized = true;

    if (optimizationList.isTzSorted()) {
        TzSorter tz_sorter(indexConverter, factories);
        common::Logger::detailed_msg("Tz-sortation has started.");
        space = tz_sorter.apply(std::move(space));
        common::Logger::verbose("Sizes of blocks:\n{}", fmt::join(sizes_of_blocks(space), ", "));
        common::Logger::detailed_msg("Tz-sortation is finished.");
    }

    common::Logger::separate(1, common::PrintLevel::detailed);

    if (optimizationList.isPositiveProjectionsEliminated()) {
        uint32_t max_ntz_proj = indexConverter.get_max_ntz_proj();

        PositiveProjectionsEliminator positiveProjectionsEliminator(max_ntz_proj);
        common::Logger::detailed_msg("Positive projections elimination has started.");
        space = positiveProjectionsEliminator.apply(std::move(space));
        common::Logger::verbose("Sizes of blocks:\n{}", fmt::join(sizes_of_blocks(space), ", "));
        common::Logger::detailed_msg("Positive projections elimination is finished.");
        common::Logger::separate(1, common::detailed);
    }

    for (size_t i = 0; i < optimizationList.getGroupsToApply().size(); ++i) {
        // Symmetrization breaks normalization, because of using integer values of coefficients.
        spaceIsNormalized = false;

        const auto& group = optimizationList.getGroupsToApply().at(i);
        Symmetrizer symmetrizer(indexConverter, group, factories);
        common::Logger::detailed_msg("Symmetrization has started.");
        space = symmetrizer.apply(std::move(space));
        common::Logger::verbose("Sizes of blocks:\n{}", fmt::join(sizes_of_blocks(space), ", "));
        common::Logger::detailed_msg("Symmetrization is finished.");
        if (i + 1 == optimizationList.getGroupsToApply().size()) {
            common::Logger::separate(1, common::detailed);
        } else {
            common::Logger::separate(2, common::verbose);
        }
    }

    if (!spaceIsNormalized) {
        for (auto& subspace : space.getBlocks()) {
            subspace.decomposition->normalize();
        }
    }

    if (optimizationList.isSSquaredTransformed()) {
        const auto number_of_mults = indexConverter.get_mults().size();
        auto group_adapter =
            spin_algebra::GroupAdapter(optimizationList.getGroupsToApply(), number_of_mults);

        S2Transformer transformer(
            indexConverter,
            factories,
            group_adapter.getOrderOfSummations(),
            group_adapter.getRepresentationMultiplier());

        common::Logger::detailed_msg("S2-transformation has started.");
        space = transformer.apply(std::move(space));
        common::Logger::detailed_msg("S2-transformation is finished.");
    }

    common::Logger::separate(0, common::detailed);

    return space;
}
}  // namespace space::optimization
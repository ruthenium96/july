#ifndef SPINNER_SUBSPACE_H
#define SPINNER_SUBSPACE_H

#include "src/entities/BlockProperties.h"
#include "src/entities/data_structures/AbstractDenseSemiunitaryMatrix.h"
#include "src/entities/data_structures/AbstractSparseSemiunitaryMatrix.h"

namespace space {
struct Subspace {
    BlockProperties properties;
    std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix> decomposition;
    std::optional<std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
        dense_semiunitary_matrix;

    explicit Subspace(
        std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>&&,
        std::optional<std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>&& =
            std::nullopt);
};
}  // namespace space

#endif  // SPINNER_SUBSPACE_H

#include "Subspace.h"

namespace space {
Subspace::Subspace(
    std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>&&
        new_basis_decomposition,
    std::optional<std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>&&
        dense_semiunitary_matrix) :
    decomposition(std::move(new_basis_decomposition)),
    dense_semiunitary_matrix(std::move(dense_semiunitary_matrix)) {}
}  // namespace space

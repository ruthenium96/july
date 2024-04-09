#include "Matrix.h"

Matrix::Matrix(std::vector<Submatrix>&& m) : blocks(std::move(m)) {}

Matrix::Matrix(
    const space::Space& space,
    const model::operators::Operator& new_operator,
    const lexicographic::IndexConverter& converter,
    const quantum::linear_algebra::FactoriesList& factories) {
    blocks.reserve(space.getBlocks().size());

    for (const auto& subspace : space.getBlocks()) {
        blocks.emplace_back(subspace, new_operator, converter, factories);
    }
}

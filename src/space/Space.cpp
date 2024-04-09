#include "Space.h"
namespace space {
Space::Space(uint32_t total_space_size, const quantum::linear_algebra::FactoriesList& factories) {
    auto identity_decomposition =
        factories.createSparseSemiunitaryMatrix(total_space_size, total_space_size);
    for (uint32_t lex = 0; lex < total_space_size; ++lex) {
        identity_decomposition->add_to_position(1.0, lex, lex);
    }
    blocks_.emplace_back(std::move(identity_decomposition));
}

Space::Space(std::vector<Subspace>&& v) {
    for (auto& subspace : v) {
        if (!subspace.decomposition->empty()) {
            blocks_.emplace_back(std::move(subspace));
        }
    }
}

std::vector<Subspace>& Space::getBlocks() {
    return blocks_;
}

const std::vector<Subspace>& Space::getBlocks() const {
    return blocks_;
}

}  // namespace space
#ifndef SPINNER_SPACE_H
#define SPINNER_SPACE_H

#include <deque>
#include <iostream>
#include <map>

#include "Subspace.h"
#include "src/common/lexicographic/IndexConverter.h"
#include "src/entities/data_structures/FactoriesList.h"

namespace space {
// Space is responsible for _sparse_ unitary transformation of Hamiltonian matrix.
// This transformation leads to block-diagonal Hamiltonian matrix.
class Space {
  public:
    explicit Space(
        uint32_t total_space_size,
        const quantum::linear_algebra::FactoriesList& factories);
    explicit Space(std::vector<Subspace>&& m);

    std::vector<Subspace>& getBlocks();
    const std::vector<Subspace>& getBlocks() const;

  private:
    std::vector<Subspace> blocks_;
};
}  // namespace space

#endif  // SPINNER_SPACE_H

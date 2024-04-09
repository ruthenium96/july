#ifndef SPINNER_BLOCKPROPERTIES_H
#define SPINNER_BLOCKPROPERTIES_H

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

struct BlockProperties {
    std::string get_representation_name() const;

    std::optional<uint32_t> n_proj = std::nullopt;
    std::optional<uint32_t> total_mult = std::nullopt;
    uint32_t dimensionality = 1;
    // TODO: double instead int for cases with TzSorter and without PositiveProjectionsEliminator:
    double degeneracy = 1;
    std::vector<uint8_t> representation;

    bool operator==(const BlockProperties&) const = default;
};

#endif  //SPINNER_BLOCKPROPERTIES_H

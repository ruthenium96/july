#include "Subspectrum.h"

#include <utility>

Subspectrum::Subspectrum(
    std::unique_ptr<quantum::linear_algebra::AbstractDenseVector> raw_data_,
    BlockProperties properties_) {
    raw_data = std::move(raw_data_);
    properties = std::move(properties_);
}

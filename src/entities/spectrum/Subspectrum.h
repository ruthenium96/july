#ifndef SPINNER_SUBSPECTRUM_H
#define SPINNER_SUBSPECTRUM_H

#include "src/entities/BlockProperties.h"
#include "src/entities/data_structures/AbstractDenseVector.h"
#include "src/entities/matrix/Submatrix.h"

struct Subspectrum {
    BlockProperties properties;
    std::unique_ptr<quantum::linear_algebra::AbstractDenseVector> raw_data;

    Subspectrum() = default;

    Subspectrum(
        std::unique_ptr<quantum::linear_algebra::AbstractDenseVector> raw_data_,
        BlockProperties properties_);
};

#endif  //SPINNER_SUBSPECTRUM_H

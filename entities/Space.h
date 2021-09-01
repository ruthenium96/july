#ifndef JULY_SPACE_H
#define JULY_SPACE_H

#include "Entity.h"
#include "Subspace.h"
#include "common/LexicographicIndexWorker.h"
#include <deque>
#include <iostream>
#include <map>

struct Space : public entities::Entity {

    explicit Space(const spaces::LexicographicIndexWorker& indexes);
    friend std::ostream& operator<<(std::ostream& os, const Space& space);

    std::deque<Subspace> blocks;
};

#endif // JULY_SPACE_H
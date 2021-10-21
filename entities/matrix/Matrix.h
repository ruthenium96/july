#ifndef JULY_MATRIX_H
#define JULY_MATRIX_H

#include "Submatrix.h"

#include <iostream>
#include <vector>


struct Matrix {
    std::vector<Submatrix> blocks;
    explicit Matrix(std::vector<Submatrix>&& m);
};

std::ostream &operator<<(std::ostream &os, const Matrix &space);


#endif //JULY_MATRIX_H
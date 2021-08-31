#include "Space.h"

void Space::print() {
    for (Subspace& Ss : blocks) {
        for (auto& m: Ss.basis) {
            for (auto& d: m) {
                std::cout << d.second << "*[" << d.first << "] ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    std::cout << "------" << std::endl;
}

Space::Space(const Spaces::Indexes& indexes) {

    blocks.emplace_back();
    Subspace& lex_block = blocks[0];
    lex_block.basis.resize(indexes.total_space_size);
    for (Lex_Index lex = 0; lex < indexes.total_space_size; ++lex) {
        lex_block.basis[lex][lex] = 1.0;
    }
}

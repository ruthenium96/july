#include "LexicographicIndexWorker.h"

//TODO using namesapace spaces;

spaces::LexicographicIndexWorker::LexicographicIndexWorker(std::vector<int> mults) : mults_(mults) {
    // auxiliary vector for uint32_t <-> std::vector<uint8_t> transformation
    cumulative_product.resize(mults.size() + 1);
    cumulative_product.back() = 1;
    std::partial_sum(mults.rbegin(), mults.rend(),
                     cumulative_product.rbegin() + 1, std::multiplies<>());

    // It is the size of our projection-based space.
    total_space_size = cumulative_product[0];
}

uint8_t spaces::LexicographicIndexWorker::lex_to_ntzproj(uint32_t lex) const {
    uint8_t ntz_proj = 0;
    for (int i = 0; i < mults_.size(); ++i) {
        ntz_proj += (lex % cumulative_product[i]) / cumulative_product[i + 1];
    }
    return ntz_proj;
}

std::vector<uint8_t> spaces::LexicographicIndexWorker::lex_to_nzs(uint32_t lex) const {
    std::vector<uint8_t> nzs(mults_.size());
    for (int i = 0; i < mults_.size(); ++i) {
        nzs[i] = (lex % cumulative_product[i]) / cumulative_product[i + 1];
    }
    return(std::move(nzs));
}

uint32_t spaces::LexicographicIndexWorker::nzs_to_lex(const std::vector<uint8_t> &nzs) const {
    uint32_t lex = 0;
    for (int i = 0; i < mults_.size(); ++i) {
        lex += nzs[i] * cumulative_product[i + 1];
    }
    return lex;
}
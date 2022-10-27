#include "SzSzTwoCenterTerm.h"

#include <utility>

namespace model::operators {

SzSzTwoCenterTerm::SzSzTwoCenterTerm(
    lexicographic::IndexConverter converter,
    std::shared_ptr<const TwoDNumericalParameters<double>> parameters) :
    converter_(std::move(converter)),
    coefficients_(std::move(parameters)) {}

void SzSzTwoCenterTerm::construct(
    UnitarySparseMatrix& matrix_in_lexicografical_basis,
    uint32_t index_of_vector,
    uint32_t center_a,
    uint32_t center_b) const {
    // todo: move this two from here.
    //  this two from summation in Submatrix: \sum_{a=1}^N \sum_{b=a+1}^N
    double factor = 2 * coefficients_->at(center_a, center_b);
    uint32_t projection_of_center_a =
        converter_.convert_lex_index_to_one_sz_projection(index_of_vector, center_a);
    uint32_t projection_of_center_b =
        converter_.convert_lex_index_to_one_sz_projection(index_of_vector, center_b);

    // Saz Sbz
    double diagonal_value = (projection_of_center_a - converter_.get_spins()[center_a])
        * (projection_of_center_b - converter_.get_spins()[center_b]) * factor;
    matrix_in_lexicografical_basis.add_to_position(
        diagonal_value,
        index_of_vector,
        index_of_vector);
}

std::unique_ptr<TwoCenterTerm> SzSzTwoCenterTerm::clone() const {
    return std::make_unique<SzSzTwoCenterTerm>(converter_, coefficients_);
}

std::shared_ptr<const TwoDNumericalParameters<double>> SzSzTwoCenterTerm::get_parameters() const {
    return coefficients_;
}

}  // namespace model::operators
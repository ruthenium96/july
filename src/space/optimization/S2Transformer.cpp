#include "S2Transformer.h"

#include <utility>

#include "wignerSymbols.h"

namespace space::optimization {

S2Transformer::S2Transformer(
    lexicographic::IndexConverter converter,
    quantum::linear_algebra::FactoriesList factories) :
    converter_(std::move(converter)),
    factories_(std::move(factories)) {
    {  // TODO: it is construction of chain-order summation. Move it to another place.
        order_of_summation_ = std::make_shared<std::vector<spin_algebra::AdditionInstruction>>();

        auto number_of_mults = converter_.get_mults().size();

        for (size_t i = 0; i < number_of_mults - 1; ++i) {
            spin_algebra::AdditionInstruction instruction;
            instruction.positions_of_summands = {i, number_of_mults + i - 1};
            instruction.position_of_sum = number_of_mults + i;
            order_of_summation_->push_back(instruction);
        }
    }

    sorted_s_squared_states_ = spin_algebra::SSquaredState::addAllMultiplicitiesAndSort(
        converter_.get_mults(),
        order_of_summation_);
}

space::Space S2Transformer::apply(Space&& space) const {
    std::vector<Subspace> vector_result;
    vector_result.reserve(space.getBlocks().size());

    for (size_t i = 0; i < space.getBlocks().size(); ++i) {
        Subspace& subspace = space.getBlocks()[i];

        // It is total_proj + max_total_spin = ntz_value
        // We need 2 * |total_proj| + 1 = 2 * |ntz_value - max_total_spin| + 1 =
        // |2 * ntz_value - 2 * max_total_spin - 1 + 1| + 1 =
        // |2 * ntz_value - max_total_proj + 1| + 1
        auto ntz_value = subspace.properties.n_proj.value();
        int a = 2 * (int)ntz_value + 1 - (int)converter_.get_max_ntz_proj();
        spin_algebra::Multiplicity current_mult = std::abs(a) + 1;

        const auto& s_squared_states = sorted_s_squared_states_.at(current_mult);

        subspace.dense_semiunitary_matrix =
            constructTransformationMatrix(s_squared_states, subspace);

        subspace.properties.total_mult = current_mult;
        subspace.properties.degeneracy *= current_mult;
        // TODO: for cases with TzSorter and without PositiveProjectionsEliminator:
        if (current_mult != 1) {
            subspace.properties.degeneracy /= 2;
        }

        vector_result.emplace_back(std::move(subspace));
    }

    return space::Space(std::move(vector_result));
}

std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>
S2Transformer::constructTransformationMatrix(
    const std::vector<spin_algebra::SSquaredState>& s_squared_states,
    const Subspace& subspace) const {
    auto number_of_sz_states = subspace.decomposition->size_cols();
    auto number_of_s2_states = s_squared_states.size();
    auto transformation_matrix =
        factories_.createDenseSemiunitaryMatrix(number_of_s2_states, number_of_sz_states);

#pragma omp parallel for shared( \
    number_of_s2_states, \
    number_of_sz_states, \
    s_squared_states, \
    subspace, \
    transformation_matrix) default(none)
    for (size_t i = 0; i < number_of_s2_states; ++i) {
        const auto& s_squared_state = s_squared_states[i];
        for (size_t j = 0; j < number_of_sz_states; ++j) {
            auto lex_index = subspace.decomposition->GetNewIterator(j)->getNext().index;

            double value = total_CG_coefficient(s_squared_state, lex_index);
            transformation_matrix->add_to_position(value, i, j);
        }
    }
    std::cout << "Transformation matrix (" << number_of_s2_states << ", " << number_of_sz_states
              << ") was constructed" << std::endl;
    return transformation_matrix;
}

double S2Transformer::total_CG_coefficient(
    const spin_algebra::SSquaredState& s_squared_state,
    uint32_t lex_index) const {
    auto number_of_mults = converter_.get_mults().size();

    double c = 1;

    std::vector<double> projections;
    projections.resize(2 * number_of_mults - 1);
    for (size_t i = 0; i < number_of_mults; ++i) {
        projections[i] = converter_.convert_lex_index_to_one_sz_projection(lex_index, i)
            - converter_.get_spins()[i];
    }
    for (size_t i = 0; i < number_of_mults - 1; ++i) {
        projections[i + number_of_mults] = projections[i] + projections[i + number_of_mults - 1];
    }

    for (size_t i = 0; i < order_of_summation_->size(); ++i) {
        size_t pos_one = order_of_summation_->at(i).positions_of_summands[0];
        size_t pos_two = order_of_summation_->at(i).positions_of_summands[1];
        size_t pos_sum = order_of_summation_->at(i).position_of_sum;

        double spin_one = ((double)s_squared_state.getMultiplicity(pos_one) - 1.0) / 2.0;
        double spin_two = ((double)s_squared_state.getMultiplicity(pos_two) - 1.0) / 2.0;
        double spin_sum = ((double)s_squared_state.getMultiplicity(pos_sum) - 1.0) / 2.0;

        double proj_one = projections[pos_one];
        double proj_two = projections[pos_two];
        c *= hashed_clebsh_gordan(spin_one, spin_two, spin_sum, proj_one, proj_two);
        if (c == 0.0) {
            return c;
        }
    }

    return c;
}

double S2Transformer::hashed_clebsh_gordan(double l1, double l2, double l3, double m1, double m2) {
    // TODO: implement hashing of CG-coefficients
    return WignerSymbols::clebschGordan(l1, l2, l3, m1, m2, m1 + m2);
}
}  // namespace space::optimization
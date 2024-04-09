#include "gtest/gtest.h"
#include "src/common/runner/Runner.h"

TEST(symbolic_worker, throw_assign_the_empty_name) {
    model::symbols::SymbolicWorker symbols(2);
    EXPECT_THROW(symbols.addSymbol("", 10), std::invalid_argument);
}

TEST(symbolic_worker, throw_add_the_same_symbol_name) {
    model::symbols::SymbolicWorker symbols(2);
    symbols.addSymbol("same", 10);
    EXPECT_THROW(symbols.addSymbol("same", 10), std::invalid_argument);
}

TEST(symbolic_worker, throw_assign_the_same_exchange) {
    model::symbols::SymbolicWorker symbols(2);
    auto J = symbols.addSymbol("same", 10, model::symbols::J);
    symbols.assignSymbolToIsotropicExchange(J, 0, 1);
    EXPECT_THROW(symbols.assignSymbolToIsotropicExchange(J, 0, 1), std::invalid_argument);
}

TEST(symbolic_worker, throw_2222_gfactor_only_one_g_was_initialized) {
    std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2};

    model::ModelInput modelInput(mults);

    double J_value = 10;
    auto J = modelInput.addSymbol("J", J_value);
    modelInput.assignSymbolToIsotropicExchange(J, 0, 1)
        .assignSymbolToIsotropicExchange(J, 1, 2)
        .assignSymbolToIsotropicExchange(J, 2, 3)
        .assignSymbolToIsotropicExchange(J, 3, 0);

    double g_value = 2.0;
    auto g = modelInput.addSymbol("g", g_value);
    modelInput.assignSymbolToGFactor(g, 0);

    EXPECT_THROW(runner::Runner runner(modelInput), std::invalid_argument);
}

TEST(symbolic_worker, throw_2222_gfactor_only_one_g_was_not_initialized) {
    std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2};

    model::ModelInput modelInput(mults);

    double J_value = 10;
    auto J = modelInput.addSymbol("J", J_value);
    modelInput.assignSymbolToIsotropicExchange(J, 0, 1)
        .assignSymbolToIsotropicExchange(J, 1, 2)
        .assignSymbolToIsotropicExchange(J, 2, 3)
        .assignSymbolToIsotropicExchange(J, 3, 0);

    double g_value = 2.0;
    auto g = modelInput.addSymbol("g", g_value);
    modelInput.assignSymbolToGFactor(g, 0)
        .assignSymbolToGFactor(g, 1)
        .assignSymbolToGFactor(g, 2);

    EXPECT_THROW(runner::Runner runner(modelInput), std::invalid_argument);
}

TEST(symbolic_worker, throw_set_new_value_to_unchangeable_symbol) {
    model::symbols::SymbolicWorker symbols_(10);

    auto unchangeable =
        symbols_.addSymbol("Unchangeable", NAN, false, std::nullopt);

    EXPECT_THROW(
        symbols_.setNewValueToChangeableSymbol(unchangeable, INFINITY),
        std::invalid_argument);
}

TEST(symbolic_worker, throw_specified_not_as_J_symbol) {
    size_t number_of_spins = 10;
    model::symbols::SymbolicWorker symbols_(number_of_spins);

    double gChangeable = 2;

    auto not_specified =
        symbols_.addSymbol("not_specified", NAN, true, std::nullopt);
    symbols_.assignSymbolToIsotropicExchange(not_specified, 1, 3);
    auto g_changeable =
        symbols_.addSymbol("gChangeable", gChangeable, true, model::symbols::g_factor);
    EXPECT_THROW(
        symbols_.assignSymbolToIsotropicExchange(g_changeable, 2, 7),
        std::invalid_argument);
}

TEST(symbolic_worker, throw_specified_not_as_g_symbol) {
    size_t number_of_spins = 10;
    model::symbols::SymbolicWorker symbols_(number_of_spins);

    double JChangeable = 10;

    auto not_specified =
        symbols_.addSymbol("not_specified", NAN, true, std::nullopt);
    symbols_.assignSymbolToGFactor(not_specified, 1);
    auto J_changeable = symbols_.addSymbol("JChangeable", JChangeable, true, model::symbols::J);
    EXPECT_THROW(symbols_.assignSymbolToGFactor(J_changeable, 2), std::invalid_argument);
}

TEST(numerical_worker, set_new_value_to_changeable_J_g) {
    size_t number_of_spins = 10;
    model::symbols::SymbolicWorker symbols_(number_of_spins);

    double JChangeable_value = 10;
    double gChangeable_value = 2;

    auto J_changeable =
        symbols_.addSymbol("JChangeable", JChangeable_value, true, model::symbols::J);
    auto g_changeable =
        symbols_.addSymbol("gChangeable", gChangeable_value, true, model::symbols::g_factor);
    symbols_.assignSymbolToIsotropicExchange(J_changeable, 2, 7);
    for (size_t i = 0; i < number_of_spins; ++i) {
        symbols_.assignSymbolToGFactor(g_changeable, i);
    }

    auto parametersHandler = model::symbols::NumericalWorker(symbols_, number_of_spins);

    auto shared_ptr_J = parametersHandler.getIsotropicExchangeParameters();
    auto shared_ptr_g = parametersHandler.getGFactorParameters();
    EXPECT_EQ(JChangeable_value, shared_ptr_J->at(2, 7));
    EXPECT_EQ(gChangeable_value, shared_ptr_g->at(7));
    parametersHandler.setNewValueToChangeableSymbol(J_changeable, 2 * JChangeable_value);
    EXPECT_EQ(2 * JChangeable_value, shared_ptr_J->at(2, 7));
    EXPECT_EQ(gChangeable_value, shared_ptr_g->at(7));
    parametersHandler.setNewValueToChangeableSymbol(g_changeable, 2 * gChangeable_value);
    EXPECT_EQ(2 * JChangeable_value, shared_ptr_J->at(2, 7));
    EXPECT_EQ(2 * gChangeable_value, shared_ptr_g->at(7));
}
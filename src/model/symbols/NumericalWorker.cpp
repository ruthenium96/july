#include "NumericalWorker.h"

#include <cmath>
#include <utility>

namespace model::symbols {
NumericalWorker::NumericalWorker(SymbolicWorker symbolicWorker, size_t number_of_spins) :
    symbolicWorker_(std::move(symbolicWorker)),
    number_of_spins_(number_of_spins) {
    if (getSymbolicWorker().isIsotropicExchangeInitialized()) {
        updateIsotropicExchangeParameters();
    }
    if (getSymbolicWorker().isZFSInitialized()) {
        updateZFSParameters();
    }
    if (getSymbolicWorker().isGFactorInitialized()) {
        updateGFactorParameters();
    }
    if (getSymbolicWorker().isThetaInitialized()) {
        updateThetaParameter();
    }
}

const SymbolicWorker& NumericalWorker::getSymbolicWorker() const {
    return symbolicWorker_;
}

void NumericalWorker::updateIsotropicExchangeParameters() {
    if (numeric_isotropic_exchanges_ == nullptr) {
        numeric_isotropic_exchanges_ =
            std::make_shared<TwoDNumericalParameters<double>>(number_of_spins_, NAN);
    }

    for (size_t i = 0; i < number_of_spins_; ++i) {
        for (size_t j = 0; j < number_of_spins_; ++j) {
            double value;
            auto mb_symbol_name = getSymbolicWorker().getIsotropicExchangeSymbolName(i, j);
            if (!mb_symbol_name.has_value()) {
                value = NAN;
            } else {
                value = getSymbolicWorker().getValueOfName(mb_symbol_name.value());
            }
            numeric_isotropic_exchanges_->at(i, j) = value;
        }
    }
}

void NumericalWorker::updateGFactorParameters() {
    if (numeric_g_factors_ == nullptr) {
        numeric_g_factors_ = std::make_shared<OneDNumericalParameters<double>>(number_of_spins_, 0);
    }
    if (numeric_g_g_.first == nullptr) {
        numeric_g_g_.first = std::make_shared<OneDNumericalParameters<double>>(number_of_spins_, 0);
        numeric_g_g_.second =
            std::make_shared<TwoDNumericalParameters<double>>(number_of_spins_, 0);
    }

    for (size_t i = 0; i < number_of_spins_; ++i) {
        auto symbol_name = getSymbolicWorker().getGFactorSymbolName(i);
        double value = getSymbolicWorker().getValueOfName(symbol_name);
        numeric_g_factors_->at(i) = value;
    }

    for (size_t i = 0; i < number_of_spins_; ++i) {
        auto symbol_name = getSymbolicWorker().getGFactorSymbolName(i);
        double value = getSymbolicWorker().getValueOfName(symbol_name);
        numeric_g_g_.first->at(i) = value * value;
    }

    for (size_t i = 0; i < number_of_spins_; ++i) {
        for (size_t j = 0; j < number_of_spins_; ++j) {
            if (i == j) {
                continue;
            }
            auto symbol_name_i = getSymbolicWorker().getGFactorSymbolName(i);
            auto symbol_name_j = getSymbolicWorker().getGFactorSymbolName(j);
            double value = getSymbolicWorker().getValueOfName(symbol_name_i)
                * getSymbolicWorker().getValueOfName(symbol_name_j);
            ;
            numeric_g_g_.second->at(i, j) = value;
        }
    }

    for (auto& [symbol_name, _] : numeric_g_g_derivatives_) {
        updateGGFactorDerivativesParameters(symbol_name);
    }
}

void NumericalWorker::updateZFSParameters() {
    if (numeric_ZFS_.first == nullptr) {
        numeric_ZFS_.first =
            std::make_shared<OneDNumericalParameters<double>>(number_of_spins_, NAN);
    }

    for (size_t i = 0; i < number_of_spins_; ++i) {
        auto mb_symbol = getSymbolicWorker().getZFSSymbolNames(i);
        if (mb_symbol.has_value()) {
            double D_value = getSymbolicWorker().getValueOfName(mb_symbol.value().D);
            numeric_ZFS_.first->at(i) = D_value;
        } else {
            numeric_ZFS_.first->at(i) = NAN;
        }
    }
}

void NumericalWorker::updateThetaParameter() {
    if (numeric_Theta_ == nullptr) {
        numeric_Theta_ = std::make_shared<double>(NAN);
    }

    auto symbol_name = getSymbolicWorker().getThetaSymbolName().value();
    double value = getSymbolicWorker().getValueOfName(symbol_name);
    *numeric_Theta_ = value;
}

void NumericalWorker::updateGGFactorDerivativesParameters(const SymbolName& symbol_name) {
    auto pair_of_pointers = numeric_g_g_derivatives_[symbol_name];
    for (size_t i = 0; i < number_of_spins_; ++i) {
        if (getSymbolicWorker().getGFactorSymbolName(i) == symbol_name) {
            pair_of_pointers.first->at(i) = 2 * numeric_g_factors_->at(i);
        }
    }

    for (size_t i = 0; i < number_of_spins_; ++i) {
        for (size_t j = 0; j < number_of_spins_; ++j) {
            if (i == j) {
                continue;
            }
            double value = 0;
            if (getSymbolicWorker().getGFactorSymbolName(i) == symbol_name) {
                value += numeric_g_factors_->at(j);
            }
            if (getSymbolicWorker().getGFactorSymbolName(j) == symbol_name) {
                value += numeric_g_factors_->at(i);
            }
            pair_of_pointers.second->at(i, j) = value;
        }
    }
}

std::shared_ptr<const TwoDNumericalParameters<double>>
NumericalWorker::getIsotropicExchangeParameters() const {
    // TODO: refactor it:
    if (numeric_isotropic_exchanges_ == nullptr) {
        throw std::invalid_argument("Isotropic exchange interaction was not initialized");
    }
    return numeric_isotropic_exchanges_;
}

std::pair<
    std::shared_ptr<OneDNumericalParameters<double>>,
    std::optional<std::shared_ptr<OneDNumericalParameters<double>>>>
NumericalWorker::getZFSParameters() const {
    if (numeric_ZFS_.first == nullptr) {
        throw std::invalid_argument("Zero field splitting was not initialized");
    }
    return numeric_ZFS_;
}

std::shared_ptr<const double> NumericalWorker::getThetaParameter() const {
    if (numeric_Theta_ == nullptr) {
        throw std::invalid_argument("Theta was not initialized");
    }
    return numeric_Theta_;
}

std::shared_ptr<const OneDNumericalParameters<double>>
NumericalWorker::getGFactorParameters() const {
    if (numeric_g_factors_ == nullptr) {
        throw std::invalid_argument("g-factors were not initialized");
    }
    return numeric_g_factors_;
}

std::pair<
    std::shared_ptr<const OneDNumericalParameters<double>>,
    std::shared_ptr<const TwoDNumericalParameters<double>>>
NumericalWorker::getGGParameters() const {
    if (numeric_g_g_.first == nullptr) {
        throw std::invalid_argument("g-factors were not initialized");
    }
    return numeric_g_g_;
}

std::shared_ptr<const TwoDNumericalParameters<double>>
NumericalWorker::constructIsotropicExchangeDerivativeParameters(const SymbolName& symbol_name) {
    if (getSymbolicWorker().getSymbolProperty(symbol_name).type_enum != SymbolTypeEnum::J) {
        throw std::invalid_argument(
            symbol_name.get_name() + " has been specified as not J parameter");
    }

    auto ptr_to_derivative =
        std::make_shared<TwoDNumericalParameters<double>>(number_of_spins_, NAN);

    for (size_t i = 0; i < number_of_spins_; ++i) {
        for (size_t j = 0; j < number_of_spins_; ++j) {
            double value;
            auto mb_symbol_name = getSymbolicWorker().getIsotropicExchangeSymbolName(i, j);
            if (mb_symbol_name.has_value() && mb_symbol_name.value() == symbol_name) {
                value = 1;
            } else {
                value = NAN;
            }
            ptr_to_derivative->at(i, j) = value;
        }
    }

    // TODO: Is it just storing of this pointer?
    numeric_isotropic_exchange_derivatives_.push_back(ptr_to_derivative);

    return ptr_to_derivative;
}

std::shared_ptr<const OneDNumericalParameters<double>>
NumericalWorker::constructZFSDerivativeParameters(const SymbolName& symbol_name) {
    if (getSymbolicWorker().getSymbolProperty(symbol_name).type_enum != SymbolTypeEnum::D) {
        throw std::invalid_argument(
            symbol_name.get_name() + " has been specified as not D parameter");
    }

    auto ptr_to_derivative =
        std::make_shared<OneDNumericalParameters<double>>(number_of_spins_, NAN);

    for (size_t i = 0; i < number_of_spins_; ++i) {
        double value;
        auto mb_symbol_name = getSymbolicWorker().getZFSSymbolNames(i);
        if (mb_symbol_name.has_value() && mb_symbol_name.value().D == symbol_name) {
            value = 1;
        } else {
            value = NAN;
        }
        ptr_to_derivative->at(i) = value;
    }

    return ptr_to_derivative;
}

std::pair<
    std::shared_ptr<const OneDNumericalParameters<double>>,
    std::shared_ptr<const TwoDNumericalParameters<double>>>
NumericalWorker::constructGGDerivativeParameters(const SymbolName& symbol_name) {
    if (getSymbolicWorker().getSymbolProperty(symbol_name).type_enum != SymbolTypeEnum::g_factor) {
        throw std::invalid_argument(
            symbol_name.get_name() + " has been specified as not g factor parameter");
    }

    if (!getSymbolicWorker().isGFactorInitialized()) {
        throw std::length_error("G factor parameters have not been initialized");
    }

    if (numeric_g_g_derivatives_.find(symbol_name) != numeric_g_g_derivatives_.end()) {
        throw std::invalid_argument(
            "Derivatives from GG for" + symbol_name.get_name() + " have been already constructed");
    }

    numeric_g_g_derivatives_[symbol_name] = {
        std::make_shared<OneDNumericalParameters<double>>(number_of_spins_, 0),
        std::make_shared<TwoDNumericalParameters<double>>(number_of_spins_, 0)};

    updateGGFactorDerivativesParameters(symbol_name);

    return numeric_g_g_derivatives_[symbol_name];
}

void NumericalWorker::setNewValueToChangeableSymbol(
    const SymbolName& symbol_name,
    double new_value) {
    symbolicWorker_.setNewValueToChangeableSymbol(symbol_name, new_value);

    auto symbol_data = getSymbolicWorker().getSymbolProperty(symbol_name);

    if (symbol_data.type_enum == symbols::SymbolTypeEnum::J) {
        updateIsotropicExchangeParameters();
        // NB: numeric_isotropic_exchange_derivatives_ does not change when J changes.
    } else if (symbol_data.type_enum == symbols::SymbolTypeEnum::g_factor) {
        updateGFactorParameters();
    } else if (symbol_data.type_enum == symbols::SymbolTypeEnum::Theta) {
        updateThetaParameter();
    } else if (symbol_data.type_enum == symbols::SymbolTypeEnum::D) {
        updateZFSParameters();
    }
    // and other things
}

}  // namespace model::symbols
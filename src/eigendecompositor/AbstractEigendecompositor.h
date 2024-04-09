#ifndef SPINNER_ABSTRACTEIGENDECOMPOSITOR_H
#define SPINNER_ABSTRACTEIGENDECOMPOSITOR_H

#include <optional>

#include "src/model/Model.h"

namespace eigendecompositor {
class AbstractEigendecompositor {
  public:
    void BuildSpectra(
        const std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
            operators,
        const std::map<
            std::pair<common::QuantityEnum, model::symbols::SymbolName>,
            std::shared_ptr<const model::operators::Operator>>& derivatives_operators,
        const space::Space& space);
    bool BuildSpectraWasCalled() const;

    // What means Spectrum's optional?
    // The getter *may* return std::nullopt if corresponding Operator has been passed into initialize.
    // What means Matrix's optional?
    // The getter returns some value if underlying wrapper actually store Matrix.
    virtual std::optional<std::reference_wrapper<const Spectrum>>
        getSpectrum(common::QuantityEnum) const = 0;
    virtual std::optional<std::reference_wrapper<const Matrix>>
        getMatrix(common::QuantityEnum) const = 0;
    virtual std::optional<std::reference_wrapper<const Spectrum>>
    getSpectrumDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const = 0;
    virtual std::optional<std::reference_wrapper<const Matrix>>
    getMatrixDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const = 0;

    virtual ~AbstractEigendecompositor() = default;

    virtual void initialize(
        std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
            operators_to_calculate,
        std::map<
            std::pair<common::QuantityEnum, model::symbols::SymbolName>,
            std::shared_ptr<const model::operators::Operator>>& derivatives_operators_to_calculate,
        uint32_t number_of_subspaces) = 0;
    virtual std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
    BuildSubspectra(size_t number_of_block, const space::Subspace& subspace) = 0;
    virtual void finalize() = 0;
  private:
    bool buildSpectraWasCalled = false;
};

}  // namespace eigendecompositor

#endif  //SPINNER_ABSTRACTEIGENDECOMPOSITOR_H

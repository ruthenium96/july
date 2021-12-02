#ifndef JULY_UNIQUEGONLYSSQUAREDMUSQUAREDWORKER_H
#define JULY_UNIQUEGONLYSSQUAREDMUSQUAREDWORKER_H

#include "MuSquaredWorker.h"

namespace magnetic_susceptibility {

// mu^2 = mu_B^2 * g_{iso}^2 * <S^2_{total}>
class UniqueGOnlySSquaredMuSquaredWorker: public MuSquaredWorker {
  public:
    UniqueGOnlySSquaredMuSquaredWorker(
        DenseVector&& energy,
        DenseVector&& degeneracy,
        DenseVector&& s_squared,
        double g_unique);

    virtual double theory_at_temperature(double temperature) const override;

    std::vector<ValueAtTemperature> calculateDerivative(
        symbols::SymbolTypeEnum symbol_type,
        DenseVector&& derivative_value) const override;

  private:
    double g_unique_;
    DenseVector s_squared_;
};

}  // namespace magnetic_susceptibility

#endif  //JULY_UNIQUEGONLYSSQUAREDMUSQUAREDWORKER_H

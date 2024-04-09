#ifndef SPINNER_PRINTINGFUNCTIONS_H
#define SPINNER_PRINTINGFUNCTIONS_H

#include <fmt/ostream.h>

#include <ostream>

#include "Quantity.h"
#include "src/eigendecompositor/AbstractEigendecompositor.h"
#include "src/entities/BlockProperties.h"
#include "src/entities/magnetic_susceptibility/assistant/ExperimentalValuesWorker.h"
#include "src/entities/matrix/Matrix.h"
#include "src/entities/matrix/Submatrix.h"
#include "src/entities/spectrum/Spectrum.h"
#include "src/entities/spectrum/Subspectrum.h"
#include "src/model/symbols/SymbolName.h"
#include "src/space/Space.h"
#include "src/space/Subspace.h"

std::ostream& operator<<(std::ostream& os, const space::Space& space);
std::ostream& operator<<(std::ostream& os, const space::Subspace& subspace);
std::ostream& operator<<(std::ostream& os, const Spectrum& spectrum);
std::ostream& operator<<(std::ostream& os, const Subspectrum& subspectrum);
template <> struct fmt::formatter<Spectrum> : fmt::ostream_formatter {};
std::ostream& operator<<(std::ostream& os, const Matrix& matrix);
std::ostream& operator<<(std::ostream& os, const Submatrix& submatrix);
std::ostream& operator<<(std::ostream& os, const BlockProperties& properties);
std::ostream& operator<<(std::ostream& os, const common::QuantityEnum& quantity);

namespace common {
void inputPrint(const std::string& input_string);
void preRegressionPrint(
    const std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>& quantities,
    const std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>& derivatives);
void postRegressionPrint(
    const std::vector<model::symbols::SymbolName>& changeable_names,
    const std::vector<double>& changeable_values,
    double rss);
void stepOfRegressionStartPrint(const std::vector<model::symbols::SymbolName>& changeable_names,
                               const std::vector<double>& changeable_values);
void stepOfRegressionFinishPrint(double loss);
void initialExperimentalValuesPrint(
    const std::vector<magnetic_susceptibility::ValueAtTemperature>& experimental_values,
    magnetic_susceptibility::ExperimentalValuesEnum experimental_values_type);
void experimentalValuesPrint(
    const std::vector<magnetic_susceptibility::ValueAtTemperature>& experimental_mu_squared,
    const std::vector<double>& weights);
void theoreticalValuesPrint(
    const std::vector<magnetic_susceptibility::ValueAtTemperature>& theor_values);
}

#endif  //SPINNER_PRINTINGFUNCTIONS_H

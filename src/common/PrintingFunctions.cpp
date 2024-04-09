#include "PrintingFunctions.h"

#include "src/common/Logger.h"

#include <magic_enum.hpp>

std::ostream& operator<<(std::ostream& os, const space::Space& space) {
    for (const space::Subspace& subspace : space.getBlocks()) {
        os << subspace;
    }
    os << "------" << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const space::Subspace& subspace) {
    os << subspace.properties;
    subspace.decomposition->print(os);
    if (subspace.dense_semiunitary_matrix.has_value()) {
        os << std::endl;
        subspace.dense_semiunitary_matrix.value()->print(os);
    }
    os << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const Spectrum& spectrum) {
    for (const Subspectrum& subspectrum : spectrum.blocks) {
        os << subspectrum;
    }
    os << "------" << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const Subspectrum& subspectrum) {
    os << subspectrum.properties << std::endl;
    subspectrum.raw_data->print(os);
    os << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
    for (const Submatrix& submatrix : matrix.blocks) {
        os << submatrix;
    }
    os << "------" << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const Submatrix& submatrix) {
    os << submatrix.properties;
    submatrix.raw_data->print(os);
    os << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const BlockProperties& properties) {
    os << "Total n-projection: ";
    if (properties.n_proj.has_value()) {
        os << properties.n_proj.value();
    } else {
        os << "none";
    }
    os << "\nTotal multiplicity: ";
    if (properties.total_mult.has_value()) {
        os << properties.total_mult.value();
    } else {
        os << "none";
    }
    os << '\n'
       << "    dimensionality: " << properties.dimensionality << '\n'
       << "        degeneracy: " << properties.degeneracy << '\n'
       << "    representation: " << properties.get_representation_name() << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const common::QuantityEnum& quantity_enum) {
    os << "Quantity type: " << magic_enum::enum_name(quantity_enum) << std::endl;
    return os;
}

namespace common {
void inputPrint(const std::string& input_string) {
    common::Logger::verbose_msg("Input file:");
    common::Logger::verbose("{}", input_string);
    common::Logger::separate(0, common::verbose);
}

void preRegressionPrint(
    const std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>& quantities,
    const std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>& derivatives) {
    common::Logger::basic_msg("Regression has started.");
    common::Logger::separate(1, common::basic);

    common::Logger::verbose_msg("Quantities for explicit calculation:");

    std::string quantities_names;
    for (const auto& [value, _] : quantities) {
        quantities_names += magic_enum::enum_name(value);
        quantities_names += ", ";
    }
    common::Logger::verbose("Quantities: {}", quantities_names);

    if (!derivatives.empty()) {
        std::string derivative_names;
        for (const auto& [pair, _] : derivatives) {
            derivative_names += "d(";
            derivative_names += magic_enum::enum_name(pair.first);
            derivative_names += ")/d" + pair.second.get_name() + ", ";
        }
        common::Logger::verbose("Derivatives: {}", derivative_names);
    }

    common::Logger::separate(1, common::verbose);
}

void postRegressionPrint(
    const std::vector<model::symbols::SymbolName>& changeable_names,
    const std::vector<double>& changeable_values,
    double rss) {
    common::Logger::basic_msg("Regression is finished. Final values:");
    for (size_t i = 0; i < changeable_names.size(); ++i) {
        common::Logger::basic("{}: {}", changeable_names[i].get_name(), changeable_values[i]);
    }
    common::Logger::separate(2, common::basic);
    common::Logger::basic("Loss function = {}", rss);
    common::Logger::separate(0, common::basic);
}

void stepOfRegressionStartPrint(
    const std::vector<model::symbols::SymbolName>& changeable_names,
    const std::vector<double>& changeable_values) {
    common::Logger::verbose_msg("New step of regression has started with values:");
    for (size_t i = 0; i < changeable_names.size(); ++i) {
        common::Logger::verbose("{}: {}", changeable_names[i].get_name(), changeable_values[i]);
    }
    common::Logger::separate(2, common::verbose);
}

void stepOfRegressionFinishPrint(double loss) {
    common::Logger::verbose("Loss function = {}", loss);
    common::Logger::separate(1, common::verbose);
}

void initialExperimentalValuesPrint(
    const std::vector<magnetic_susceptibility::ValueAtTemperature>& experimental_values,
    magnetic_susceptibility::ExperimentalValuesEnum experimental_values_type) {
    common::Logger::debug_msg("Initial experimental values in {}:",
                              magic_enum::enum_name(experimental_values_type));
    for (size_t i = 0; i < experimental_values.size(); ++i) {
        common::Logger::debug(
            "{:.8e}    {:.8e}",
            experimental_values.at(i).temperature,
            experimental_values.at(i).value);
    }
    common::Logger::separate(0, common::debug);
}

void experimentalValuesPrint(
    const std::vector<magnetic_susceptibility::ValueAtTemperature>& experimental_mu_squared,
    const std::vector<double>& weights) {
    common::Logger::verbose_msg("Experimental values, corrected by ratio and in mu-squared, and weights:");
    for (size_t i = 0; i < experimental_mu_squared.size(); ++i) {
        common::Logger::verbose(
            "{:.8e}    {:.8e}    {:.8e}",
            experimental_mu_squared.at(i).temperature,
            experimental_mu_squared.at(i).value,
            weights.at(i));
    }
    common::Logger::separate(0, common::verbose);
}

// todo: also print magnetic_susceptibility::ExperimentalValuesEnum experimental_values_type
void theoreticalValuesPrint(
    const std::vector<magnetic_susceptibility::ValueAtTemperature>& theor_values) {
    common::Logger::basic_msg("\nValue, Temperature");
    for (auto [temperature, value] : theor_values) {
        common::Logger::basic("{:.8e}    {:.8e}", temperature, value);
    }
    common::Logger::separate(0, common::basic);
}

}
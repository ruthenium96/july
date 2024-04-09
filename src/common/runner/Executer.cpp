#include "Executer.h"

#include "src/common/PrintingFunctions.h"

namespace runner {
void Executer::execute(input::Parser parser) {
    auto optimization_list = parser.getOptimizationList().value();
    auto factoriesList = parser.getFactoriesList().value();

    for (const auto& model_input : parser.getModelInputs()) {
        auto runner = runner::Runner(model_input, optimization_list, factoriesList);

        if (parser.getTemperaturesForSimulation().has_value()) {
            std::vector<magnetic_susceptibility::ValueAtTemperature> theor_values;
            for (auto temp : parser.getTemperaturesForSimulation().value()) {
                double value = runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(temp);
                theor_values.push_back({temp, value});
            }

            common::theoreticalValuesPrint(theor_values);
        } else if (parser.getNonlinearSolver().has_value()) {
            runner.initializeExperimentalValues(parser.getExperimentalValuesWorker().value());

            runner.minimizeResidualError(parser.getNonlinearSolver().value());
            auto theor_final = runner.getMagneticSusceptibilityController().getTheoreticalValues();

            common::theoreticalValuesPrint(theor_final);
        }
    }
}
}  // namespace runner
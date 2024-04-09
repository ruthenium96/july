#include "ControlParser.h"

#include "Tools.h"

#ifdef _Eigen_BUILT
    #include "src/entities/data_structures/eigen/EigenFactories.h"
#endif
#ifdef _Arma_BUILT
    #include "src/entities/data_structures/arma/ArmaFactories.h"
#endif

template<>
std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory> YAML::Node::as() const {
    auto dense_algebra_package_string = as<std::string>();

    if (dense_algebra_package_string == "default") {
        return quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory::defaultFactory();
    } else if (dense_algebra_package_string == "arma") {
#ifdef _Arma_BUILT
        return std::make_shared<quantum::linear_algebra::ArmaDenseTransformAndDiagonalizeFactory>();
#else
        throw std::invalid_argument("Arma was not found, thus cannot be used");
#endif
    } else if (dense_algebra_package_string == "eigen") {
#ifdef _Eigen_BUILT
        return std::make_shared<quantum::linear_algebra::EigenDenseTransformAndDiagonalizeFactory>();
#else
        throw std::invalid_argument("Eigen was not found, thus cannot be used");
#endif
    } else {
        throw std::invalid_argument("Incorrect control::dense_algebra_package: " +
                                    dense_algebra_package_string);
    }
}


namespace input {
ControlParser::ControlParser(YAML::Node control_node) {
    print_level_ = extractValue<common::PrintLevel>(control_node, "print_level");
    common::Logger::set_level(print_level_.value());

    constructFactoriesList(control_node);

    throw_if_node_is_not_empty(control_node);
}

void ControlParser::constructFactoriesList(YAML::Node& control_node) {
    auto densePrecision =
        extractValue<quantum::linear_algebra::Precision>(control_node, "dense_precision");

    auto denseFactory = extractValue<
        std::shared_ptr<quantum::linear_algebra::AbstractDenseTransformAndDiagonalizeFactory>
        >(control_node, "dense_algebra_package");
    denseFactory->setPrecision(densePrecision);

    auto sparseFactory =
        quantum::linear_algebra::AbstractSparseTransformFactory::defaultSparseFactory();

    factoriesList_ = quantum::linear_algebra::FactoriesList(denseFactory, sparseFactory);
}

common::PrintLevel ControlParser::getPrintLevel() const {
    return print_level_.value();
}

const std::optional<quantum::linear_algebra::FactoriesList>& ControlParser::getFactoriesList() const {
    return factoriesList_;
}
}  // namespace input
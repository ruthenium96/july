include(FetchContent)

FetchContent_Declare(
        magic_enum
        GIT_REPOSITORY https://github.com/Neargye/magic_enum.git
        GIT_TAG "e55b9b54d5cf61f8e117cafb17846d7d742dd3b4"
)
FetchContent_MakeAvailable(magic_enum)
include_directories(${magic_enum_SOURCE_DIR}/)

message("Trying to find system yaml-cpp...")
find_package(yaml-cpp)
if (yaml-cpp_FOUND)
    include_directories(${YAML_CPP_SOURCE_DIR} ${YAML_CPP_SOURCE_DIR}/include/ ${YAML_CPP_SOURCE_DIR}/yaml-cpp)
    message(" *** System yaml-cpp was found")
else ()
    include(FetchContent)
    FetchContent_Declare(
            yaml-cpp
            GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
            GIT_TAG "0.8.0"
    )
    FetchContent_MakeAvailable(yaml-cpp)
    include_directories(${YAML_CPP_SOURCE_DIR} ${YAML_CPP_SOURCE_DIR}/include/ ${YAML_CPP_SOURCE_DIR}/yaml-cpp)
    message(" *** Yaml-cpp was downloaded")
endif ()

add_library(input
        Parser.cpp Parser.h
        Tools.cpp Tools.h
        ModelInputParser.cpp ModelInputParser.h
        OptimizationsParser.cpp OptimizationsParser.h
        JobParser.cpp JobParser.h
        ControlParser.cpp ControlParser.h)

add_definitions(-D_LBFGSppAdapter_BUILT=${Eigen3_FOUND}) # necessary for using of lbfgs-cpp
add_definitions(-D_Eigen_BUILT=${Eigen3_FOUND}) # necessary for construction of FactoriesList
add_definitions(-D_Arma_BUILT=${ARMADILLO_FOUND}) # necessary for construction of FactoriesList

target_link_libraries(input yaml-cpp::yaml-cpp magic_enum model common nonlinear_solver data_structures)

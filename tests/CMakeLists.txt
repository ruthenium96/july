message("Trying to find system googletest...")
find_package(GTest)
if (GTest_FOUND)
    include_directories(${GTest_INCLUDE_DIRS})
    message(" *** System GTest was found")
else ()
    include(FetchContent)
    FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG "release-1.11.0"
    )
    FetchContent_MakeAvailable(googletest)
    include_directories(${CMAKE_SOURCE_DIR}/libs/googletest/googletest/include)
    message(" *** GTest was downloaded")
endif ()

add_library(tests_tools
        tools/MeanAndDeviation.cpp
        tools/GenerateSameMatrix.cpp)

add_executable(spinner_time_test
        time_tests/main.cpp)

add_executable(spinner_test
        main.cpp
        integration_tests/spectrum_final_equivalence_tests.cpp
        integration_tests/simple_analytical_dependencies_tests.cpp
        unit_tests/group_test.cpp
        unit_tests/consistentModelOptimizationList_tests.cpp
        unit_tests/symmetrizer_tests.cpp
        unit_tests/non_abelian_simplifier_tests.cpp
        unit_tests/interactions_test.cpp
        unit_tests/hamiltonian_operator_tests.cpp
        unit_tests/symbols_tests.cpp
        unit_tests/parser_tests.cpp
        non_hamiltonian_operators_tests.cpp
        unit_tests/magnetic_susceptibility_tests.cpp
        integration_tests/spectrum_builder_tests.cpp)

find_package(Armadillo 10.6)
find_package(Eigen3 3.3)

if (ARMADILLO_FOUND)
    target_sources(spinner_test PRIVATE unit_tests/linear_algebra/armaDoubleIndividual_tests.cpp)
    target_sources(spinner_test PRIVATE unit_tests/linear_algebra/armaFloatIndividual_tests.cpp)
endif ()
if (Eigen3_FOUND)
    target_sources(spinner_test PRIVATE unit_tests/linear_algebra/eigenDoubleIndividual_tests.cpp)
    target_sources(spinner_test PRIVATE unit_tests/linear_algebra/eigenFloatIndividual_tests.cpp)
endif ()

if (Eigen3_FOUND)
    target_sources(tests_tools PRIVATE tools/concreteSolverConstructors/create_LBFGSpp.cpp)
endif ()
target_sources(tests_tools PRIVATE
        tools/concreteSolverConstructors/create_stlbfgs.cpp
        tools/concreteSolverConstructors/create_optimNM.cpp
        )

if (Eigen3_FOUND)
    target_sources(spinner_test PRIVATE integration_tests/fitting_theoretical_curves/LBFGSpp_tests.cpp)
endif ()
target_sources(spinner_test PRIVATE
        integration_tests/fitting_theoretical_curves/stlbfgs_tests.cpp
        integration_tests/fitting_theoretical_curves/optimNM_tests.cpp
        )

if (Eigen3_FOUND)
    target_sources(spinner_test PRIVATE unit_tests/nonlinear_solver/LBFGSpp_tests.cpp)
endif ()
target_sources(spinner_test PRIVATE
        unit_tests/nonlinear_solver/stlbfgs_tests.cpp
        unit_tests/nonlinear_solver/optimNM_tests.cpp
        )

if (ARMADILLO_FOUND)
    target_sources(spinner_time_test PRIVATE
            time_tests/PerformanceTest.cpp
            )
    target_sources(tests_tools PRIVATE
            tools/AllSparseSemiunitaryMatrixFactories.cpp
            )
endif ()

if (ARMADILLO_FOUND AND Eigen3_FOUND)
    target_sources(tests_tools PRIVATE
            tools/AllSymmetricMatrixFactories.cpp
            )
    target_sources(spinner_test PRIVATE unit_tests/linear_algebra/armaEigenComparison_tests.cpp)
    target_sources(spinner_time_test PRIVATE
            time_tests/linear_algebra/armaEigenPerfomance_tests.cpp
            )
endif ()

target_link_libraries(tests_tools gtest gtest_main)

target_link_libraries(spinner_test gtest gtest_main tests_tools
        group space model spectrum common matrix magnetic_susceptibility nonlinear_solver spin_algebra eigendecompositor
        input)

target_link_libraries(spinner_time_test gtest gtest_main tests_tools
        group space model spectrum common matrix magnetic_susceptibility nonlinear_solver eigendecompositor)
add_library(eigendecompositor
        ExactEigendecompositor.cpp
        ImplicitSSquareEigendecompositor.cpp ImplicitSSquareEigendecompositor.h
        OneSymbolInHamiltonianEigendecompositor.cpp OneSymbolInHamiltonianEigendecompositor.h
        AbstractEigendecompositor.cpp
        ExplicitQuantitiesEigendecompositor.cpp ExplicitQuantitiesEigendecompositor.h
        EigendecompositorConstructor.cpp EigendecompositorConstructor.h)

target_link_libraries(eigendecompositor model space matrix spectrum)

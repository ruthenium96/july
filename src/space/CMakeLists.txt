add_library(space
        ../entities/BlockProperties.h ../entities/BlockProperties.cpp
        Space.h Space.cpp
        Subspace.h Subspace.cpp
        optimization/Symmetrizer.cpp optimization/Symmetrizer.h
        optimization/TzSorter.cpp optimization/TzSorter.h
        optimization/PositiveProjectionsEliminator.cpp optimization/PositiveProjectionsEliminator.h
        optimization/NonAbelianSimplifier.cpp optimization/NonAbelianSimplifier.h
        optimization/OptimizedSpaceConstructor.cpp optimization/OptimizedSpaceConstructor.h
        optimization/S2Transformer.cpp optimization/S2Transformer.h)

target_link_libraries(space spin_algebra data_structures common)
add_library(matrix
        matrix/Matrix.cpp matrix/Matrix.h
        matrix/Submatrix.cpp matrix/Submatrix.h)

add_library(spectrum
        spectrum/Spectrum.cpp spectrum/Spectrum.h spectrum/Subspectrum.cpp spectrum/Subspectrum.h)

add_subdirectory(data_structures)

add_library(magnetic_susceptibility
        magnetic_susceptibility/MagneticSusceptibilityController.h
        magnetic_susceptibility/assistant/EnsembleAverager.cpp
        magnetic_susceptibility/assistant/EnsembleAverager.h
        magnetic_susceptibility/assistant/ExperimentalValuesWorker.cpp
        magnetic_susceptibility/assistant/ExperimentalValuesWorker.h
        magnetic_susceptibility/worker/UniqueGOnlySSquaredWorker.cpp
        magnetic_susceptibility/worker/UniqueGOnlySSquaredWorker.h
        magnetic_susceptibility/worker/CurieWeissWorker.cpp
        magnetic_susceptibility/worker/CurieWeissWorker.h
        magnetic_susceptibility/worker/BasicWorker.cpp
        magnetic_susceptibility/worker/BasicWorker.h
        magnetic_susceptibility/MagneticSusceptibilityController.cpp
        magnetic_susceptibility/worker/AbstractWorker.h
        magnetic_susceptibility/worker/GSzSquaredWorker.cpp
        magnetic_susceptibility/worker/GSzSquaredWorker.h
        magnetic_susceptibility/assistant/WeightingScheme.cpp
        magnetic_susceptibility/assistant/WeightingScheme.h
        magnetic_susceptibility/worker/WorkerConstructor.cpp magnetic_susceptibility/worker/WorkerConstructor.h)

target_link_libraries(magnetic_susceptibility common)
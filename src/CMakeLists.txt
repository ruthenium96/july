# Fast hashmap:
FetchContent_Declare(
        emhash
        GIT_REPOSITORY https://github.com/ktprime/emhash.git
        GIT_TAG "2e3edc378373b331afa8bd570d657bf3f1229b86"
)
FetchContent_MakeAvailable(emhash)
include_directories(${emhash_SOURCE_DIR}/)

add_subdirectory(common)
add_subdirectory(eigendecompositor)
add_subdirectory(entities)
add_subdirectory(group)
add_subdirectory(input)
add_subdirectory(model)
add_subdirectory(nonlinear_solver)
add_subdirectory(space)
add_subdirectory(spin_algebra)

add_executable(spinner_main main.cpp)
target_link_libraries(spinner_main input)
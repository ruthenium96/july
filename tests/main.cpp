#include <gtest/gtest.h>

#include "src/common/Logger.h"

int main(int argc, char** argv) {
    common::Logger::set_level(common::PrintLevel::off);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
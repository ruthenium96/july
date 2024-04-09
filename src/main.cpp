#include <iostream>

#include "src/common/Logger.h"
#include "src/common/runner/Executer.h"
#include "src/input/Parser.h"

int main(int argc, char *argv[]) {
    common::Logger::set_pattern("[%C/%m/%d %H:%M:%S.%e] %v");
    if (argc == 1) {
        common::Logger::error_msg("Missing input file.");
        return -1;
    }
    if (argc > 2) {
        common::Logger::error_msg("Too much arguments.");
        return -1;
    }
    common::Logger::basic_msg("The calculations have started");
    common::Logger::separate(0, common::basic);
    try {
        auto parser = input::Parser(argv[1]);
        runner::Executer::execute(parser);
    } catch (const std::exception& e) {
        common::Logger::error_msg("An error occurred during program execution:\n{}", e.what());
        return -1;
    }
    common::Logger::basic_msg("The calculations have been completed successfully.");
    return 0;
}
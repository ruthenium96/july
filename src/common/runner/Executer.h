#ifndef SPINNER_EXECUTER_H
#define SPINNER_EXECUTER_H

#include "Runner.h"
#include "src/input/Parser.h"

namespace runner {

class Executer {
  public:
    static void execute(input::Parser parser);
};

}  // namespace runner

#endif  //SPINNER_EXECUTER_H

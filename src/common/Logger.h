#ifndef SPINNER_LOGGER_H
#define SPINNER_LOGGER_H

#include "PrintingFunctions.h"
#include <spdlog/sinks/stdout_sinks.h>
#include <spdlog/spdlog.h>

namespace common {

enum PrintLevel {
    trace, debug, verbose, detailed, basic, error, off
};

class Logger {
  public:
    template<typename... Args>
    static void set_pattern(Args&&... args) {
        logger_with_message->set_pattern(args...);
        logger_without_message->set_pattern("%v");
    }
    static void set_level(PrintLevel level) {
        auto levelSpdlog = static_cast<spdlog::level::level_enum>(level);
        logger_with_message->set_level(levelSpdlog);
        logger_without_message->set_level(levelSpdlog);
    }
    static PrintLevel get_level() {
        return static_cast<PrintLevel>(logger_with_message->level());
    }

    static void separate(size_t level, PrintLevel printLevel) {
        size_t number = 64 >> level;
        std::string s(number, '*');
        log(printLevel, "{}", s);
    };

    template <typename... Args>
    static void log(PrintLevel level, fmt::format_string<Args...> s, Args&&... args) {
        logger_without_message->log(
            static_cast<spdlog::level::level_enum>(level),
            s, std::forward<Args>(args)...);
    }
    template <typename... Args>
    static void trace(fmt::format_string<Args...> s, Args&&... args) {
        log(common::PrintLevel::trace, s, std::forward<Args>(args)...);
    }
    template <typename... Args>
    static void debug(fmt::format_string<Args...> s, Args&&... args) {
        log(common::PrintLevel::debug, s, std::forward<Args>(args)...);
    }
    template <typename... Args>
    static void verbose(fmt::format_string<Args...> s, Args&&... args) {
        log(common::PrintLevel::verbose, s, std::forward<Args>(args)...);
    }
    template <typename... Args>
    static void detailed(fmt::format_string<Args...> s, Args&&... args) {
        log(common::PrintLevel::detailed, s, std::forward<Args>(args)...);
    }
    template <typename... Args>
    static void basic(fmt::format_string<Args...> s, Args&&... args) {
        log(common::PrintLevel::basic, s, std::forward<Args>(args)...);
    }
    template <typename... Args>
    static void error(fmt::format_string<Args...> s, Args&&... args) {
        log(common::PrintLevel::error, s, std::forward<Args>(args)...);
    }

    template <typename... Args>
    static void log_msg(PrintLevel level, fmt::format_string<Args...> s, Args&&... args) {
        logger_with_message->log(
            static_cast<spdlog::level::level_enum>(level),
            s, std::forward<Args>(args)...);
    }
    template <typename... Args>
    static void debug_msg(fmt::format_string<Args...> s, Args&&... args) {
        log_msg(common::PrintLevel::debug, s, std::forward<Args>(args)...);
    }
    template <typename... Args>
    static void verbose_msg(fmt::format_string<Args...> s, Args&&... args) {
        log_msg(common::PrintLevel::verbose, s, std::forward<Args>(args)...);
    }
    template <typename... Args>
    static void detailed_msg(fmt::format_string<Args...> s, Args&&... args) {
        log_msg(common::PrintLevel::detailed, s, std::forward<Args>(args)...);
    }
    template <typename... Args>
    static void basic_msg(fmt::format_string<Args...> s, Args&&... args) {
        log_msg(common::PrintLevel::basic, s, std::forward<Args>(args)...);
    }
    template <typename... Args>
    static void error_msg(fmt::format_string<Args...> s, Args&&... args) {
        log_msg(common::PrintLevel::error, s, std::forward<Args>(args)...);
    }

  private:
    inline static std::shared_ptr<spdlog::logger>
        logger_with_message = spdlog::stdout_logger_mt("with_message");
    inline static std::shared_ptr<spdlog::logger>
        logger_without_message = spdlog::stdout_logger_mt("without_message");
};

class SpdlogStream : public std::stringbuf {
  public:
    explicit SpdlogStream(PrintLevel level): level_(level) {};

    int sync() override {
        auto str = this->str();
        common::Logger::log(level_, "{}", str);
        this->str("");
        return 0;
    }
  private:
    PrintLevel level_;
};

class SpdlogStreamMsg : public std::stringbuf {
  public:
    explicit SpdlogStreamMsg(PrintLevel level): level_(level) {};

    int sync() override {
        auto str = this->str();
        common::Logger::log_msg(level_, "{}", str);
        this->str("");
        return 0;
    }
  private:
    PrintLevel level_;
};
}

#endif  //SPINNER_LOGGER_H
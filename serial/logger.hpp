// logger.hpp
#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <string>

class Logger {
public:
    Logger(std::string program);
    void minfo(const std::string &msg);
    void ginfo(const std::string &msg);
    void cinfo(const std::string &msg);
    void warning(const std::string &msg);
    void error(const std::string &msg);
    void output(const std::string &msg);

private:
    std::string program_;
};

#endif // LOGGER_HPP

// logger.cpp
#include "logger.hpp"
#include "const.hpp"
#include <iostream>

using namespace std;

Logger::Logger(string program) : program_(move(program)) {}

void Logger::minfo(const string &msg)
{
    string buf = CYAN_COLOR + "[INF:" + program_ + "] " + RESET_COLOR + msg + '\n';
    cerr << buf;
}

void Logger::ginfo(const string &msg)
{
    string buf = MAGENTA_COLOR + "[INF:" + program_ + "] " + RESET_COLOR + msg + '\n';
    cerr << buf;
}

void Logger::cinfo(const string &msg)
{
    string buf = YELLOW_COLOR + "[INF:" + program_ + "] " + RESET_COLOR + msg + '\n';
    cerr << buf;
}

void Logger::warning(const string &msg)
{
    string buf = YELLOW_COLOR + "[WARN:" + program_ + "] " + RESET_COLOR + msg + '\n';
    cerr << buf;
}

void Logger::error(const string &msg)
{
    string buf = RED_COLOR + "[ERR:" + program_ + "] " + RESET_COLOR + msg + '\n';
    cerr << buf;
}

void Logger::output(const string &msg)
{
    string buf = GREEN_COLOR + "[OUT:" + program_ + "] " + RESET_COLOR + msg + '\n';
    cerr << buf;
}

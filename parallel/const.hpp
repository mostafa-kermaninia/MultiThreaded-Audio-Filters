// const.hpp
#ifndef CONST_HPP
#define CONST_HPP

#include <string>

// Paths
const std::string PIPE_PATH = "./pipes/";
const std::string EXE_CITY = "./city.out";
const std::string EXE_GOODS = "./goods.out";
const std::string EXE_MAIN = "./main.out";

// Named pipe prefixes
const std::string GOODS_PIPE_PREFIX = "/tmp/";

// Constants
const int MAX_BUFFER_SIZE = 2048;

// Colors for logging
const std::string RESET_COLOR = "\033[0m";
const std::string RED_COLOR = "\033[1;31m";
const std::string GREEN_COLOR = "\033[1;32m";
const std::string YELLOW_COLOR = "\033[1;33m";
const std::string BLUE_COLOR = "\033[1;34m";
const std::string CYAN_COLOR = "\033[1;36m";
const std::string MAGENTA_COLOR = "\033[1;35m";

#endif // CONST_HPP

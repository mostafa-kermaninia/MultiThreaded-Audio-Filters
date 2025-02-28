# Compiler used
CXX = g++

# Compilation flags
CXXFLAGS = -O2 -Wall -I. -std=c++11

# Linking flags
LDFLAGS =
LIBS = -lsndfile -lm

# Source and object files
SOURCES = main.cpp logger.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# Final output
TARGET = VoiceFilters.out

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS) $(LIBS) -o $(TARGET)

# Separate pattern rules for each source file if you want them to
# recompile when headers change. We assume main.cpp includes some headers,
# and logger.cpp includes logger.hpp and const.hpp.
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

logger.o: logger.cpp logger.hpp const.hpp
	$(CXX) $(CXXFLAGS) -c logger.cpp -o logger.o

clean:
	rm -f $(OBJECTS) $(TARGET)

CXX := g++
CXXFLAGS := -fdiagnostics-color=always -g3 -Wall -fopenmp
INCLUDES := -I/usr/include/eigen3/ -I/usr/include/ -I/usr/include/eigen3/Eigen

# Update the library paths accordingly
LDFLAGS := -L/usr/lib
LDLIBS := -lgmsh -larmadillo -llapack -lopenblas -lgomp

SRC := $(wildcard *.cpp) $(wildcard src/*.cpp)
HEADER := $(wildcard src/*.h)
OBJ := $(SRC:.cpp=.o)

TARGET := TOTEM

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ) $(HEADER)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

%.o: %.cpp $(HEADER)
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDES)

clean:
	rm -f $(OBJ) $(TARGET)

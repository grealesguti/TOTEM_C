CXX := C:/msys64/mingw64/bin/g++.exe
CXXFLAGS := -fdiagnostics-color=always -g3 -Wall -fopenmp
INCLUDES =  -I/C:/msys64/usr/bin -I/C:/msys64/usr/include -I/C:/msys64/mingw64/include/ -I/C:/msys64/mingw64/include/armadillo_bits -I/C:/msys64/mingw64/include/openblas
LDFLAGS = -L/C:/msys64/mingw64/lib/ -I/C:/msys64/usr/lib -L/C:/msys64/mingw64/lib/gmsh-4.11.1.dist-info -L/C:/msys64/mingw64/lib/LAPACK
LDLIBS = -lgmsh -larmadillo -llapack -lopenblas -lgomp

SRC := $(wildcard *.cpp) $(wildcard src/*.cpp)
HEADER := $(wildcard src/*.h)
OBJ := $(SRC:.cpp=.o)

TARGET := TOTEM.exe

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ) $(HEADER)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

%.o: %.cpp $(HEADER)
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDES)

clean:
	rm -f $(OBJ) $(TARGET)

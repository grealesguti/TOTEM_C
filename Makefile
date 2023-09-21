SRC_DIR     = src
HEADER_DIR  = include
BUILD_DIR   = build
BIN_DIR     = bin

CXX := g++
CXXFLAGS := -fdiagnostics-color=always -g3 -Wall -fopenmp
INCLUDES := -I/usr/include/eigen3 -I$(HEADER_DIR)

# Update the library paths accordingly
LDFLAGS := -L/usr/lib
LDLIBS  := -lgmsh -larmadillo -llapack -lopenblas -lgomp

SRC     := $(wildcard $(SRC_DIR)/*.cpp)
HEADER  := $(wildcard $(HEADER_DIR)/*.hpp)
OBJ     := $(SRC:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

TARGET := TOTEM

.PHONY: all clean

all: $(BIN_DIR)/$(TARGET)

# Link
$(BIN_DIR)/$(TARGET): $(OBJ) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

# Compile
$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDES)

# Create folder
$(BIN_DIR) $(BUILD_DIR):
	mkdir -p $@

clean:
	rm -f $(OBJ) $(BIN_DIR)/$(TARGET)

CXX = g++
CXXFLAGS = -std=c++20 -O3 -I ./eigen-3.4.0
TARGET = Dcb2D.exe
SRC = Dcb2D.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@
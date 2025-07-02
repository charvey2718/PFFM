CXX = g++
CXXFLAGS = -O3 -fno-math-errno -DNDEBUG -march=native -I ./eigen-3.4.0
TARGETI = Dcb2D-ModeI.exe
TARGETII = Dcb2D-ModeII.exe
SRCI = Dcb2D-ModeI.cpp
SRCII = Dcb2D-ModeII.cpp

all: $(TARGETI) $(TARGETII)

modeI: $(TARGETI)

modeII: $(TARGETII)

$(TARGETI): $(SRCI)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(TARGETII): $(SRCII)
	$(CXX) $(CXXFLAGS) $^ -o $@
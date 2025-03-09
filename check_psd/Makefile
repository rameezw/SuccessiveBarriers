INC_DIRS=/opt/homebrew/include/eigen3/
LIB_DIRS=/opt/homebrew/lib/
LIBS=-lmpfr 
CXXFLAGS=-std=c++14 -g -I$(INC_DIRS) -I/opt/homebrew/include -I/Users/srirams/tools/mpreal -L$(LIB_DIRS) 
CXX=g++
eigen-library-test: eigen-library-test.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^  $(LIBS)

clean:
	rm -f eigen-library-test

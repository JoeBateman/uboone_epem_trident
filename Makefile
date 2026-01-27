CXX = g++
# HEADER = include/form_factors.h
CFLAGS = -I include
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

TEG_v2: TEG_v2.cpp
	$(CXX) -o TEG_v2 TEG_v2.cpp $(CFLAGS) $(ROOTCFLAGS) $(ROOTLIBS)

clean:
	rm -f TEG_v2
#########################################################################
# Makefile for convertSensitivityFilesToFITS
#########################################################################
CXX           = g++
CXXFLAGS      = $(shell root-config --cflags) -I./inc/ -I/include/ -g
LD            = g++
OutPutOpt     = -o
LDFLAGS       = $(shell root-config --ldflags)
GLIBS         = $(shell root-config --libs) -lcfitsio

VPATH = src:inc

CONVERTFITS=		./obj/convertSensitivityFilesToFITS.o \
			        ./obj/VDL3IRFs.o

./obj/convertSensitivityFilesToFITS.o: 	./src/convertSensitivityFilesToFITS.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

convertSensitivityFilesToFITS:	$(CONVERTFITS)
	$(LD) $(LDFLAGS) $^ $(GLIBS) $(OutPutOpt) ./$@
	@echo "$@ done"

./obj/%.o:	%.cpp %.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	-rm -f ./obj/*.o

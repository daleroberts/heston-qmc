CXX = g++
LDFLAGS  ?= -L/usr/local/lib -L. -lQuantLib -lgfortran -lRmath
CXXFLAGS ?= -I/usr/local/include  -I/Library/Frameworks/R.framework/Headers -O3 -ftree-vectorize -std=c++0x
CXXFLAGS += -fopenmp
#CXXFLAGS += -DDEBUG
 
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

all: table1 table2 table3 checkomp

main: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)

asianSVJ: asianSVJ.o SVJtrajectory.o exact.o utils.o bessel.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp %.hpp %.h
	$(CXX) $(CXXFLAGS) -o $@ -c $<
    
LowDiscrepancy.o: LowDiscrepancy.f
	gfortran -O3 -o $@ -c $<

main.o: main.cpp asian.h
	$(CXX) $(CXXFLAGS) -o $@ -c main.cpp

test.o: test.cpp asianSVJ.h bridges.h asian.h
	$(CXX) $(CXXFLAGS) -o $@ -c $<

table2.o: table2.cpp bridges.h asian.h
	$(CXX) $(CXXFLAGS) -o $@ -c $<

table3.o: table3.cpp bridges.h asianSVJ.h
	$(CXX) $(CXXFLAGS) -o $@ -c $<

table4.o: table4.cpp bridges.h asian.h
	$(CXX) $(CXXFLAGS) -o $@ -c $<

eurosbridge.o: eurosbridge.cpp bridges.h asian.h
	$(CXX) $(CXXFLAGS) -o $@ -c $<

test: test.o utils.o bessel.o scramblepoints.o exact.o Sobolpoints.o LowDiscrepancy.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

table1: table1.o euro.o exact.o utils.o bessel.o Sobolpoints.o scramblepoints.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

table2: table2.o exact.o utils.o bessel.o Sobolpoints.o scramblepoints.o LowDiscrepancy.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

table3: table3.o euro.o exact.o utils.o bessel.o Sobolpoints.o scramblepoints.o LowDiscrepancy.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

table4: table4.o exact.o utils.o bessel.o Sobolpoints.o scramblepoints.o LowDiscrepancy.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

eurosbridge: eurosbridge.o exact.o utils.o bessel.o Sobolpoints.o scramblepoints.o LowDiscrepancy.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

checkomp: checkomp.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

hestcall: hestcall.o ql_imps.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

mccall: mccall.o ql_imps.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

mctrials: mctrials.o ql_imps.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

testread: testread.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

intbridge: intbridge.o exact.o utils.o bessel.o LowDiscrepancy.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o table1 table2 table3 table4 checkomp hestcall intbridge mccall mctrials testread eurosbridge test

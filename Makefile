
############# EDITABLE PART #############

# compiler options
CXX = g++
CXXFLAGS += -Wall -std=c++14 -I${mkEigenInc} ##-I/usr/include/mysql -fopenmp

# name of executable
EXEC=main##OPT_JR_CPP

#source directory
SRCDIR = src/
#object file directory (it expects an existing directory)
ODIR= obj/

#libraries
##LDLIBS += -lmysqlclient





############# NOT EDITABLE PART #############


#src files
SRCS= $(wildcard $(SRCDIR)*.cpp)
#header files
HEADERS=$(wildcard $(SRCDIR)*.hh)

#object files
_OBJ = $(SRCS:.cpp=.o)
OBJS = $(patsubst $(SRCDIR)%, $(ODIR)%, $(_OBJ))

.PHONY= all clean distclean

.DEFAULT_GOAL  = all

all: $(EXEC)

$(EXEC): $(OBJS) libdpp.a
	$(CXX) -o $@ $^ $(LDLIBS) $(CXXFLAGS) -L. -ldpp

libdpp.a : $(ODIR)*.o
	$(AR) crv libdpp.a $(ODIR)*.o

$(ODIR)%.o: $(SRCDIR)%.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

clean:
	-rm -f $(OBJS)

distclean: clean
	-rm libdpp.a
	-rm -f $(EXEC)

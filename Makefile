#
# 'make depend' uses makedepend to automatically generate dependencies 
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

# define the C++ compiler to use
CXX = g++

# define any compile-time flags
CXXFLAGS= -O2 -Wall -fPIC `root-config --cflags`
ROOTFLAGS  := $(shell root-config --cflags)

SRC_DIR := src
INC_DIR := headers
BIN_DIR := bin
OBJ_DIR := obj

# define any directories containing header files other than /usr/include
#INCLUDES = -I./Analysis -I../Anl
INCLUDES = -I$(INC_DIR)


# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
#LFLAGS = -L/home/newhall/lib  -L../lib
LDFLAGS = `root-config --libs`
ROOTLIBS := $(shell root-config --libs)

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
#LIBS = -lmylib -lm

# define the C source files
#SRCS = emitter.c error.c init.c lexer.c main.c symbol.c parser.c

# MUST HAVE main() in the very first .cpp code
       
ASYM_SRCS := $(SRC_DIR)/asym_main.cpp $(SRC_DIR)/asym_func.cpp $(SRC_DIR)/config_manager.cpp
KIN_SRCS := $(SRC_DIR)/kin_average.cpp $(SRC_DIR)/config_manager.cpp

# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
#OBJS = $(SRCS:.c=.o)
ASYM_OBJS := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(ASYM_SRCS))
KIN_OBJS := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(KIN_SRCS))

# define the executable file 
EXEC := asym_calc kin_average

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all:    $(EXEC)
	@echo asym_calc has been compiled

asym_calc: $(ASYM_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(ROOTLIBS)

kin_average: $(KIN_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(ROOTLIBS)
# $@ = target, $^ = all dependencies

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<  -o $@ $(ROOTFLAGS)

clean:
	rm -f $(OBJ_DIR)/*.o $(EXEC)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
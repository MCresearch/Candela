##############################################################################
# Candela MAKEFILE
##############################################################################

D_SRC=./src
D_OBJ=./obj

SRC_CPP = $(wildcard $(D_SRC)/*.cpp)
OBJ_CPP = $(addprefix $(D_OBJ)/, $(patsubst %.cpp, %.o, $(notdir $(SRC_CPP))))

TARGET=candela

#CC=icc -g
#CC=mpiicc 
#CC=mpicxx
#CC=CC
#HONG=-D__MPI
#CC=mpiicc
CC=g++
OPTION=-O3 

${TARGET}:$(OBJ_CPP)
	$(CC) -o $@ $^

$(D_OBJ)/%.o: $(D_SRC)/%.cpp
	$(CC) -c $< -o $@

.PHONY: clean

clean:
	rm -f $(D_OBJ)/* $(TARGET) 

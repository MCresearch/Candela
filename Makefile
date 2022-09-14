##############################################################################
# Candela MAKEFILE
##############################################################################

#--------------------------------------
#----------- Please set ---------------
#-------------------------------------- 
#mpi version: mpiicc/mpiicpc/mpicxx/mpicc 
CC=mpiicpc 
#serial version: icc/icpc/g++
# CC=g++
#--------------------------------------



OPTION=-O3 
D_SRC=./src
D_OBJ=./obj
D_BIN=./bin

SRC_CPP = $(wildcard $(D_SRC)/*.cpp)
OBJ_CPP = $(addprefix $(D_OBJ)/, $(patsubst %.cpp, %.o, $(notdir $(SRC_CPP))))
TARGET=$(D_BIN)/candela

ifeq ($(findstring mpi, $(CC)), mpi)
	HONG=-D__MPI
endif


${TARGET}:$(OBJ_CPP)
	@if [ ! -d $(D_BIN) ]; then mkdir $(D_BIN); fi
	$(CC) $(OPTION) -o $@ $^

$(D_OBJ)/%.o: $(D_SRC)/%.cpp
	@-if [ ! -d $(D_OBJ) ]; then mkdir $(D_OBJ); fi
	$(CC) $(OPTION) $(HONG) -c $< -o $@

.PHONY: clean

clean:
	@rm -rf $(D_OBJ) $(D_BIN)

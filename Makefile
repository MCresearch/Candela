##############################################################################
# Candela MAKEFILE
##############################################################################

#------------------------------------------------------------------
#--------------------------- Please set ---------------------------
#------------------------------------------------------------------
# mpi version: mpiicc/mpiicpc/mpicxx/mpicc 
CC=mpiicpc 
# serial version: icc/icpc/g++
# CC=g++

# Compile integrate-test version, CC must use g++/mpicxx/mpicc
# It can check the leak of memory.
TEST=OFF
#------------------------------------------------------------------



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

ifeq ($(TEST), ON)
#CC must be a gnu compiler, or else it will fail.
	OPTION += -fsanitize=address -fno-omit-frame-pointer
endif


${TARGET}:$(OBJ_CPP)
	@if [ ! -d $(D_BIN) ]; then mkdir $(D_BIN); fi
	$(CC) $(OPTION) -o $@ $^

$(D_OBJ)/%.o: $(D_SRC)/%.cpp $(D_OBJ)/README
	$(CC) $(OPTION) $(HONG) -c $< -o $@
	
$(D_OBJ)/README:
	@if [ ! -d $(D_OBJ) ]; then mkdir $(D_OBJ); fi
	@echo "This is a temporary dirctory to store obj files." > $(D_OBJ)/README


.PHONY: clean

clean:
	@rm -rf $(D_OBJ) $(D_BIN)

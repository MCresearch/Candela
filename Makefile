##############################################################################
# Candela MAKEFILE
##############################################################################

#------------------------------------------------------------------
#--------------------------- Please set ---------------------------
#------------------------------------------------------------------
# mpi version: mpiicc/mpiicpc/mpicxx
CXX=mpiicpc 
# serial version: icc/icpc/g++
# CXX=g++

# Compile integrate-test version, CXX must use g++/mpicxx
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

ifeq ($(findstring mpi, $(CXX)), mpi)
    HONG = -D__MPI
    MPICOMPILE = ON
endif

ifeq ($(TEST), ON)
#CXX must be a gnu compiler, or else it will fail.
    OPTION = -g -fsanitize=address -fno-omit-frame-pointer
	ifeq ($(findstring mpi, $(CXX)), mpi)
        CXX = mpicxx
	else
        CXX = g++
	endif
endif


${TARGET}:$(OBJ_CPP)
	@if [ ! -d $(D_BIN) ]; then mkdir $(D_BIN); fi
	$(CXX) $(OPTION) -o $@ $^
	@if [ $(TEST) == ON ]; then cd test;./Autotest.sh $(MPICOMPILE);cd ..; fi

$(D_OBJ)/%.o: $(D_SRC)/%.cpp $(D_OBJ)/readme.log
	$(CXX) $(OPTION) $(HONG) -c $< -o $@
	
$(D_OBJ)/readme.log:
	@if [ ! -d $(D_OBJ) ]; then mkdir $(D_OBJ); fi
	@echo "This is a temporary dirctory to store obj files." > $(D_OBJ)/readme.log


.PHONY: clean, test
test:
	@cd test;./Autotest.sh $(MPICOMPILE);cd ..;

clean:
	@rm -rf $(D_OBJ) $(D_BIN)

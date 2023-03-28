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

#openmp
OPENMP=ON

# Compile integrate-test version, CXX must use g++/mpicxx
# It can check the leak of memory.
TEST=OFF
#------------------------------------------------------------------



CFLAGS=-O3 
FFLAGS=-O3 -cpp
D_SRC=./src
D_OBJ=./obj
D_BIN=./bin

SRC_CPP = $(wildcard $(D_SRC)/*.cpp)
SRC_F90 = $(wildcard $(D_SRC)/*.f90)
OBJ_CPP = $(addprefix $(D_OBJ)/, $(patsubst %.cpp, %.o, $(notdir $(SRC_CPP))))
OBJ_F90 = $(addprefix $(D_OBJ)/, $(patsubst %.f90, %.o, $(notdir $(SRC_F90))))
TARGET=$(D_BIN)/candela

ifeq ($(findstring mpii, $(CXX)), mpii)
    FC = ifort
else
	ifeq ($(findstring mpi, $(CXX)), mpi)
        FC = gfortran
    else
	    ifeq ($(findstring i, $(CXX)), i)
            FC = ifort
		else
            FC = gfortran
		endif
	endif
endif


ifeq ($(findstring mpi, $(CXX)), mpi)
    HONG = -D__MPI
    MPICOMPILE = ON
endif

ifeq ($(TEST), ON)
    DEBUG = ON
endif

ifeq ($(DEBUG), ON)
#CXX must be a gnu compiler, or else it will fail.
    CFLAGS = -g -fsanitize=address -fno-omit-frame-pointer
    FFLAGS = -g -fsanitize=address -fno-omit-frame-pointer -cpp
	ifeq ($(findstring mpi, $(CXX)), mpi)
        CXX = mpicxx
	else
        CXX = g++
	endif
    FC = gfortran
endif

ifeq ($(OPENMP), ON)
    CFLAGS += -fopenmp
    FFLAGS += -fopenmp
endif

${TARGET}:$(OBJ_CPP) $(OBJ_F90)
	@if [ ! -d $(D_BIN) ]; then mkdir $(D_BIN); fi
	$(CXX) $(CFLAGS) -o $@ $^
	@if [ $(TEST) == ON ]; then cd test;./Autotest.sh $(MPICOMPILE);cd ..; fi

$(D_OBJ)/%.o: $(D_SRC)/%.cpp $(D_OBJ)/readme.log
	$(CXX) $(CFLAGS) $(HONG) -c $< -o $@

$(D_OBJ)/%.o: $(D_SRC)/%.f90 $(D_OBJ)/readme.log
	$(FC) $(FFLAGS) -c $< -o $@
	
$(D_OBJ)/readme.log:
	@if [ ! -d $(D_OBJ) ]; then mkdir $(D_OBJ); fi
	@echo "This is a temporary dirctory to store obj files." > $(D_OBJ)/readme.log


.PHONY: clean, test
test:
	@cd test;./Autotest.sh $(MPICOMPILE);cd ..;

clean:
	@rm -rf $(D_OBJ) $(D_BIN)

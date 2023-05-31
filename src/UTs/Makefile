include ../../Makefile.vars
#------------------------------------------------------------------
D_SRC=../
D_OBJ=../../obj

CXX=mpicxx
FC = gfortran
NVCC=nvcc

CFLAGS = -g -fsanitize=address -fno-omit-frame-pointer
FFLAGS = -g -fsanitize=address -fno-omit-frame-pointer -cpp
NVFLAGS= -g -Xcompiler="-fsanitize=address -fno-omit-frame-pointer" -lasan

LINK_LIBS = -L${GTEST_DIR}/lib -lgtest -lpthread


SRC_CPP = $(wildcard $(D_SRC)*.cpp *.cpp)
SRC_F90 = $(wildcard $(D_SRC)*.f90)
SRC_CU = $(wildcard $(D_SRC)*.cu)
OBJ_CPP = $(addprefix $(D_OBJ)/, $(filter-out main.o, $(patsubst %.cpp, %.o, $(notdir $(SRC_CPP)))))
OBJ_F90 = $(addprefix $(D_OBJ)/, $(patsubst %.f90, %.o, $(notdir $(SRC_F90))))
OBJ_CU = $(addprefix $(D_OBJ)/, $(patsubst %.cu, %.o, $(notdir $(SRC_CU))))
OBJS = $(OBJ_CPP) $(OBJ_F90)

TARGET=unittest.exe

ifeq ($(findstring mpi, $(CXX)), mpi)
    HONG = -D__MPI
endif

ifeq ($(OPENMP), ON)
    CFLAGS += -fopenmp
    FFLAGS += -fopenmp
endif

#CUDA
CUDA_INC = $(CUDA_ROOT)/include
CUDA_LIB = -L$(CUDA_ROOT)/lib64 -lcudart -lcublas
ifeq ($(CUDA), ON)
    OBJS += $(OBJ_CU)
    LINK_LIBS += $(CUDA_LIB)
    HONG += -DUSE_CUDA
endif

${TARGET}:$(OBJS)
	$(CXX) $(CFLAGS) -o $@ $^ $(LINK_LIBS)

$(D_OBJ)/%.o: %.cpp $(D_OBJ)/readme.log
	$(CXX) $(CFLAGS) $(HONG) -c $< -o $@

$(D_OBJ)/%.o: $(D_SRC)/%.cpp $(D_OBJ)/readme.log
	$(CXX) $(CFLAGS) $(HONG) -c $< -o $@

$(D_OBJ)/%.o: $(D_SRC)/%.f90 $(D_OBJ)/readme.log
	$(FC) $(FFLAGS) -c $< -o $@

$(D_OBJ)/%.o: $(D_SRC)/%.cu $(D_OBJ)/readme.log
	$(NVCC) $(NVFLAGS) -c $< -o $@
	
$(D_OBJ)/readme.log:
	@if [ ! -d $(D_OBJ) ]; then mkdir $(D_OBJ); fi
	@echo "This is a temporary dirctory to store obj files." > $(D_OBJ)/readme.log

test:
	@${MAKE} ${TARGET}
	@OMP_NUM_THREADS=2 ASAN_OPTIONS=protect_shadow_gap=0 ./${TARGET}
clean:
	@rm -rf ${TARGET} $(D_OBJ)

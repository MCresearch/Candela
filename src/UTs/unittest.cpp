#include "mpi.h"
#include "unittest.h"
namespace UT
{
    int NPROC = 0, RANK = 0;
}
int main(int argc, char **argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &UT::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &UT::RANK);
#endif
    int result = 0;
    testing::AddGlobalTestEnvironment(new TestEnv);
    testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
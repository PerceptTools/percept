#include <gtest/gtest.h>
#include <cstdlib>
#include <mpi.h>

int main(int argc, char **argv) { \
  if ( MPI_SUCCESS != MPI_Init( & argc , & argv ) ) { \
    std::cerr << "MPI_Init FAILED" << std::endl ; \
    std::abort(); \
  } \
  std::cout << "Running main() from gtest_main.cc\n"; \
  testing::InitGoogleTest(&argc, argv); \
  bool result = RUN_ALL_TESTS(); \
  MPI_Finalize(); \
  return result; \
}

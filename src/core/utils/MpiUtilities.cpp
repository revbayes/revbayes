#include "MpiUtilities.h"

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"

#include <iostream>
#include <sstream>

#ifdef RB_MPI
#include <mpi.h>
#endif

void RevBayesCore::MpiUtilities::DebugWait(int rank)
{
    
#ifdef RB_MPI
    char	a;
    if (rank == 0) {
    	scanf("%c", &a);
    	printf("%d: Starting now\n", rank);
    }
    
    MPI_Bcast(&a, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
    printf("%d: Starting now\n", rank);
#endif
}

void RevBayesCore::MpiUtilities::DebugMsg(const std::stringstream& s)
{
#ifdef RB_MPI
#ifdef DEBUG_MPI_MCA
    int pid = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    std::cout << pid << "   before: " << s.str() << "\n";
    MPI_Barrier();
    std::cout << pid << "   after:  " << s.str() << "\n";
#endif
#endif
}

void RevBayesCore::MpiUtilities::DebugMsg(const std::string& s) {
#ifdef RB_MPI
#ifdef DEBUG_MPI_MCA
    int pid = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    std::cout << pid << "   before: " << s << "\n";
    MPI_Barrier();
    std::cout << pid << "   after:  " << s << "\n";
#endif
#endif
}

void RevBayesCore::MpiUtilities::DebugMsg(const std::string& s, int x) {
#ifdef RB_MPI
#ifdef DEBUG_MPI_MCA
    std::stringstream ss;
    ss << s << " " << x;
    int pid = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    std::cout << pid << "   before: " << ss.str() << "\n";
    MPI_Barrier();
    std::cout << pid << "   after:  " << ss.str() << "\n";
#endif
#endif
}

void RevBayesCore::MpiUtilities::DebugMsg(const std::string& s, double x) {
#ifdef RB_MPI
#ifdef DEBUG_MPI_MCA
    std::stringstream ss;
    ss << s << " " << x;
    int pid = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    std::cout << pid << "   before: " << ss.str() << "\n";
    MPI_Barrier();
    if (pid == 0) std::cout << "\n";
    MPI_Barrier();
    std::cout << pid << "   after:  " << ss.str() << "\n";
#endif
#endif
}

void RevBayesCore::MpiUtilities::DebugMsgPid(const std::string& s, int p)
{
#ifdef RB_MPI
#ifdef DEBUG_MPI_MCA
    int pid = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Barrier();
    if (pid == p) {
        std::cout << s;
    }
    MPI_Barrier();
#endif
#endif
}


// NOTE: This does more than just synchronize all the copies of the global RNG.
//       It also resets them to the common starting seed.
#ifdef RB_MPI
void RevBayesCore::MpiUtilities::synchronizeRNG( const MPI_Comm &analysis_comm )
#else
void RevBayesCore::MpiUtilities::synchronizeRNG( void )
#endif
{
    unsigned int seed = 0;

    int process_id = 0;
    #ifdef RB_MPI
    MPI_Comm_rank(analysis_comm, &process_id);
    #endif

    // sync the random number generators
    if ( process_id == 0 )
    {
        seed = RevBayesCore::GLOBAL_RNG->getNewSeed();
        seed = RevBayesCore::GLOBAL_RNG->getSeed();
    }

    #ifdef RB_MPI
    MPI_Bcast(&seed, 1, MPI_INT, 0, analysis_comm);
    #endif

    
    RevBayesCore::GLOBAL_RNG->setSeed( seed );

}

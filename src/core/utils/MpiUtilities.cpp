#include "MpiUtilities.h"

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"

#include <iostream>
#include <sstream>

#ifdef RB_MPI
#include <mpi.h>
#endif

void RevBayesCore::MpiUtilities::DebugWait(int rank)
{
    
#ifdef RB_MPI
    char	a;
    if (rank == 0)
    {
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
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << pid << "   after:  " << s.str() << "\n";
#endif
#endif
}

#include <chrono>
#include <thread>
void RevBayesCore::MpiUtilities::DebugMsg(const std::string& s)
{
#ifdef RB_MPI
#ifdef DEBUG_MPI_MCA
    int pid = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    std::cout << pid << "   before: " << s << "\n";
    std::cout.flush();
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    using namespace std::this_thread;     // sleep_for, sleep_until
    using namespace std::chrono_literals; // ns, us, ms, s, h, etc.
    using std::chrono::system_clock;

    sleep_for(1s);
    
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << pid << "   after:  " << s << "\n";

    sleep_for(1s);
#endif
#endif
}

void RevBayesCore::MpiUtilities::DebugMsg(const std::string& s, int x)
{
#ifdef RB_MPI
#ifdef DEBUG_MPI_MCA
    std::stringstream ss;
    ss << s << " " << x;
    int pid = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    std::cout << pid << "   before: " << ss.str() << "\n";
    std::cout.flush();
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    using namespace std::this_thread;     // sleep_for, sleep_until
    using namespace std::chrono_literals; // ns, us, ms, s, h, etc.
    using std::chrono::system_clock;

    sleep_for(1s);
    
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << pid << "   after:  " << ss.str() << "\n";
    
    sleep_for(1s);
#endif
#endif
}

void RevBayesCore::MpiUtilities::DebugMsg(const std::string& s, double x)
{
#ifdef RB_MPI
#ifdef DEBUG_MPI_MCA
    std::stringstream ss;
    ss << s << " " << x;
    int pid = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    std::cout << pid << "   before: " << ss.str() << "\n";
    MPI_Barrier(MPI_COMM_WORLD);
    if (pid == 0) std::cout << "\n";
    MPI_Barrier(MPI_COMM_WORLD);
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
    MPI_Barrier(MPI_COMM_WORLD);
    if (pid == p)
    {
        std::cout << s;
    }
    MPI_Barrier(MPI_COMM_WORLD);
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
        seed = int(RevBayesCore::GLOBAL_RNG->uniform01() * RbConstants::Integer::max);
    }

    #ifdef RB_MPI
    MPI_Bcast(&seed, 1, MPI_INT, 0, analysis_comm);
    #endif

    
    RevBayesCore::GLOBAL_RNG->setSeed( seed );

}

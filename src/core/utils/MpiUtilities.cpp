#include "MpiUtilities.h"
#include <iostream>
#include <sstream>
#include <string>
#include <cstdio>
#include <stdio.h>

#ifdef RB_MPI
#include <mpi.h>
#endif

void RevBayesCore::MpiUtilities::DebugWait(int rank) {
    
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

void RevBayesCore::MpiUtilities::DebugMsg(const std::stringstream& s) {
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

void RevBayesCore::MpiUtilities::DebugMsgPid(const std::string& s, int p) {
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

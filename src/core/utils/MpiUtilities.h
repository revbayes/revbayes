#ifndef MpiUtilities_H

#ifdef RB_MPI
#include <mpi.h>
#endif

#include <sstream>

namespace RevBayesCore {
    
    namespace MpiUtilities {
        
        void DebugWait(int rank);
        void DebugMsg(const std::stringstream& s);
        void DebugMsg(const std::string& s);
        void DebugMsg(const std::string& s, int x);
        void DebugMsg(const std::string& s, double x);
        void DebugMsgPid(const std::string& s, int pid);

#ifdef RB_MPI
        void synchronizeRNG(const MPI_Comm &analysis_comm);
#else
        void synchronizeRNG(void);
#endif


    }
}

#endif /* MpiUtilities_H */


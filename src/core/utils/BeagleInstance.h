#ifndef BEAGLE_INSTANCE_H
#define BEAGLE_INSTANCE_H

#if defined (RB_BEAGLE)

#include "RbSettings.h"
#include "RlUserInterface.h"
#include "BeagleUtilities.h"
#include "libhmsbeagle/beagle.h"

#include <vector>


namespace RevBayesCore {

    class BeagleInstance {

    public:

        BeagleInstance ( void );        

        ~BeagleInstance ( void );
        
        // overloaded operators
        BeagleInstance&     operator=(const BeagleInstance &i);

        
        void                createBEAGLE( int  b_tipCount
                                        , int  b_partialsBufferCount
                                        , int  b_compactBufferCount
                                        , int  b_stateCount
                                        , int  b_patternCount
                                        , int  b_eigenBufferCount
                                        , int  b_matrixBufferCount
                                        , int  b_categoryCount
                                        , int  b_scaleBufferCount
                                        );
        void                freeBEAGLE ( void );
        
        // getters
        int                 getResourceID ( void ) { return handle; };
        double              getStoredLnLikelihood ( void ) { return stored_lnLikelihood; };

        // setters
        void                setStoredLnLikelihood ( double v ) { this->stored_lnLikelihood = v; };
        void                setCPUThreadCount(size_t n);

        // functions
        long                generateSettingsBitmap( void );
        void                pushBeagleOp( BeagleOperation );
        void                pushNode( int );
        void                pushBranch( double );
        void                clearQueues( void );

    private:

        int                 handle;
        int                 resourcenumber;
        std::string         resourcename;
        unsigned            nstates;
        unsigned            nratecateg;
        unsigned            npatterns;
        unsigned            nthreads;
        unsigned            partial_offset;
        unsigned            tmatrix_offset;

        // Collection 'queues' to be filled during tree traversal
        std::vector<BeagleOperation>   b_ops;
        std::vector<int>               b_node_indices;
        std::vector<double>            b_branch_lengths;

        // Useful when we need to collect likelihoods accross multiple instances
        double                         stored_lnLikelihood;

        // Sebastian: not used (yet)
        //std::vector<unsigned>       subsets;
    };

}

#endif

#endif //-- end BEAGLE_INSTANCE_H

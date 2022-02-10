#ifndef BEAGLE_INSTANCE_H
#define BEAGLE_INSTANCE_H

#include "RbSettings.h"
#include "RlUserInterface.h"
#include "BeagleUtilities.h"
#include "libhmsbeagle/beagle.h"


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
        void                freeBEAGLE(void);
        
        // getters
        int                 getResourceID( void ) { return handle; };
        
        // setters
        void                setCPUThreadCount(size_t n);
        
    private:

        int                         handle;
        int                         resourcenumber;
        std::string                 resourcename;
        unsigned                    nstates;
        unsigned                    nratecateg;
        unsigned                    npatterns;
        unsigned                    nthreads;
        unsigned                    partial_offset;
        unsigned                    tmatrix_offset;
        
        // Sebastian: not used (yet)
//        std::vector<unsigned>       subsets;
    };

}

#endif //-- end BEAGLE_INSTANCE_H

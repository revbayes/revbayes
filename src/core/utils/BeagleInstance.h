#ifndef BEAGLE_INSTANCE_H
#define BEAGLE_INSTANCE_H

#include "RbSettings.h"
#include "RlUserInterface.h"
#include "BeagleUtilities.h"
#include "libhmsbeagle/beagle.h"


namespace RevBayesCore {

    class BeagleInstance {

    private:

        int resourceID = -1;

    protected:

        static BeagleInstance *beagle_singleton;

        BeagleInstance ( int  b_resource
                       , bool b_use_cpu_threading
                       , bool b_use_scaling
                       , int  b_tipCount
                       , int  b_partialsBufferCount
                       , int  b_compactBufferCount
                       , int  b_stateCount
                       , int  b_patternCount
                       , int  b_eigenBufferCount
                       , int  b_matrixBufferCount
                       , int  b_categoryCount
                       , int  b_scaleBufferCount
                       );

    public:

        ~BeagleInstance ( void );

        static BeagleInstance *getInstance  ( int  b_resource
                                            , bool b_use_cpu_threading
                                            , bool b_use_scaling
                                            , int  b_tipCount
                                            , int  b_partialsBufferCount
                                            , int  b_compactBufferCount
                                            , int  b_stateCount
                                            , int  b_patternCount
                                            , int  b_eigenBufferCount
                                            , int  b_matrixBufferCount
                                            , int  b_categoryCount
                                            , int  b_scaleBufferCount
                                            );

        static int            getResourceID ( void ) { return (beagle_singleton ? beagle_singleton->resourceID : -1); };

        //static bool           hasInstance   ( void ) { return (beagle_singleton ? true : false); };

    };

}

#endif //-- end BEAGLE_INSTANCE_H

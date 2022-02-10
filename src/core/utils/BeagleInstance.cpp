#include <iostream>

#include "RbException.h"
#include "RbSettings.h"
#include "RlUserInterface.h"
#include "BeagleUtilities.h"
#include "libhmsbeagle/beagle.h"

#include "BeagleInstance.h"

using namespace RevBayesCore;

#define RB_BEAGLE_INFO


BeagleInstance::BeagleInstance ( void ) :
    handle(-1),
    resourcenumber(-1),
    resourcename(""),
    nstates(0),
    nratecateg(0),
    npatterns(0),
    nthreads(1),
    partial_offset(0),
    tmatrix_offset(0)
{
    // we don't do anything yet
}




BeagleInstance::~BeagleInstance ( )
{
    // free our BEAGLE instance
    freeBEAGLE();
}


BeagleInstance& BeagleInstance::operator=(const BeagleInstance &i)
{
    if ( this != &i )
    {
        // reset internal flags
        handle              = -1;
        resourcenumber      = -1;
        resourcename        = "";
        nstates             = 0;
        nratecateg          = 0;
        npatterns           = 0;
        partial_offset      = 0;
        tmatrix_offset      = 0;
    }
    
    return *this;
}



void BeagleInstance::createBEAGLE(  int  b_tipCount
                                  , int  b_partialsBufferCount
                                  , int  b_compactBufferCount
                                  , int  b_stateCount
                                  , int  b_patternCount
                                  , int  b_eigenBufferCount
                                  , int  b_matrixBufferCount
                                  , int  b_categoryCount
                                  , int  b_scaleBufferCount
                                  )
{
    std::stringstream ss;
    RBOUT ( "Initializing BEAGLE" );
    
#if defined ( RB_BEAGLE_INFO )
        ss << std::endl;
        ss << "Using BEAGLE library v" << beagleGetVersion();
        ss << " for parallel likelihood evaluation (https://beagle-dev.github.io/)";
        ss << std::endl;
#if defined ( RB_USE_EIGEN3 )
        ss << "Using Eigen3 library for eigensystem calculations." << std::endl;
#else
        ss << "Using Revbayes internal methods for eigensystem calculations." << std::endl;
#endif /* RB_USE_EIGEN3 */
#endif /* RB_BEAGLE_INFO */
    
    int  b_resource            = 0;
    bool b_use_cpu_threading   = RbSettings::userSettings().getBeagleMaxCPUThreads() != 1
                               ? true : false;
    bool b_use_scaling         = RbSettings::userSettings().getBeagleScalingMode()   != "none"
                               ? true : false;
    b_use_scaling = true;
    
    BeagleInstanceDetails b_return_info;

    int* b_resourceList        = &b_resource;
    int  b_resourceCount       = 1;
    long b_requirementFlags    = 0;
//    long b_preferenceFlags     = ( RbSettings::userSettings().getBeagleUseDoublePrecision()
//                                 ? BEAGLE_FLAG_PRECISION_DOUBLE
//                                 : BEAGLE_FLAG_PRECISION_SINGLE
//                                 )
//                               | ( b_use_cpu_threading
//                                 ? BEAGLE_FLAG_THREADING_CPP
//                                 : 0
//                                 )
//                               | BEAGLE_FLAG_SCALING_MANUAL
//                               | BEAGLE_FLAG_SCALERS_LOG;

    long b_preferenceFlags     = BEAGLE_FLAG_PRECISION_SINGLE | BEAGLE_FLAG_THREADING_CPP;
    
    if ( RbSettings::userSettings().getBeagleAuto() == true )
    {
        ss << "Running benchmarks to automatically select fastest BEAGLE resource... ";
        ss << std::endl;

        // select fastest resource
        BeagleBenchmarkedResourceList* rBList;
        rBList = beagleGetBenchmarkedResourceList( b_tipCount
                                                 , b_compactBufferCount
                                                 , b_stateCount
                                                 , b_patternCount
                                                 , b_categoryCount
                                                 , NULL  // resourceList
                                                 , 0     // resourceCount
                                                 , b_preferenceFlags
                                                 , b_requirementFlags
                                                 , b_eigenBufferCount
                                                 , 1     // partitionCount
                                                 , 0     // calculateDerivatives
                                                 , 0     // benchmarkFlags
                                                 );

        if (rBList != NULL)
        {
            b_resource = rBList->list[0].number;
            ss << "\tUsing resource " << rBList->list[0].number << ": " << rBList->list[0].name;
            if ( rBList->list[0].number != 0 )
            {
                ss << " (" << rBList->list[0].performanceRatio << "x CPU)";
            }
            ss << std::endl;
        } else {
            ss << "\tResource benchmarking failed, using resource "
               << b_resource << ": " << rBList->list[0].name
               << std::endl;
        }
    }
    
    handle = beagleCreateInstance( b_tipCount               // tips
                                 , b_partialsBufferCount    // partials
                                 , b_compactBufferCount     // sequences
                                 , b_stateCount             // states
                                 , b_patternCount           // patterns (total across all subsets that use this instance)
                                 , b_eigenBufferCount       // subset models (one for each distinct eigen decomposition)
                                 , b_matrixBufferCount      // transition matrices (one for each node in each subset)
                                 , b_categoryCount          // rate categories
                                 , b_scaleBufferCount       // scale buffers
                                 , b_resourceList           // resource restrictions
                                 , b_resourceCount          // length of resource list
                                 , b_preferenceFlags        // preferred flags
                                 , b_requirementFlags       // required flags
                                 , &b_return_info           // pointer for details
                                 );

    
//    int inst = beagleCreateInstance(
//         _ntaxa,                           // tips
//         npartials,                        // partials
//         nsequences,                       // sequences
//         nstates,                          // states
//         num_patterns,                     // patterns (total across all subsets that use this instance)
//         num_subsets,                      // models (one for each distinct eigen decomposition)
//         num_subsets*num_transition_probs, // transition matrices (one for each node in each subset)
//         ngammacat,                        // rate categories
//         0,                                // scale buffers
//         NULL,                             // resource restrictions
//         0,                                // length of resource list
//         preferenceFlags,                  // preferred flags
//         requirementFlags,                 // required flags
//         &instance_details);               // pointer for details
    
    
    RBOUT ( "Initialized BEAGLE " + std::to_string(handle) );

    if ( handle < 0 )
    {
        ss << "Failed to start BEAGLE instance. "
           << "Reverting to RevBayes likelihood calculator." << std::endl;
        RbSettings::userSettings().setUseBeagle(false);
    }
    else
    {
#if defined ( RB_BEAGLE_INFO )

        ss << "BEAGLE resource "        << b_return_info.resourceNumber     << std::endl;
        ss << "\t" << "Rsrc Name : "    << b_return_info.resourceName       << std::endl;
        ss << "\t" << "Impl Name : "    << b_return_info.implName           << std::endl;
        ss << "\t" << "Flags :";
        ss << BeagleUtilities::printBeagleFlags(b_return_info.flags)        << std::endl;

        ss << "BEAGLE instance "               << handle                    << std::endl;
        ss << "\t" << "tipCount            : " << b_tipCount                << std::endl;
        ss << "\t" << "partialsBufferCount : " << b_partialsBufferCount     << std::endl;
        ss << "\t" << "compactBufferCount  : " << b_compactBufferCount      << std::endl;
        ss << "\t" << "stateCount          : " << b_stateCount              << std::endl;
        ss << "\t" << "patternCount        : " << b_patternCount            << std::endl;
        ss << "\t" << "eigenBufferCount    : " << b_eigenBufferCount        << std::endl;
        ss << "\t" << "matrixBufferCount   : " << b_matrixBufferCount       << std::endl;
        ss << "\t" << "categoryCount       : " << b_categoryCount           << std::endl;
        ss << "\t" << "scaleBufferCount    : " << b_scaleBufferCount        << std::endl;
        ss << "\t" << "resource            : " << b_resource                << std::endl;
        ss << std::endl;
#endif /* RB_BEAGLE_INFO */

        if ( b_use_cpu_threading )
        {
            this->setCPUThreadCount( RbSettings::userSettings().getBeagleMaxCPUThreads() );
        }

        ss << std::endl;

    }

    ss << std::endl;
    RBOUT(ss.str());
}



void BeagleInstance::freeBEAGLE( void )
{
    
    RBOUT ( "Finalizing BEAGLE " + std::to_string(handle) );
    // finalize BEAGLE
    int code = beagleFinalizeInstance( handle );
    if (code != 0)
    {
        std::stringstream ss;
        ss << "Likelihood failed to finalize BeagleLib instance. BeagleLib error code was ";
        ss << code << " (" << BeagleUtilities::printErrorCode(code) <<  ").";
        
        throw RbException( ss.str() );
    }
    
    // reset internal flags
    handle              = -1;
    resourcenumber      = -1;
    resourcename        = "";
    nstates             = 0;
    nratecateg          = 0;
    npatterns           = 0;
    partial_offset      = 0;
    tmatrix_offset      = 0;
}



void BeagleInstance::setCPUThreadCount(size_t n)
{
    nthreads = n;
    beagleSetCPUThreadCount( handle, n);
}

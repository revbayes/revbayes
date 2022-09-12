#if defined (RB_BEAGLE)

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
        handle         = -1;
        resourcenumber = -1;
        resourcename   = "";
        nstates        = 0;
        nratecateg     = 0;
        npatterns      = 0;
        partial_offset = 0;
        tmatrix_offset = 0;
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

    size_t b_max_cpu_threads  = RbSettings::userSettings().getBeagleMaxCPUThreads();

    int    b_resource         = 0;
    int  * b_resourceList     = &b_resource;
    int    b_resourceCount    = 1;

    long  b_requirementFlags  = 0;  //-- We dont make any hard requirements, only soft preferences.
    long  b_preferenceFlags   = this->generateSettingsBitmap();
    //long  b_requirementFlags  = this->generateSettingsBitmap();
    //long  b_preferenceFlags   = 0; 
    //long  b_requirementFlags  = this->generateSettingsBitmap();
    //long  b_preferenceFlags   = this->generateSettingsBitmap();

    BeagleInstanceDetails b_return_info;

    if (RbSettings::userSettings().getBeagleDevice() == "auto")
    {
#if defined ( RB_BEAGLE_INFO )
	ss << "Running benchmarks to automatically select fastest BEAGLE resource... ";
	ss << std::endl;
#endif 
    
	BeagleBenchmarkedResourceList * rBList;
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

    b_resource = (rBList != NULL) ? rBList->list[0].number : 0;

#if defined ( RB_BEAGLE_INFO )
    if ( rBList != NULL ) {
	ss << "\tUsing resource " << rBList->list[0].number << ": " << rBList->list[0].name;
	if ( rBList->list[0].number != 0 ) {
	    ss << " (" << rBList->list[0].performanceRatio << "x CPU)";
	}
	ss << std::endl;
    } else {
	ss << "\tResource benchmarking failed, using resource "
	   << b_resource << ": " << rBList->list[0].name
	   << std::endl;
    }
#endif //-- RB_BEAGLE_INFO
    }
    else
    {
        b_resource = RbSettings::userSettings().getBeagleResource();
    }

    handle = beagleCreateInstance( b_tipCount              // tips
                                 , b_partialsBufferCount   // partials
                                 , b_compactBufferCount    // sequences
                                 , b_stateCount            // states
                                 , b_patternCount          // patterns (total across all subsets that use this instance)
                                 , b_eigenBufferCount      // subset models (one for each distinct eigen decomposition)
                                 , b_matrixBufferCount     // transition matrices (one for each node in each subset)
                                 , b_categoryCount         // rate categories
                                 , b_scaleBufferCount      // scale buffers
                                 , b_resourceList          // resource restrictions
                                 , b_resourceCount         // length of resource list
                                 , b_preferenceFlags       // preferred flags
                                 , b_requirementFlags      // required flags
                                 , &b_return_info          // pointer for details
                                 );

    if ( handle < 0 ) {
        RbSettings::userSettings().setUseBeagle(false);
        ss << "Failed to start BEAGLE instance. "
           << "Reverting to RevBayes likelihood calculator." << std::endl;
    }
    else {
	//--- Remaining beagle post-configuration ---

	//-- Set threads count
    	this->setCPUThreadCount( b_max_cpu_threads );
	
#if defined ( RB_BEAGLE_INFO )
	ss << "\t" << "resource            : " << b_return_info.resourceNumber << std::endl;
	ss << "\t" << "Rsrc Name           : " << b_return_info.resourceName   << std::endl;
	ss << "\t" << "Impl Name           : " << b_return_info.implName       << std::endl;
	ss << "\t" << "tipCount            : " << b_tipCount                   << std::endl;
	ss << "\t" << "partialsBufferCount : " << b_partialsBufferCount        << std::endl;
	ss << "\t" << "compactBufferCount  : " << b_compactBufferCount         << std::endl;
	ss << "\t" << "stateCount          : " << b_stateCount                 << std::endl;
	ss << "\t" << "patternCount        : " << b_patternCount               << std::endl;
	ss << "\t" << "eigenBufferCount    : " << b_eigenBufferCount           << std::endl;
	ss << "\t" << "matrixBufferCount   : " << b_matrixBufferCount          << std::endl;
	ss << "\t" << "categoryCount       : " << b_categoryCount              << std::endl;
	ss << "\t" << "scaleBufferCount    : " << b_scaleBufferCount           << std::endl;
        if ( b_max_cpu_threads > 1 ) {
            ss << "\t" << "max_threads         : " << b_max_cpu_threads << std::endl;
        }
#endif /* RB_BEAGLE_INFO */
    }

    RBOUT(ss.str());
}


void BeagleInstance::freeBEAGLE( void )
{
#if defined ( RB_BEAGLE_INFO )
    std::stringstream ss;
    ss << "Finalizing BEAGLE " << std::to_string(handle) << std::endl << std::endl;
#endif /* RB_BEAGLE_INFO */

    // finalize BEAGLE
    int code = beagleFinalizeInstance( handle );
    if (code != 0)
    {
        std::stringstream err_ss;
        err_ss << "Likelihood failed to finalize BeagleLib instance. BeagleLib error code was ";
        err_ss << code << " (" << BeagleUtilities::printErrorCode(code) <<  ").";
        throw RbException( err_ss.str() );
    }
    
    // reset internal flags
    handle          = -1;
    resourcenumber  = -1;
    resourcename    = "";
    nstates         = 0;
    nratecateg      = 0;
    npatterns       = 0;
    partial_offset  = 0;
    tmatrix_offset  = 0;

#if defined ( RB_BEAGLE_INFO )
    RBOUT(ss.str());
#endif /* RB_BEAGLE_INFO */
}


long BeagleInstance::generateSettingsBitmap(void)
{
    //-- Read in user settings
    bool              b_use_double      = RbSettings::userSettings().getBeagleUseDoublePrecision();
    size_t            b_max_cpu_threads = RbSettings::userSettings().getBeagleMaxCPUThreads();
    const std::string b_scaling_mode    = RbSettings::userSettings().getBeagleScalingMode();
    const std::string b_device          = RbSettings::userSettings().getBeagleDevice();

    //-- Beagle flags that we will configure
    long b_flag_device        = 0;
    long b_flag_precision     = 0;
    long b_flag_vectorization = 0;
    long b_flag_threading     = 0;
    long b_flag_scaling       = 0;
    long b_flag_gpu_ops       = 0;

    //-- Set flags related to each device type
    if (b_device == "cpu") {
    	b_flag_device = BEAGLE_FLAG_PROCESSOR_CPU | BEAGLE_FLAG_FRAMEWORK_CPU;
    	b_flag_vectorization = BEAGLE_FLAG_VECTOR_NONE;
    } else if (b_device == "cpu_sse") {
    	b_flag_device = BEAGLE_FLAG_PROCESSOR_CPU | BEAGLE_FLAG_FRAMEWORK_CPU;
    	b_flag_vectorization = BEAGLE_FLAG_VECTOR_SSE;
    } else if (b_device == "cpu_avx") {
    	b_flag_device = BEAGLE_FLAG_PROCESSOR_CPU | BEAGLE_FLAG_FRAMEWORK_CPU;
    	b_flag_vectorization = BEAGLE_FLAG_VECTOR_AVX;
    } else if (b_device == "gpu_opencl") {
    	b_flag_device = BEAGLE_FLAG_PROCESSOR_GPU | BEAGLE_FLAG_FRAMEWORK_OPENCL;
        b_flag_gpu_ops = BEAGLE_FLAG_PARALLELOPS_GRID | BEAGLE_FLAG_COMPUTATION_SYNCH;
    } else if (b_device == "gpu_cuda") {
    	b_flag_device = BEAGLE_FLAG_PROCESSOR_GPU | BEAGLE_FLAG_FRAMEWORK_CUDA;
        b_flag_gpu_ops = BEAGLE_FLAG_PARALLELOPS_GRID;
    }

    // Configure floating point precision.
    if (b_use_double) {
    	b_flag_precision = BEAGLE_FLAG_PRECISION_DOUBLE;
    } else {
    	b_flag_precision = BEAGLE_FLAG_PRECISION_SINGLE;
    }

    // Configure scaling. We error check when the passing the setting from revlang,
    // validity check not needed here.
    if (b_scaling_mode == "dynamic") {
	b_flag_scaling = BEAGLE_FLAG_SCALING_DYNAMIC;
    } else if (b_scaling_mode == "auto") {
        b_flag_scaling = BEAGLE_FLAG_SCALING_AUTO;    
    } else if (b_scaling_mode == "always") {
	b_flag_scaling = BEAGLE_FLAG_SCALING_ALWAYS; 
    } else if (b_scaling_mode == "manual") {
	b_flag_scaling = BEAGLE_FLAG_SCALING_MANUAL;
    }

    // Configure threading. We error check when the passing the setting from revlang,
    // validity check not needed here.
    if (b_device == "cpu" || b_device == "cpu_sse" || b_device == "cpu_avx") {
      if (b_max_cpu_threads != 1) {
        b_flag_threading = BEAGLE_FLAG_THREADING_CPP;
      } else {
    	b_flag_threading = BEAGLE_FLAG_THREADING_NONE;
      }
    }

    return ( b_flag_device
	   | b_flag_precision 
	   | b_flag_vectorization 
	   | b_flag_threading 
	   | b_flag_scaling
	   | b_flag_gpu_ops 
	   );
}


void BeagleInstance::setCPUThreadCount(size_t n)
{
    nthreads = n;
    beagleSetCPUThreadCount( handle, n );
}

#endif


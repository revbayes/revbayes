/**
 * @file
 * This file contains implementations of helper functions for BEAGLE.
 *
 * @brief Implementation of helper functions for the BEAGLE library.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author Daniel L. Ayres
 * @license GPL version 3
 *
 * $Id$
 */

#if defined( RB_BEAGLE )

#include "BeagleUtilities.h"

#include <sstream>
#include <vector>

#include "libhmsbeagle/beagle.h"

/** outputs BEAGLE capability flags */
std::string RevBayesCore::BeagleUtilities::printBeagleFlags( long inFlags )
{
    std::vector<std::string> printFlags;

    if (inFlags & BEAGLE_FLAG_PRECISION_SINGLE   ) printFlags.push_back( "PRECISION_SINGLE"   );
    if (inFlags & BEAGLE_FLAG_PRECISION_DOUBLE   ) printFlags.push_back( "PRECISION_DOUBLE"   );
    if (inFlags & BEAGLE_FLAG_COMPUTATION_SYNCH  ) printFlags.push_back( "COMPUTATION_SYNCH"  );
    if (inFlags & BEAGLE_FLAG_COMPUTATION_ASYNCH ) printFlags.push_back( "COMPUTATION_ASYNCH" );
    if (inFlags & BEAGLE_FLAG_EIGEN_REAL         ) printFlags.push_back( "EIGEN_REAL"         );
    if (inFlags & BEAGLE_FLAG_EIGEN_COMPLEX      ) printFlags.push_back( "EIGEN_COMPLEX"      );
    if (inFlags & BEAGLE_FLAG_SCALING_MANUAL     ) printFlags.push_back( "SCALING_MANUAL"     );
    if (inFlags & BEAGLE_FLAG_SCALING_AUTO       ) printFlags.push_back( "SCALING_AUTO"       );
    if (inFlags & BEAGLE_FLAG_SCALING_ALWAYS     ) printFlags.push_back( "SCALING_ALWAYS"     );
    if (inFlags & BEAGLE_FLAG_SCALING_DYNAMIC    ) printFlags.push_back( "SCALING_DYNAMIC"    );
    if (inFlags & BEAGLE_FLAG_SCALERS_RAW        ) printFlags.push_back( "SCALERS_RAW"        );
    if (inFlags & BEAGLE_FLAG_SCALERS_LOG        ) printFlags.push_back( "SCALERS_LOG"        );
    if (inFlags & BEAGLE_FLAG_INVEVEC_STANDARD   ) printFlags.push_back( "INVEVEC_STANDARD"   );
    if (inFlags & BEAGLE_FLAG_INVEVEC_TRANSPOSED ) printFlags.push_back( "INVEVEC_TRANSPOSED" );
    if (inFlags & BEAGLE_FLAG_VECTOR_SSE         ) printFlags.push_back( "VECTOR_SSE"         );
    if (inFlags & BEAGLE_FLAG_VECTOR_AVX         ) printFlags.push_back( "VECTOR_AVX"         );
    if (inFlags & BEAGLE_FLAG_VECTOR_NONE        ) printFlags.push_back( "VECTOR_NONE"        );
    if (inFlags & BEAGLE_FLAG_THREADING_CPP      ) printFlags.push_back( "THREADING_CPP"      );
    if (inFlags & BEAGLE_FLAG_THREADING_OPENMP   ) printFlags.push_back( "THREADING_OPENMP"   );
    if (inFlags & BEAGLE_FLAG_THREADING_NONE     ) printFlags.push_back( "THREADING_NONE"     );
    if (inFlags & BEAGLE_FLAG_PROCESSOR_CPU      ) printFlags.push_back( "PROCESSOR_CPU"      );
    if (inFlags & BEAGLE_FLAG_PROCESSOR_GPU      ) printFlags.push_back( "PROCESSOR_GPU"      );
    if (inFlags & BEAGLE_FLAG_PROCESSOR_FPGA     ) printFlags.push_back( "PROCESSOR_FPGA"     );
    if (inFlags & BEAGLE_FLAG_PROCESSOR_CELL     ) printFlags.push_back( "PROCESSOR_CELL"     );
    if (inFlags & BEAGLE_FLAG_PROCESSOR_PHI      ) printFlags.push_back( "PROCESSOR_PHI"      );
    if (inFlags & BEAGLE_FLAG_PROCESSOR_OTHER    ) printFlags.push_back( "PROCESSOR_OTHER"    );
    if (inFlags & BEAGLE_FLAG_FRAMEWORK_CUDA     ) printFlags.push_back( "FRAMEWORK_CUDA"     );
    if (inFlags & BEAGLE_FLAG_FRAMEWORK_OPENCL   ) printFlags.push_back( "FRAMEWORK_OPENCL"   );
    if (inFlags & BEAGLE_FLAG_FRAMEWORK_CPU      ) printFlags.push_back( "FRAMEWORK_CPU"      );
    if (inFlags & BEAGLE_FLAG_PARALLELOPS_STREAMS) printFlags.push_back( "PARALLELOPS_STREAMS");
    if (inFlags & BEAGLE_FLAG_PARALLELOPS_GRID   ) printFlags.push_back( "PARALLELOPS_GRID"   );

    std::stringstream ss;

    int flagCount = 0;
    for ( std::vector<std::string>::iterator itFlag = printFlags.begin(); itFlag != printFlags.end(); ++itFlag)
    {
        if ( flagCount % 6 == 0 && flagCount > 0 )
        {
            ss << std::endl << "\t      ";
        }
        ss << " " << *itFlag;
        flagCount++;
    }

    return ss.str();
}


std::string RevBayesCore::BeagleUtilities::printErrorCode(int code )
{
    
    switch( code )
    {
        case  0 : return "success";
        case -1 : return "unspecified error";
        case -2 : return "not enough memory could be allocated";
        case -3 : return "unspecified exception";
        case -4 : return "the instance index is out of range, or the instance has not been created";
        case -5 : return "one of the indices specified exceeded the range of the array";
        case -6 : return "no resource matches requirements";
        case -7 : return "no implementation matches requirements";
        case -8 : return "floating-point range exceeded";
    }

    return "Unknown BEAGLE error code";
}


#endif /* RB_BEAGLE */

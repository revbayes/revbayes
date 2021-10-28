#include <iostream>

#include "RbSettings.h"
#include "RlUserInterface.h"
#include "BeagleUtilities.h"
#include "libhmsbeagle/beagle.h"

#include "BeagleInstance.h"

using namespace RevBayesCore;


BeagleInstance *BeagleInstance::beagle_singleton = NULL;


BeagleInstance *
BeagleInstance::getInstance ( int  b_resource
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
                            )
{
    if ( !beagle_singleton )
    {
        beagle_singleton = new BeagleInstance( b_resource
                                             , b_use_cpu_threading
                                             , b_use_scaling
                                             , b_tipCount
                                             , b_partialsBufferCount
                                             , b_compactBufferCount
                                             , b_stateCount
                                             , b_patternCount
                                             , b_eigenBufferCount
                                             , b_matrixBufferCount
                                             , b_categoryCount
                                             , b_scaleBufferCount
                                             );
    }
    return beagle_singleton;
}


BeagleInstance::BeagleInstance ( int  b_resource
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
                               )
{

    std::stringstream ss;

    if ( RbSettings::userSettings().getUseBeagle() == true )
    {
        ss << std::endl;
        ss << "Using BEAGLE library v" << beagleGetVersion();
        ss << " for parallel likelihood evaluation (https://beagle-dev.github.io/)";
        ss << std::endl;

        //int  b_resource            = RbSettings::userSettings().getBeagleResource();
        //bool b_use_cpu_threading   = RbSettings::userSettings().getBeagleMaxCPUThreads() != 1
        //                           ? true : false;
        //bool b_use_scaling         = RbSettings::userSettings().getBeagleScalingMode() != "none"
        //                           ? true : false;

        BeagleInstanceDetails b_return_info;

        int* b_resourceList        = &b_resource;
        int  b_resourceCount       = 1;
        long b_requirementFlags    = 0;
        long b_preferenceFlags     = ( RbSettings::userSettings().getBeagleUseDoublePrecision()
                                     ? BEAGLE_FLAG_PRECISION_DOUBLE
                                     : BEAGLE_FLAG_PRECISION_SINGLE
                                     )
                                   | ( b_use_cpu_threading
                                     ? BEAGLE_FLAG_THREADING_CPP
                                     : 0
                                     );


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
                ss << "Using resource " << rBList->list[0].number << ": " << rBList->list[0].name;
                if ( rBList->list[0].number != 0 )
                {
                    ss << " (" << rBList->list[0].performanceRatio << "x CPU)";
                }
                ss << std::endl;
            } else {
                ss << "Resource benchmarking failed, using resource "
                   << b_resource << ": " << rBList->list[0].name
                   << std::endl;
            }
        }

        this->resourceID = beagleCreateInstance( b_tipCount
                                               , b_partialsBufferCount
                                               , b_compactBufferCount
                                               , b_stateCount
                                               , b_patternCount
                                               , b_eigenBufferCount
                                               , b_matrixBufferCount
                                               , b_categoryCount
                                               , b_scaleBufferCount
                                               , b_resourceList
                                               , b_resourceCount
                                               , b_preferenceFlags
                                               , b_requirementFlags
                                               , &b_return_info
                                               );

        #if defined ( RB_BEAGLE_INFO )
            ss << "BEAGLE parameters"                                       << std::endl;
            ss << "\t" << "tipCount            : " << b_tipCount            << std::endl;
            ss << "\t" << "partialsBufferCount : " << b_partialsBufferCount << std::endl;
            ss << "\t" << "compactBufferCount  : " << b_compactBufferCount  << std::endl;
            ss << "\t" << "stateCount          : " << b_stateCount          << std::endl;
            ss << "\t" << "patternCount        : " << b_patternCount        << std::endl;
            ss << "\t" << "eigenBufferCount    : " << b_eigenBufferCount    << std::endl;
            ss << "\t" << "matrixBufferCount   : " << b_matrixBufferCount   << std::endl;
            ss << "\t" << "categoryCount       : " << b_categoryCount       << std::endl;
            ss << "\t" << "scaleBufferCount    : " << b_scaleBufferCount    << std::endl;
            ss << "\t" << "resource            : " << b_resource            << std::endl;
            ss << std::endl;
            ss << "BEAGLE instance: " << this->resourceID << std::endl;
        #endif /* RB_BEAGLE_INFO */

        if ( this->resourceID < 0 )
        {
            ss << "Failed to start BEAGLE instance. "
               << "Reverting to RevBayes likelihood calculator." << std::endl;
            RbSettings::userSettings().setUseBeagle(false);
        }
        else
        {
            ss << "Using BEAGLE resource " << b_return_info.resourceNumber << std::endl;
            ss << "\t" << "Rsrc Name : "   << b_return_info.resourceName   << std::endl;
            ss << "\t" << "Impl Name : "   << b_return_info.implName       << std::endl;
            ss << "\t" << "Flags :";
            ss << BeagleUtilities::printBeagleFlags(b_return_info.flags);
            ss << std::endl;

            if ( b_use_cpu_threading )
            {
                beagleSetCPUThreadCount( this->resourceID
                                       , RbSettings::userSettings().getBeagleMaxCPUThreads()
                                       );
            }

        }
        ss << std::endl;

    }

    ss << std::endl;
    RBOUT(ss.str());
}


BeagleInstance::~BeagleInstance ( )
{
    RBOUT ( "Finalizing BEAGLE" );
    beagleFinalizeInstance(beagle_singleton->getResourceID());
    BeagleInstance::beagle_singleton = NULL;
    //this->resourceID = -1;
}

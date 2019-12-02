#ifndef AddContinuousDistribution_h
#define AddContinuousDistribution_h

#include "RlContinuousDistribution.h"
#include "DistributionFunctionCdf.h"
#include "DistributionFunctionPdf.h"
#include "DistributionFunctionQuantileContinuous.h"
#include "DistributionFunctionQuantilePositiveContinuous.h"
#include "DistributionFunctionRv.h"
#include "RlPositiveContinuousDistribution.h"
#include "RlTypedDistribution.h"
#include "StringUtilities.h"
#include "Workspace.h"

/**
 * Generic function
 * @param value
 * @return
 */
template <typename T>
void AddContinuousDistribution( RevLanguage::ContinuousDistribution *dist)
{
    
    RevLanguage::Workspace::globalWorkspace().addDistribution( dist );
    RevLanguage::Workspace::globalWorkspace().addFunction( new DistributionFunctionPdf<T>( dist->clone() ) );
    RevLanguage::Workspace::globalWorkspace().addFunction( new DistributionFunctionCdf( dist->clone() ) );
    RevLanguage::Workspace::globalWorkspace().addFunction( new DistributionFunctionQuantileContinuous( dist->clone() ) );
    RevLanguage::Workspace::globalWorkspace().addFunction( new DistributionFunctionRv<T>( dist->clone() ) );
    
}


/**
 * Generic function
 * @param dist The distribution that we are going to add
 */
template <typename T>
void AddContinuousDistribution( RevLanguage::PositiveContinuousDistribution *dist)
{
    
    RevLanguage::Workspace::globalWorkspace().addDistribution( dist );
    RevLanguage::Workspace::globalWorkspace().addFunction( new DistributionFunctionPdf<T>( dist->clone() ) );
    RevLanguage::Workspace::globalWorkspace().addFunction( new DistributionFunctionCdf( dist->clone() ) );
    RevLanguage::Workspace::globalWorkspace().addFunction( new DistributionFunctionQuantilePositiveContinuous( dist->clone() ) );
    RevLanguage::Workspace::globalWorkspace().addFunction( new DistributionFunctionRv<T>( dist->clone() ) );
    
}


#endif

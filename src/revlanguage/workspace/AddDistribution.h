#ifndef AddDistribution_h
#define AddDistribution_h

#include "DistributionFunctionPdf.h"
#include "DistributionFunctionRv.h"
#include "RlTypedDistribution.h"
#include "StringUtilities.h"
#include "Workspace.h"

/**
 * Generic function
 * @param value
 * @return
 */
template <typename T>
void AddDistribution( RevLanguage::TypedDistribution<T> *dist )
{

    RevLanguage::Workspace::globalWorkspace().addDistribution( dist );
    RevLanguage::Workspace::globalWorkspace().addFunction( new DistributionFunctionPdf<T>( dist->clone() ) );
    RevLanguage::Workspace::globalWorkspace().addFunction( new DistributionFunctionRv<T>( dist->clone() ) );

}


#endif

#ifndef DistributionTypeUtilities_H
#define DistributionTypeUtilities_H

#include "Distribution.h"
#include "RbException.h"
#include "TypedDistribution.h"

namespace RevBayesCore {

    /** Recover a distribution's known value type, checking the invariant in debug builds. */
    template <typename valueType>
    TypedDistribution<valueType>* assumeDistributionOf(Distribution *distribution)
    {
#ifndef NDEBUG
        TypedDistribution<valueType> *typedDistribution = dynamic_cast<TypedDistribution<valueType>*>( distribution );
        if ( typedDistribution == NULL )
        {
            throw RbException("A distribution does not produce the expected value type.");
        }
        return typedDistribution;
#else
        return static_cast<TypedDistribution<valueType>*>( distribution );
#endif
    }

    /** Recover a distribution's known value type, checking the invariant in debug builds. */
    template <typename valueType>
    const TypedDistribution<valueType>* assumeDistributionOf(const Distribution *distribution)
    {
#ifndef NDEBUG
        const TypedDistribution<valueType> *typedDistribution = dynamic_cast<const TypedDistribution<valueType>*>( distribution );
        if ( typedDistribution == NULL )
        {
            throw RbException("A distribution does not produce the expected value type.");
        }
        return typedDistribution;
#else
        return static_cast<const TypedDistribution<valueType>*>( distribution );
#endif
    }

}

#endif

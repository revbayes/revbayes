#ifndef DistributionTypeUtilities_H
#define DistributionTypeUtilities_H

#include "Distribution.h"
#include "RbException.h"
#include "TypedDistribution.h"

namespace RevBayesCore {

    /** Recover TypedDistribution<valueType>, checking the invariant in debug builds.
     * The object must inherit from this exact specialization; value-type inheritance is insufficient.
     */
    template <typename valueType>
    TypedDistribution<valueType>* assumeTypedDistribution(Distribution *distribution)
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

    /** Recover TypedDistribution<valueType>, checking the invariant in debug builds.
     * The object must inherit from this exact specialization; value-type inheritance is insufficient.
     */
    template <typename valueType>
    const TypedDistribution<valueType>* assumeTypedDistribution(const Distribution *distribution)
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

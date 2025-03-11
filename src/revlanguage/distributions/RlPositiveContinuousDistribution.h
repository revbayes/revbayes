#ifndef RlPositiveContinuousDistribution_H
#define RlPositiveContinuousDistribution_H

#include <cmath>
#include <iosfwd>
#include <string>
#include <vector>

#include "ContinuousDistribution.h"
#include "RealPos.h"
#include "RlTypedDistribution.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RevPtr.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlStochasticNode.h"
#include "RlTypedFunction.h"
#include "StochasticNode.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

namespace RevLanguage {
class TypeSpec;
    
    class PositiveContinuousDistribution : public TypedDistribution<RealPos> {
        
    public:
        virtual                                         ~PositiveContinuousDistribution(void);                                                                  //!< Destructor
        PositiveContinuousDistribution(const PositiveContinuousDistribution &x);                                                                //!< Copy constuctor
        
        virtual RealPos*                                createRandomVariable(void) const;                                                   //!< Create a random variable from this distribution        
        
        // Basic utility functions you have to override
        virtual PositiveContinuousDistribution*         clone(void) const = 0;                                                              //!< Clone object
        static const std::string&                       getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                             //!< Get class type spec
        
        
        // Distribution functions you have to override
        virtual RevBayesCore::ContinuousDistribution*   createDistribution(void) const = 0;                                                 //!< Create a random variable from this distribution
        
        
    protected:
        PositiveContinuousDistribution(void);                                                                                                 //!< Basic constructor
        
    };
    
}

#endif


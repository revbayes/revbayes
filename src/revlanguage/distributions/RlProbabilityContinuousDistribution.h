#ifndef RlProbabilityContinuousDistribution_H
#define RlProbabilityContinuousDistribution_H

#include <math.h>
#include <iosfwd>
#include <string>
#include <vector>

#include "ContinuousDistribution.h"
#include "Probability.h"
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
    
    class ProbabilityContinuousDistribution : public TypedDistribution<Probability> {
        
    public:
        virtual                                         ~ProbabilityContinuousDistribution(void);                                                                  //!< Destructor
        ProbabilityContinuousDistribution(const ProbabilityContinuousDistribution &x);                                                      //!< Copy constuctor
        
        virtual Probability*                            createRandomVariable(void) const;                                                   //!< Create a random variable from this distribution
        
        // Basic utility functions you have to override
        virtual ProbabilityContinuousDistribution*      clone(void) const = 0;                                                              //!< Clone object
        static const std::string&                       getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                             //!< Get class type spec
        
        
        // Distribution functions you have to override
        virtual RevBayesCore::ContinuousDistribution*   createDistribution(void) const = 0;                                                 //!< Create a random variable from this distribution
        
        
    protected:
        ProbabilityContinuousDistribution(void);                                                                                            //!< Basic constructor
        
    };
    
}

#endif

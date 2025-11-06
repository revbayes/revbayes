#ifndef Dist_phyloDistanceGamma_H
#define Dist_phyloDistanceGamma_H

#include <cmath>
#include <iosfwd>
#include <string>
#include <vector>

#include "RlTypedDistribution.h"
#include "RlDistanceMatrix.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DistanceMatrix.h"
#include "DynamicNode.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
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
    
    class Dist_phyloDistanceGamma : public TypedDistribution< DistanceMatrix > {
        
    public:
        Dist_phyloDistanceGamma( void );
        virtual ~Dist_phyloDistanceGamma();
        
        // Basic utility functions
        Dist_phyloDistanceGamma*                        clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::TypedDistribution< RevBayesCore::DistanceMatrix >*      createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       tree;
        RevPtr<const RevVariable>                       varianceMatrix;
        RevPtr<const RevVariable>                       distanceMatrix;
		RevPtr<const RevVariable>                       names;
		

		
    };
    
}


#endif

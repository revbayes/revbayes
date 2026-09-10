#ifndef Dist_PhyloMultivariateBrownianREML_H
#define Dist_PhyloMultivariateBrownianREML_H

#include <cmath>
#include <iosfwd>
#include <string>
#include <vector>

#include "ContinuousCharacterData.h"
#include "RlContinuousCharacterData.h"
#include "RlTypedDistribution.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
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

namespace RevBayesCore { class RbHelpReference; }

namespace RevLanguage {
class TypeSpec;
    
    class Dist_PhyloMultivariateBrownianREML :  public TypedDistribution< ContinuousCharacterData > {
        
    public:
        Dist_PhyloMultivariateBrownianREML( void );
        virtual ~Dist_PhyloMultivariateBrownianREML();
        
        // Basic utility functions
        Dist_PhyloMultivariateBrownianREML*             clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                        getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::TypedDistribution< RevBayesCore::ContinuousCharacterData >*      createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       tree;
        RevPtr<const RevVariable>                       branchRates;
        RevPtr<const RevVariable>                       site_rates;
//        RevPtr<const RevVariable>                       nSites;
        RevPtr<const RevVariable>                       rate_matrix;
        
        
    };
    
}

#endif

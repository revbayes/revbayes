#ifndef Dist_phyloCTMCDASequence_H
#define Dist_phyloCTMCDASequence_H

#include <cmath>
#include <iosfwd>
#include <string>
#include <vector>

#include "AbstractHomologousDiscreteCharacterData.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlTypedDistribution.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "MethodTable.h"
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
    
    class Dist_phyloCTMCDASequence : public TypedDistribution< AbstractHomologousDiscreteCharacterData > {
        
    public:
        Dist_phyloCTMCDASequence( void );
        virtual ~Dist_phyloCTMCDASequence();
        
        // Basic utility functions
        Dist_phyloCTMCDASequence*                        clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        virtual MethodTable                             getDistributionMethods( void ) const;                                                                       //!< Get the member methods
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::TypedDistribution< RevBayesCore::AbstractHomologousDiscreteCharacterData >*      createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       tree;
        RevPtr<const RevVariable>                       q;
        RevPtr<const RevVariable>                       clado_probs;
        RevPtr<const RevVariable>                       rate;
        RevPtr<const RevVariable>                       siteRates;
        RevPtr<const RevVariable>                       rootFrequencies;
        RevPtr<const RevVariable>                       pInv;
        RevPtr<const RevVariable>                       nSites;
        RevPtr<const RevVariable>                       type;
        RevPtr<const RevVariable>                       storeInternalNodes;
        RevPtr<const RevVariable>                       treatAmbiguousAsGap;
        RevPtr<const RevVariable>                       coding;
        
        
    };
    
}


#endif

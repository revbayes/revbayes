#ifndef Dist_PhyloMultivariateBrownianMultiSampleREML_H
#define Dist_PhyloMultivariateBrownianMultiSampleREML_H

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

namespace RevLanguage {
class TypeSpec;
    
    class Dist_PhyloMultivariateBrownianMultiSampleREML : public TypedDistribution< ContinuousCharacterData > {
        
    public:
        Dist_PhyloMultivariateBrownianMultiSampleREML( void );
        virtual ~Dist_PhyloMultivariateBrownianMultiSampleREML();
        
        // Basic utility functions
        Dist_PhyloMultivariateBrownianMultiSampleREML*                                  clone(void) const;                                                                      //!< Clone the object
        static const std::string&                                                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                                                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::TypedDistribution< RevBayesCore::ContinuousCharacterData >*       createDistribution(void) const;
        
    protected:
        
        void                                                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                                                       tree;
        RevPtr<const RevVariable>                                                       branch_rates;
        RevPtr<const RevVariable>                                                       within_species_variances;
        RevPtr<const RevVariable>                                                       within_species_variances2;
        RevPtr<const RevVariable>                                                       rate_matrix;
        RevPtr<const RevVariable>                                                       taxa;                                                                                   //!< The taxons
        
        
    };
    
}

#endif

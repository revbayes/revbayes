#ifndef Dist_PhyloOrnsteinUhlenbeckThreePoint_H
#define Dist_PhyloOrnsteinUhlenbeckThreePoint_H

#include "ContinuousCharacterData.h"
#include "RlContinuousCharacterData.h"
#include "RlTypedDistribution.h"
#include "Tree.h"

namespace RevLanguage {
    
    class Dist_PhyloOrnsteinUhlenbeckThreePoint :  public TypedDistribution< ContinuousCharacterData > {
        
    public:
        Dist_PhyloOrnsteinUhlenbeckThreePoint( void );
        virtual ~Dist_PhyloOrnsteinUhlenbeckThreePoint();
        
        // Basic utility functions
        Dist_PhyloOrnsteinUhlenbeckThreePoint*          clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
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
        RevPtr<const RevVariable>                       alpha;
        RevPtr<const RevVariable>                       theta;
        RevPtr<const RevVariable>                       sigma;
        RevPtr<const RevVariable>                       root_states;
        RevPtr<const RevVariable>                       n_sites;
        
        
    };
    
}


#endif



#ifndef Dist_PhyloNodeStateOU_H
#define Dist_PhyloNodeStateOU_H

#include <iostream>


#include "PhyloNodeStateOU.h"
#include "ModelVector.h"
#include "RlTypedDistribution.h"
#include "RealPos.h"

namespace RevLanguage {

    class Dist_PhyloNodeStateOU : public TypedDistribution< ModelVector<Real> > {
        
    public:
                                                        Dist_PhyloNodeStateOU( void ) {};
        virtual                                        ~Dist_PhyloNodeStateOU() {};
        
        // Basic utility functions
        Dist_PhyloNodeStateOU*                          clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                        getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::PhyloNodeStateOU*                 createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       tree;
        RevPtr<const RevVariable>                       root_state;
        RevPtr<const RevVariable>                       sigma;
        RevPtr<const RevVariable>                       alpha;
        RevPtr<const RevVariable>                       theta;
        
    };
    
}

#endif

#ifndef Dist_PhyloBrownian_H
#define Dist_PhyloBrownian_H

#include <iostream>


#include "PhyloBrownianProcess.h"
#include "ModelVector.h"
#include "RlTypedDistribution.h"
#include "RlSimplex.h"
#include "Real.h"

namespace RevLanguage {

    class Dist_PhyloBrownian :  public TypedDistribution< ModelVector<Real> > {
        
    public:
                                                        Dist_PhyloBrownian( void ) {};
        virtual                                        ~Dist_PhyloBrownian() {};
        
        // Basic utility functions
        Dist_PhyloBrownian*                             clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                        getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::PhyloBrownianProcess*             createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       tree;
        RevPtr<const RevVariable>                       sigma;
        RevPtr<const RevVariable>                       drift;
        
    };
    
}

#endif /* defined(__revbayes__Dist_brownian__) */

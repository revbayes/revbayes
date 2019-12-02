#ifndef Dist_PhyloMvtBrownian_H
#define Dist_PhyloMvtBrownian_H

#include <iostream>


#include "PhyloMultivariateBrownianProcess.h"
#include "RlTypedDistribution.h"
#include "MatrixReal.h"
#include "Real.h"

namespace RevLanguage {

    class Dist_PhyloMvtBrownian :  public TypedDistribution< ModelVector< ModelVector<Real> >  > {
        
    public:
                                                        Dist_PhyloMvtBrownian( void ) {};
        virtual                                        ~Dist_PhyloMvtBrownian() {};
        
        // Basic utility functions
        Dist_PhyloMvtBrownian*                          clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
                
        // Distribution functions you have to override
        RevBayesCore::PhyloMultivariateBrownianProcess* createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       tree;
        RevPtr<const RevVariable>                       sigma;
 //       RevPtr<const RevVariable>                          rootval;
        
    };
    
}

#endif /* defined(__revbayes__Dist_mvtBrownian__) */

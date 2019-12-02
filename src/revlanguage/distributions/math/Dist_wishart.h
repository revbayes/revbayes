#ifndef Dist_wishart_H
#define Dist_wishart_H

#include <iostream>

#include "RlMatrixRealSymmetric.h"
#include "RlTypedDistribution.h"
#include "WishartDistribution.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the wishart distribution.
     *
     * The RevLanguage wrapper of the wishart distribution simply
     * manages the interactions through the Rev with our core.
     * That is, the internal distribution object can be constructed and hooked up
     * in a model graph.
     * See the WishartDistribution for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Nicolas Lartillot)
     * @since 2014-03-27, version 1.0
     *
     */
    class Dist_wishart : public TypedDistribution<MatrixRealSymmetric> {
        
    public:
        Dist_wishart( void );
        virtual ~Dist_wishart();
        
        // Basic utility functions
        Dist_wishart*                                   clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::WishartDistribution*              createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
 //       RevPtr<const RevVariable>                       omega;
        RevPtr<const RevVariable>                       kappa;
        RevPtr<const RevVariable>                       df;
        RevPtr<const RevVariable>                       dim;
    };
    
}
#endif /* defined(__revbayes__Dist_wishart__) */

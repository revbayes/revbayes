#ifndef Dist_decomposedInverseWishart_H
#define	Dist_decomposedInverseWishart_H


#include <iostream>

#include "RlMatrixRealSymmetric.h"
#include "RlTypedDistribution.h"
#include "DecomposedInverseWishartDistribution.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the decomposed inverse Wishart distribution.
     *
     * The RevLanguage wrapper of the decomposed inverse Wishart distribution simply
     * manages the interactions through the Rev with our core.
     * That is, the internal distribution object can be constructed and hooked up
     * in a model graph.
     * See the DecomposedInverseWishartDistribution.h for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (John Huelsenbeck)
     * @since 2014-11-21, version 1.0
     *
     */
    class Dist_decomposedInverseWishart :  public TypedDistribution<MatrixReal> {
        
    public:
                                                            Dist_decomposedInverseWishart(void);
        virtual                                            ~Dist_decomposedInverseWishart(void);
        
        // Basic utility functions
        Dist_decomposedInverseWishart*                      clone(void) const;                                                                      //!< Clone the object
        static const std::string&                           getClassType(void);                                                                     //!< Get class name
        static const TypeSpec&                              getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                         getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                     getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                  getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                                printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        // Distribution functions you have to override
        RevBayesCore::DecomposedInverseWishartDistribution* createDistribution(void) const;
        
    protected:
        void                                                setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        RevPtr<const RevVariable>                           sigma;
        RevPtr<const RevVariable>                           diagonal;
        RevPtr<const RevVariable>                           kappa;
        RevPtr<const RevVariable>                           df;
        RevPtr<const RevVariable>                           dim;
    };
    
}

#endif	/* DIST_INVERSEWISHART_H */


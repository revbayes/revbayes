#ifndef Dist_LKJPartial_H
#define Dist_LKJPartial_H

#include <iostream>

#include "RlMatrixRealSymmetric.h"
#include "RlTypedDistribution.h"
#include "LKJPartialDistribution.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the LKJPartial distribution.
     *
     * The RevLanguage wrapper of the LKJPartial distribution simply
     * manages the interactions through the Rev with our core.
     * That is, the internal distribution object can be constructed and hooked up
     * in a model graph.
     * See the LKJPartialDistribution for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Nicolas Lartillot)
     * @since 2014-03-27, version 1.0
     *
     */
    class Dist_LKJPartial : public TypedDistribution<MatrixRealSymmetric> {
        
    public:
        Dist_LKJPartial( void );
        virtual ~Dist_LKJPartial();
        
        // Basic utility functions
        Dist_LKJPartial*                                clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                        getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::LKJPartialDistribution*           createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        RevPtr<const RevVariable>                       eta;
        RevPtr<const RevVariable>                       dim;

    };
    
}
#endif /* defined(__revbayes__Dist_LKJPartial__) */

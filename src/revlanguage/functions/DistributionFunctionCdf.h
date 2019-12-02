#ifndef DistributionFunctionCdf_H
#define DistributionFunctionCdf_H

#include "Probability.h"
#include "RlContinuousDistribution.h"
#include "RlTypedFunction.h"
#include "RlPositiveContinuousDistribution.h"

#include <string>
#include <vector>

namespace RevLanguage {
    
    
    /**
     * \brief Wraper for the cummulative distribution function.
     *
     * This class manages the function call from Rev to the cdf (cummulative distribution function),
     * of a distribution, for example, "pnorm(x,m,s)". Internally we simply call cdf
     * for the instance of the distribution.
     *
     * \copyright (c) Copyright 2009-2013 (GPL version 3)
     * \author The RevBayes Development Core Team (Sebastian Hoehna)
     *
     */
    class DistributionFunctionCdf : public TypedFunction<Probability> {
        
    public:
        DistributionFunctionCdf(ContinuousDistribution *d); 
        DistributionFunctionCdf(PositiveContinuousDistribution *d);                                                                     //!< Object constructor
        DistributionFunctionCdf(const DistributionFunctionCdf& obj);                                                                    //!< Copy constructor

        virtual                                         ~DistributionFunctionCdf(void);                                                                    //!< Copy constructor

        // overloaded operators
        DistributionFunctionCdf&                        operator=(const DistributionFunctionCdf& c);
        
        // Basic utility functions
        DistributionFunctionCdf*                        clone(void) const;                                                              //!< Clone the object
        static const std::string&                       getClassType(void);                                                             //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                         //!< Get class type spec
        std::vector<std::string>                        getFunctionNameAliases(void) const;                                             //!< Get the aliases of the name of the function in Rev
        std::string                                     getFunctionName(void) const;                                                    //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                                        //!< Get language type of the object
        
        // Regular functions
        RevBayesCore::TypedFunction<double>*            createFunction(void) const;                                                     //!< Create internal function object
        const ArgumentRules&                            getArgumentRules(void) const;                                                   //!< Get argument rules
        
    protected:
        
        ArgumentRules                                   argRules;                                                                       //!< Member rules converted to reference rules
        ContinuousDistribution*                         templateObject;                                                                 //!< The template object
        PositiveContinuousDistribution*                 templateObjectPositive;                                                         //!< The template object
        
    };
    
}


#endif

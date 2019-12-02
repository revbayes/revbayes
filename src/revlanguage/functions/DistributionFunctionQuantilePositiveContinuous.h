/**
 * @file
 * This file contains the declaration of the QuantilePositiveContinuous function of a distribution.
 *
 * @brief Declaration of DistributionFunctionQuantilePositiveContinuous
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-06-01 14:40:37 +0200 (Fri, 01 Jun 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: ConstructorFunction.h 1602 2012-06-01 12:40:37Z hoehna $
 */

#ifndef DistributionFunctionQuantilePositiveContinuous_H
#define DistributionFunctionQuantilePositiveContinuous_H

#include "RlContinuousDistribution.h"
#include "RlTypedFunction.h"
#include "RlPositiveContinuousDistribution.h"

#include <string>
#include <vector>

namespace RevLanguage {
    
    class DistributionFunctionQuantilePositiveContinuous : public TypedFunction<RealPos> {
        
    public:
        DistributionFunctionQuantilePositiveContinuous(PositiveContinuousDistribution *d);                                                                             //!< Object constructor
        DistributionFunctionQuantilePositiveContinuous(const DistributionFunctionQuantilePositiveContinuous& obj);                                                                    //!< Copy constructor
        virtual ~DistributionFunctionQuantilePositiveContinuous(void);
        
        // overloaded operators
        DistributionFunctionQuantilePositiveContinuous&         operator=(const DistributionFunctionQuantilePositiveContinuous& c);
        
        // Basic utility functions
        DistributionFunctionQuantilePositiveContinuous*         clone(void) const;                                                              //!< Clone the object
        static const std::string&                               getClassType(void);                                                             //!< Get Rev type
        static const TypeSpec&                                  getClassTypeSpec(void);                                                         //!< Get class type spec
        std::vector<std::string>                                getFunctionNameAliases(void) const;                                             //!< Get the aliases of the name of the function in Rev
        std::string                                             getFunctionName(void) const;                                                    //!< Get the primary name of the function in Rev
        const TypeSpec&                                         getTypeSpec(void) const;                                                        //!< Get language type of the object
        
        // Regular functions
        RevBayesCore::TypedFunction<double>*                    createFunction(void) const;                                                     //!< Create a function object
        const ArgumentRules&                                    getArgumentRules(void) const;                                                   //!< Get argument rules
        
    protected:
        
        ArgumentRules                                           argRules;                                                                       //!< Member rules converted to reference rules
        PositiveContinuousDistribution*                         templateObjectPositive;                                                               //!< The template object
        
    };
    
}


#endif

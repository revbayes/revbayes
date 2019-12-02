/**
 * @file
 * This file contains the declaration of the QuantileContinuous function of a distribution.
 *
 * @brief Declaration of DistributionFunctionQuantileContinuous
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-06-01 14:40:37 +0200 (Fri, 01 Jun 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: ConstructorFunction.h 1602 2012-06-01 12:40:37Z hoehna $
 */

#ifndef DistributionFunctionQuantileContinuous_H
#define DistributionFunctionQuantileContinuous_H

#include "RlContinuousDistribution.h"
#include "RlTypedFunction.h"

#include <string>
#include <vector>

namespace RevLanguage {
    
    class DistributionFunctionQuantileContinuous : public TypedFunction<Real> {
        
    public:
        DistributionFunctionQuantileContinuous(ContinuousDistribution *d);
        DistributionFunctionQuantileContinuous(const DistributionFunctionQuantileContinuous& obj);                                                                    //!< Copy constructor
        virtual ~DistributionFunctionQuantileContinuous(void);
        
        // overloaded operators
        DistributionFunctionQuantileContinuous&         operator=(const DistributionFunctionQuantileContinuous& c);
        
        // Basic utility functions
        DistributionFunctionQuantileContinuous*         clone(void) const;                                                              //!< Clone the object
        static const std::string&                       getClassType(void);                                                             //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                         //!< Get class type spec
        std::vector<std::string>                        getFunctionNameAliases(void) const;                                             //!< Get the aliases of the name of the function in Rev
        std::string                                     getFunctionName(void) const;                                                    //!< Get the primary name of the function in Rev
        const TypeSpec&                                 getTypeSpec(void) const;                                                        //!< Get language type of the object
        
        // Regular functions
        RevBayesCore::TypedFunction<double>*            createFunction(void) const;                                                     //!< Create a function object
        const ArgumentRules&                            getArgumentRules(void) const;                                                   //!< Get argument rules
        
    protected:
        
        ArgumentRules                                   argRules;                                                                       //!< Member rules converted to reference rules
        ContinuousDistribution*                         templateObject;                                                                 //!< Member rules converted to reference rules
        
    };
    
}


#endif

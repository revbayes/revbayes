/**
 * @file
 * This file contains the declaration of Distribution, which is
 * the interface and abstract base class for RevLanguage distribution objects.
 *
 * @brief Declaration of Distribution
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-08-03 17:51:49 +0200 (Fri, 03 Aug 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @interface Distribution
 * @package distributions
 * @since Version 1.0, 2012-08-06
 *
 * $Id: Function.h 1734 2012-08-03 15:51:49Z hoehna $
 */

#ifndef RlContinuousDistribution_H
#define RlContinuousDistribution_H

#include <cmath>
#include <iosfwd>
#include <string>
#include <vector>

#include "ContinuousDistribution.h"
#include "Real.h"
#include "RlTypedDistribution.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RevPtr.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlStochasticNode.h"
#include "RlTypedFunction.h"
#include "StochasticNode.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

namespace RevLanguage {
class TypeSpec;
    
    class ContinuousDistribution : public TypedDistribution<Real> {
        
    public:
        virtual                                         ~ContinuousDistribution(void);                                                                  //!< Destructor
        ContinuousDistribution(const ContinuousDistribution &x);                                                                //!< Copy constuctor
        
        virtual Real*                                   createRandomVariable(void) const;                                                   //!< Create a random variable from this distribution        
        
        // Basic utility functions you have to override
        virtual ContinuousDistribution*                 clone(void) const = 0;                                                              //!< Clone object
        static const std::string&                       getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                             //!< Get class type spec
        
        
        // Distribution functions you have to override
        virtual RevBayesCore::ContinuousDistribution*   createDistribution(void) const = 0;                                                 //!< Create a random variable from this distribution
        
        
    protected:
        ContinuousDistribution(void);                                                                                                 //!< Basic constructor
        
    };
    
}

#endif


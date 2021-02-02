/**
 * @file
 * This file contains the declaration and implementation of Dist_empiricalTree, which
 * is the Rev language wrapper for the EmpiricalTreeDistribution 
 * class. The distribution is constructed from an input tree trace.
 * 
 *
 * @brief Declaration and implementation of the Dist_empiricalTree class
 *
 * (c) Copyright 2014-
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @since Version 1.0
 *
 */

#ifndef Dist_empiricalTree_H
#define Dist_empiricalTree_H

#include "EmpiricalDistribution.h"
#include "Tree.h"
#include "RlTree.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {
    
	class Dist_empiricalTree : public TypedDistribution<Tree> {
        
    public:
		Dist_empiricalTree( void );                                                                                                                      //!< Constructor
        virtual ~Dist_empiricalTree();                                                                                                                   //!< Virtual destructor
        
        // Basic utility functions
        Dist_empiricalTree*                                      clone(void) const;                                                                      //!< Clone the object
        static const std::string&                                getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                                   getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                              getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        virtual MethodTable                                      getDistributionMethods( void ) const;                                                                       //!< Get the member methods
        const TypeSpec&                                          getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                       getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // overridden function from TypedDistribution
		RevBayesCore::EmpiricalDistribution<RevBayesCore::Tree>* createDistribution(void) const;

    protected:
        
        void                                                     setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                                density;
        RevPtr<const RevVariable>                                trace;
    };
	
}


#endif

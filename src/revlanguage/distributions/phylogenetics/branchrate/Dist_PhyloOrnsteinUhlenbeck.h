//
//  Dist_branchHeterogeneousOrnsteinUhlenbeck.h
//  RevBayes
//
//  Created by Bastien on 05/10/14.
//  Copyright (c) 2014 Bastien Boussau. All rights reserved.
//

#ifndef Dist_PhyloOrnsteinUhlenbeck_H
#define Dist_PhyloOrnsteinUhlenbeck_H


#include <iostream>


#include "PhyloOrnsteinUhlenbeckProcess.h"
#include "ModelVector.h"
#include "RlTypedDistribution.h"
#include "RlSimplex.h"
#include "Real.h"

namespace RevLanguage {
    
    class Dist_PhyloOrnsteinUhlenbeck :  public TypedDistribution< ModelVector<Real> > {
        
    public:
                                                        Dist_PhyloOrnsteinUhlenbeck( void ) {};
        virtual                                        ~Dist_PhyloOrnsteinUhlenbeck() {};
        
        // Basic utility functions
        Dist_PhyloOrnsteinUhlenbeck*                    clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                        getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getMemberRules(void) const;                                                             //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        // Distribution functions you have to override
        RevBayesCore::PhyloOrnsteinUhlenbeckProcess*    createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       tree;
        RevPtr<const RevVariable>                       sigma;
        RevPtr<const RevVariable>                       mean;
        RevPtr<const RevVariable>                       phi;
        RevPtr<const RevVariable>                       root;
        
    };
    
}


#endif	/* DIST_BRANCHHETEROGENEOUSORNSTEINUHLENBECK_H */


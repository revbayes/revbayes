//
//  RlWhiteNoiseProcess.h
//  revbayes
//
//  Created by Nicolas Lartillot on 2014-03-22.
//  Copyright (c) 2014 revbayes team. All rights reserved.
//

#ifndef Dist_PhyloWhiteNoise_H
#define Dist_PhyloWhiteNoise_H

#include <iostream>


#include "PhyloWhiteNoiseProcess.h"
#include "ModelVector.h"
#include "RlTypedDistribution.h"
#include "RlSimplex.h"
#include "RealPos.h"

namespace RevLanguage {

    class Dist_PhyloWhiteNoise :  public TypedDistribution< ModelVector<RealPos>  > {
        
    public:
                                                        Dist_PhyloWhiteNoise( void );
        virtual                                        ~Dist_PhyloWhiteNoise();
        
        // Basic utility functions
        Dist_PhyloWhiteNoise*                           clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        void                                            printValue(std::ostream& o) const;                                                      //!< Print the general information on the function ('usage')
        
        
        
        // Distribution functions you have to override
        RevBayesCore::PhyloWhiteNoiseProcess*           createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       tree;
        RevPtr<const RevVariable>                       sigma;
        
    };
    
}

#endif /* defined(Dist_PhyloWhiteNoise_H) */

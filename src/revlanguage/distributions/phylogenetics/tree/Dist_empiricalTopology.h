#ifndef Dist_empiricalTopology_H
#define Dist_empiricalTopology_H

#include "EmpiricalTopologyDistribution.h"
#include "RlBranchLengthTree.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the empirical topology distribution
     *
     * The RevLanguage wrapper of the empirical topology distribution connects
     * the variables/parameters of the distribution and creates the internal UniformTopologyDistribution object.
     * Please read the UniformTopologyDistribution.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *
     */
    class Dist_empiricalTopology : public TypedDistribution<BranchLengthTree> {
        
    public:
        Dist_empiricalTopology( void );                                                                                                                   //!< Default constructor
        
        // Basic utility functions
        Dist_empiricalTopology*                         clone(void) const;                                                                              //!< Clone the object
        static const std::string&                       getClassType(void);                                                                             //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                         //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                        //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                        //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                                  //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::EmpiricalTopologyDistribution*    createDistribution(void) const;                                                                 //!< Create the internal distribution object
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);               //!< Set member variable
        
    private:
        
        RevPtr<const RevVariable>                       tree_trace;                                                                                     //!< The tree samples
    };
    
}

#endif

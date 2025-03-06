#ifndef Dist_ConstrainedTopology_H
#define Dist_ConstrainedTopology_H

#include "TopologyConstrainedTreeDistribution.h"
#include "RlTypedDistribution.h"
#include "RlTimeTree.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the Dist_ConstrainedTopology
     *
     * The RevLanguage wrapper of the constrained-topology distribution connects
     * the variables/parameters of the process and creates the internal Dist_ConstrainedTopology object.
     * Please read the Dist_ConstrainedTopology.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *
     */
    class Dist_ConstrainedTopology : public TypedDistribution<TimeTree> {
        
    public:
                                                            Dist_ConstrainedTopology(void);
        
        // Basic utility functions
        Dist_ConstrainedTopology*                           clone(void) const;                                                                      //!< Clone the object
        static const std::string&                           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                              getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                            getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                         getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                     getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                  getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        // Distribution functions you have to override
        RevBayesCore::TopologyConstrainedTreeDistribution*  createDistribution(void) const;
        
    protected:
        void                                                setConstParameter(const std::string& name, const RevPtr<const RevVariable>& var);       //!< Set member variable
        
    private:
        RevPtr<const RevVariable>                           base_distribution;
        RevPtr<const RevVariable>                           constraints;
        RevPtr<const RevVariable>                           backbone;
        RevPtr<const RevVariable>                           invert_constraint;                                                                      //!< Boolean indicating whether topologies are rooted
        RevPtr<const RevVariable>                           initial_tree;
        RevPtr<const RevVariable>                           age_check_precision;                                                                    //!< Number of decimal places to use when checking the initial tree against taxon ages
    };
}

#endif

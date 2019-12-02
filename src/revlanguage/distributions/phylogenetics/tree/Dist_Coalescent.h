#ifndef Dist_Coalescent_H
#define Dist_Coalescent_H

#include "ConstantPopulationCoalescent.h"
#include "RlTypedDistribution.h"
#include "RlTimeTree.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the Dist_Coalescent
     *
     * The RevLanguage wrapper of the constant population size coalescent process connects
     * the variables/parameters of the process and creates the internal ConstantPopulationCoalescent object.
     * Please read the ConstantPopulationCoalescent.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-03-05, version 1.0
     *
     */
    class Dist_Coalescent :  public TypedDistribution<TimeTree> {
        
    public:
        Dist_Coalescent( void );
        
        // Basic utility functions
        Dist_Coalescent*                                    clone(void) const;                                                                      //!< Clone the object
        static const std::string&                           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                              getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                         getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                     getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                  getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::ConstantPopulationCoalescent*         createDistribution(void) const;                                                         //!< Create an internal object of the diveristy-dependent pure-birth process.
        
    protected:
        
        void                                                setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        // members
        RevPtr<const RevVariable>                           taxa;                                                                                   //!< The taxon names that will be applied to the initally simulated tree
        RevPtr<const RevVariable>                           constraints;                                                                            //!< Topological constraints that will be used for calibrations
        RevPtr<const RevVariable>                           theta;
        
    };
    
}

#endif

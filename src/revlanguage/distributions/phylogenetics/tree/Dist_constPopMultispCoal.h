#ifndef Dist_constPopMultispCoal_H
#define Dist_constPopMultispCoal_H

#include "MultispeciesCoalescent.h"
#include "RlTimeTree.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the MultispeciesCoalescent process
     *
     * The RevLanguage wrapper of the Multispecies Coalescent with constant population process connects
     * the variables/parameters of the process and creates the internal MultispeciesCoalescent object.
     * Please read the MultispeciesCoalescent.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *
     */
    class Dist_constPopMultispCoal :  public  TypedDistribution<TimeTree> {

    public:
        Dist_constPopMultispCoal( void );

        // Basic utility functions
        Dist_constPopMultispCoal*                           clone(void) const;                                                                      //!< Clone the object
        static const std::string&                           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                              getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                         getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                     getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                  getParameterRules(void) const;                                                          //!< Get member rules (const)


        // Distribution functions you have to override
        RevBayesCore::MultispeciesCoalescent*               createDistribution(void) const;

    protected:

        void                                                setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable


    private:

        RevPtr<const RevVariable>                           Ne;                                                                                     //!< The population size
        RevPtr<const RevVariable>                           species_tree;                                                                            //!< The species tree
        RevPtr<const RevVariable>                           taxa;                                                                                   //!< The taxons


    };

}

#endif

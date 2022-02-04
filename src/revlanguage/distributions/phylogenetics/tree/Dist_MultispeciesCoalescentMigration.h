#ifndef Dist_MultispeciesCoalescentMigration_H
#define Dist_MultispeciesCoalescentMigration_H

#include "MultispeciesCoalescentMigration.h"
#include "RlTimeTree.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the MultispeciesCoalescentMigration process
     *
     * The RevLanguage wrapper of the Multispecies Coalescent with migration process connects
     * the variables/parameters of the process and creates the internal MultispeciesCoalescentMigration object.
     * Please read the MultispeciesCoalescentMigration.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *
     */
    class Dist_MultispeciesCoalescentMigration : public TypedDistribution<TimeTree> {

    public:
        Dist_MultispeciesCoalescentMigration( void );

        // Basic utility functions
        Dist_MultispeciesCoalescentMigration*                   clone(void) const;                                                                      //!< Clone the object
        static const std::string&                               getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                                  getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                             getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        std::vector<std::string>                                getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        const TypeSpec&                                         getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                      getParameterRules(void) const;                                                          //!< Get member rules (const)


        // Distribution functions you have to override
        RevBayesCore::MultispeciesCoalescentMigration*          createDistribution(void) const;

    protected:

        void                                                    setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable


    private:

        RevPtr<const RevVariable>                               Ne;                                                                                     //!< The population size
        RevPtr<const RevVariable>                               species_tree;
        RevPtr<const RevVariable>                               Q;                                                                                     //!< The population size
        RevPtr<const RevVariable>                               delta;                                                                             //!< The species tree
        RevPtr<const RevVariable>                               taxa;                                                                                   //!< The taxons


    };

}

#endif

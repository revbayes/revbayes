#ifndef Dist_MPEST_H
#define Dist_MPEST_H

#include "MPEST.h"
#include "RlTimeTree.h"
#include "RlTypedDistribution.h"
#include "RlRootedTripletDistribution.h"


namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the MPEST process
     *
     * The RevLanguage wrapper of the Multispecies Coalescent with constant population process connects
     * the variables/parameters of the process and creates the internal MPEST object.
     * Please read the MPEST.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *
     */
    class Dist_MPEST :  public  TypedDistribution<RootedTripletDistribution> {
        
    public:
        Dist_MPEST( void );
        
        // Basic utility functions
        Dist_MPEST*                                         clone(void) const;                                                                      //!< Clone the object
        static const std::string&                           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                              getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                         getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                     getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                  getMemberRules(void) const;                                                             //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::MPEST*                                createDistribution(void) const;
        
    protected:
        
        void                                                setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                              geneTrees;                                                                           //!< The population size
        RevPtr<const RevVariable>                              speciesTree;                                                                         //!< The species tree
        RevPtr<const RevVariable>                              useSpecies;                                                                          //!< The taxons
        
        
    };
    
}

#endif

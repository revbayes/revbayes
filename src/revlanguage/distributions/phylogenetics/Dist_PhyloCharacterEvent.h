#ifndef Dist_PhyloCharacterEvent_H
#define Dist_PhyloCharacterEvent_H

#include "PhyloCharacterEventDistribution.h"
#include "RlCharacterHistory.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {
    
    /**
     * @file
     * This file contains the declaration of the Rev character event over a phylogeny distribution
     *
     *
     * @author Sebastian Hoehna
     */
    class Dist_PhyloCharacterEvent : public TypedDistribution<CharacterHistory> {
        
    public:
        Dist_PhyloCharacterEvent( void );                                                                                               //!< Constructor
        virtual ~Dist_PhyloCharacterEvent();                                                                                            //!< Virtual destructor
        
        // Basic utility functions
        Dist_PhyloCharacterEvent*                                       clone(void) const;                                                                      //!< Clone the object
        static const std::string&                                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        virtual MethodTable                                             getDistributionMethods( void ) const;                                                   //!< Get the member methods
        const TypeSpec&                                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::PhyloCharacterEventDistribution*                  createDistribution(void) const;                                                         //!< Create the core object corresponding to this wrapper
        
    protected:
        
        void                                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                                       tree;
        RevPtr<const RevVariable>                                       root_values;
        RevPtr<const RevVariable>                                       base_distribution;
        RevPtr<const RevVariable>                                       event_rate;
        RevPtr<const RevVariable>                                       names;
        
    };
    
}

#endif


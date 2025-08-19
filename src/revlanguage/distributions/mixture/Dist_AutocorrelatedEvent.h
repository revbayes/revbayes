#ifndef Dist_AutocorrelatedEvent_h
#define Dist_AutocorrelatedEvent_h

#include "AutocorrelatedEventDistribution.h"
#include "RlMultiValueEvent.h"
#include "Natural.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the AutocorrelatedEvent distribution.
     *
     * The RevLanguage wrapper of the AutocorrelatedEvent distribution simply
     * manages the interactions through the Rev with our core.
     * That is, the internal distribution object can be constructed and hooked up
     * in a model graph.
     * See the Dist_AutocorrelatedEvent.h for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-18, version 1.0
     *
     */
    class Dist_AutocorrelatedEvent : public TypedDistribution<MultiValueEvent> {
        
    public:
        Dist_AutocorrelatedEvent( void );                                                                                                                           //!< Default constructor
        
        // Basic utility functions
        Dist_AutocorrelatedEvent*                       clone(void) const;                                                                              //!< Clone the object
        static const std::string&                       getClassType(void);                                                                             //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                         //!< Get class type spec
        std::vector<std::string>                        getDistributionFunctionAliases(void) const;                                                     //!< Get the alternative names used for the constructor function in Rev.
        std::string                                     getDistributionFunctionName(void) const;                                                        //!< Get the Rev-name for this distribution.
        MethodTable                                     getDistributionMethods( void ) const;
        const TypeSpec&                                 getTypeSpec(void) const;                                                                        //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                                  //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::AutocorrelatedEventDistribution*  createDistribution(void) const;                                                                 //!< Create the internal distribution object
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);               //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       event_number_prior;
        RevPtr<const RevVariable>                       value_priors;
        RevPtr<const RevVariable>                       names;
        RevPtr<const RevVariable>                       min_elements;
        RevPtr<const RevVariable>                       autocorrelation_types;
        RevPtr<const RevVariable>                       autocorrelation_dep_var;
        RevPtr<const RevVariable>                       autocorrelation_sigmas;
        RevPtr<const RevVariable>                       name_of_var_to_sort_by;

    };
    
}


#endif

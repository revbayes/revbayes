#ifndef Dist_markovTimes_h
#define Dist_markovTimes_h

#include "RlOrderedEventTimes.h"
#include "ModelVector.h"
#include "MarkovTimesDistribution.h"
#include "Natural.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the MultiValueEvent distribution.
     *
     * The RevLanguage wrapper of the MultiValueEvent distribution simply
     * manages the interactions through the Rev with our core.
     * That is, the internal distribution object can be constructed and hooked up
     * in a model graph.
     * See the Dist_markovTimes.h for more details.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-18, version 1.0
     *
     */
    class Dist_markovTimes : public TypedDistribution< RlOrderedEventTimes > {
        
    public:
        Dist_markovTimes( void );                                                                                                                           //!< Default constructor
        
        // Basic utility functions
        Dist_markovTimes*                               clone(void) const;                                                                              //!< Clone the object
        static const std::string&                       getClassType(void);                                                                             //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                         //!< Get class type spec
        std::vector<std::string>                        getDistributionFunctionAliases(void) const;                                                     //!< Get the alternative names used for the constructor function in Rev.
        std::string                                     getDistributionFunctionName(void) const;                                                        //!< Get the Rev-name for this distribution.
        MethodTable                                     getDistributionMethods( void ) const;
        const TypeSpec&                                 getTypeSpec(void) const;                                                                        //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                                  //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::MarkovTimesDistribution*          createDistribution(void) const;                                                                 //!< Create the internal distribution object
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);               //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       rate;
        RevPtr<const RevVariable>                       age;

    };
    
}


#endif

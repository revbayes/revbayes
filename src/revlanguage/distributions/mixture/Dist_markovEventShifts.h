#ifndef Dist_markovEventShifts_H
#define Dist_markovEventShifts_H

#include "OrderedEventTimes.h"
#include "RlOrderedEventTimes.h"
#include "MarkovEventShiftsDistribution.h"
#include "RealPos.h"
#include "RlTypedDistribution.h"
#include "TypeSpec.h"
#include "OrderedEvents.h"
#include "RlOrderedEvents.h"

namespace RevLanguage {
    
    template <typename valType>
    class Dist_markovEventShifts : public TypedDistribution< RlOrderedEvents<valType> >{
        
    public:
                                                        Dist_markovEventShifts( void );
        virtual                                        ~Dist_markovEventShifts();
        
        // Basic utility functions
        Dist_markovEventShifts*                         clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        MethodTable                                     getDistributionMethods( void ) const;
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)

        
        // Distribution functions you have to override
        RevBayesCore::MarkovEventShiftsDistribution<typename valType::valueType>* createDistribution(void) const;
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       event_times;
        RevPtr<const RevVariable>                       initial_event;
        RevPtr<const RevVariable>                       base_distribution;
        
    };
    
}


#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MethodTable.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"


template <typename valType>
RevLanguage::Dist_markovEventShifts<valType>::Dist_markovEventShifts() : TypedDistribution< RlOrderedEvents<valType> >(),
	event_times( NULL ),
	base_distribution( NULL )
{
    
}


template <typename valType>
RevLanguage::Dist_markovEventShifts<valType>::~Dist_markovEventShifts()
{
    
}



template <typename valType>
RevLanguage::Dist_markovEventShifts<valType>* RevLanguage::Dist_markovEventShifts<valType>::clone( void ) const
{
    return new Dist_markovEventShifts(*this);
}


template <typename valType>
RevBayesCore::MarkovEventShiftsDistribution<typename valType::valueType>* RevLanguage::Dist_markovEventShifts<valType>::createDistribution( void ) const
{

    // get the parameters
    const Distribution& rlDistribution									= static_cast<const Distribution &>( base_distribution->getRevObject() );
    RevBayesCore::TypedDistribution<typename valType::valueType>* g0    = static_cast<RevBayesCore::TypedDistribution<typename valType::valueType>* >( rlDistribution.createDistribution() );
    RevBayesCore::TypedDagNode<typename valType::valueType>*      init  = static_cast<const valType &>( initial_event->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::OrderedEventTimes>*  oet   = static_cast<const RlOrderedEventTimes &>( event_times->getRevObject() ).getDagNode();

    RevBayesCore::MarkovEventShiftsDistribution<typename valType::valueType>* d = new RevBayesCore::MarkovEventShiftsDistribution<typename valType::valueType>(oet, init, g0);

    return d;
}



/* Get Rev type of object */
template <typename valType>
const std::string& RevLanguage::Dist_markovEventShifts<valType>::getClassType(void)
{
    
    static std::string rev_type = "Dist_markovEventShifts";
    
	return rev_type;
}

/* Get class type spec describing type of object */
template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_markovEventShifts<valType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< valType >::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
template <typename valType>
std::string RevLanguage::Dist_markovEventShifts<valType>::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "MarkovEventShifts";
    
    return d_name;
}

template <typename valType>
MethodTable RevLanguage::Dist_markovEventShifts<valType>::getDistributionMethods( void ) const
{
	MethodTable methods = TypedDistribution< RlOrderedEvents<valType> >::getDistributionMethods();
    return methods;
}


/** Return member rules (no members) */
template <typename valType>
const RevLanguage::MemberRules& RevLanguage::Dist_markovEventShifts<valType>::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {

		dist_member_rules.push_back( new ArgumentRule( "eventTimes"      , RlOrderedEventTimes::getClassTypeSpec(),        "The times of the events.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    	dist_member_rules.push_back( new ArgumentRule( "initialEvent",     valType::getClassTypeSpec(),                    "The intitial event.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    	dist_member_rules.push_back( new ArgumentRule( "baseDistribution", TypedDistribution<valType>::getClassTypeSpec(), "The base distribution from which rate multipliers are drawn for each event.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

		rules_set = true;
    }
    
    return dist_member_rules;
}


template <typename valType>
const RevLanguage::TypeSpec& RevLanguage::Dist_markovEventShifts<valType>::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}



/** Set a member variable */
template <typename valType>
void RevLanguage::Dist_markovEventShifts<valType>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "eventTimes" )
    {
    	event_times = var;
    }
    else if ( name == "initialEvent" )
    {
    	initial_event = var;
    }
    else if ( name == "baseDistribution" )
    {
        base_distribution = var;
    }
    else
    {
        TypedDistribution< RlOrderedEvents<valType> >::setConstParameter(name, var);
    }
}


#endif

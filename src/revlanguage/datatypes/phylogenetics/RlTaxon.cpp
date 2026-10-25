#include <cstddef>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "ConstantNode.h"
#include "RlTaxon.h"
#include "RlString.h"
#include "RealPos.h"
#include "TypeSpec.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "RbHelpReference.h"
#include "RbException.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "StringUtilities.h"
#include "Taxon.h"
#include "TimeInterval.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevLanguage { class Argument; }

using namespace RevLanguage;

/** Default constructor */
Taxon::Taxon(void) : ModelObject<RevBayesCore::Taxon>()
{
    
    initMethods();

}

/** Construct from core Taxon */
Taxon::Taxon(RevBayesCore::Taxon *c) : ModelObject<RevBayesCore::Taxon>( c )
{
    
    initMethods();

}

/** Construct from core Taxon */
Taxon::Taxon(const RevBayesCore::Taxon &t) : ModelObject<RevBayesCore::Taxon>( new RevBayesCore::Taxon( t ) )
{
    
    initMethods();

}

/** Construct from DAG node */
Taxon::Taxon(RevBayesCore::TypedDagNode<RevBayesCore::Taxon> *n) : ModelObject<RevBayesCore::Taxon>( n )
{
    
    initMethods();

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Taxon* Taxon::clone(void) const
{
    
    return new Taxon(*this);
}


void Taxon::constructInternalObject( void ) 
{
    // Validate the age range before replacing the existing taxon.
    std::string taxonName = static_cast<const RlString &>( (taxon)->getRevObject() ).getValue() ;
    std::string taxonSpecies = static_cast<const RlString &>( (species)->getRevObject() ).getValue() ;
    taxonSpecies = ( taxonSpecies == "taxonName" ) ? taxonName : taxonSpecies;
    // A nonnumeric default distinguishes omitted ages from an explicitly supplied zero.
    const bool has_min = minage->getRevObject().isType(RealPos::getClassTypeSpec());
    const bool has_max = maxage->getRevObject().isType(RealPos::getClassTypeSpec());
    // Supplying only a maximum explicitly establishes an interval with lower endpoint zero.
    double minAge = has_min ? static_cast<const RealPos&>(minage->getRevObject()).getValue() : 0.0;
    double maxAge = minAge;

    if( maxage->getRevObject().isType( RealPos::getClassTypeSpec() ) )
    {
        maxAge = static_cast<const RealPos &>( maxage->getRevObject() ).getValue();
    }

    if (std::isnan(minAge) or std::isnan(maxAge))
        throw RbException() << "Age bounds for taxon '" << taxonName << "' must not be NaN.";
    if (minAge > maxAge)
        throw RbException() << "Invalid age range for taxon '" << taxonName << "': minimum age ("
                            << minAge << ") exceeds maximum age (" << maxAge << ").";
    RevBayesCore::TimeInterval range(minAge, maxAge);

    // Release the previous node only after the new range has passed validation.
    if (dag_node != NULL and dag_node->decrementReferenceCount() == 0)
        delete dag_node;
    
    dag_node = new RevBayesCore::ConstantNode<RevBayesCore::Taxon>("", new RevBayesCore::Taxon( taxonName ) );
    
    dag_node->getValue().setSpeciesName( taxonSpecies );
    
    if (has_min or has_max)
        dag_node->getValue().setAgeRange(range);

    dag_node->incrementReferenceCount();
    
}


/* Map calls to member methods */
RevLanguage::RevPtr<RevLanguage::RevVariable> Taxon::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "getName")
    {
        found = true;
        
        std::string n = this->dag_node->getValue().getName();
        return RevPtr<RevVariable>( new RevVariable( new RlString( n ) ) );
    }
    else if (name == "getSpeciesName")
    {
        found = true;
        
        std::string n = this->dag_node->getValue().getSpeciesName();
        return RevPtr<RevVariable>( new RevVariable( new RlString( n ) ) );
    }
    if (name == "getAge")
    {
        found = true;
        
        // Compatibility: retain Rev's getAge() as a lower-bound alias pending a usage review.
        // COMPAT(age-range): Preserve scalar getters and endpoint validation for unset ranges.
        const auto interval = this->dag_node->getValue().getAgeRange().value_or(RevBayesCore::TimeInterval(0.0, 0.0));
        double a = interval.getMin();
        return RevPtr<RevVariable>( new RevVariable( new RealPos( a ) ) );
    }
    if (name == "getMinAge")
    {
        found = true;

        // COMPAT(age-range): Preserve scalar getters and endpoint validation for unset ranges.
        const auto interval = this->dag_node->getValue().getAgeRange().value_or(RevBayesCore::TimeInterval(0.0, 0.0));
        double a = interval.getMin();
        return RevPtr<RevVariable>( new RevVariable( new RealPos( a ) ) );
    }
    if (name == "getMaxAge")
    {
        found = true;

        // COMPAT(age-range): Preserve scalar getters and endpoint validation for unset ranges.
        const auto interval = this->dag_node->getValue().getAgeRange().value_or(RevBayesCore::TimeInterval(0.0, 0.0));
        double a = interval.getMax();
        return RevPtr<RevVariable>( new RevVariable( new RealPos( a ) ) );
    }
    if (name == "setMinAge")
    {
        found = true;
        double a = static_cast<const RealPos&>( args[0].getVariable()->getRevObject() ).getValue();

        if (std::isnan(a))
            throw RbException() << "Minimum age for taxon '" << dag_node->getValue().getName()
                                << "' must not be NaN.";
        // COMPAT(age-range): Preserve scalar getters and endpoint validation for unset ranges.
        const auto interval = dag_node->getValue().getAgeRange().value_or(RevBayesCore::TimeInterval(0.0, 0.0));
        if (a > interval.getMax())
            throw RbException() << "Cannot set minimum age of taxon '" << dag_node->getValue().getName()
                                << "' to " << a << ": it exceeds the current maximum age ("
                                << interval.getMax() << ").";

        dag_node->getValue().setAgeRange(RevBayesCore::TimeInterval(a, interval.getMax()));
        return NULL;
    }
    if (name == "setMaxAge")
    {
        found = true;
        double a = static_cast<const RealPos&>( args[0].getVariable()->getRevObject() ).getValue();

        if (std::isnan(a))
            throw RbException() << "Maximum age for taxon '" << dag_node->getValue().getName()
                                << "' must not be NaN.";
        // COMPAT(age-range): Preserve scalar getters and endpoint validation for unset ranges.
        const auto interval = dag_node->getValue().getAgeRange().value_or(RevBayesCore::TimeInterval(0.0, 0.0));
        if (a < interval.getMin())
            throw RbException() << "Cannot set maximum age of taxon '" << dag_node->getValue().getName()
                                << "' to " << a << ": it is below the current minimum age ("
                                << interval.getMin() << ").";

        dag_node->getValue().setAgeRange(RevBayesCore::TimeInterval(interval.getMin(), a));
        return NULL;
    }
    
    if (name == "setAgeRange")
    {
        found = true;
        double minAge = static_cast<const RealPos&>(args[0].getVariable()->getRevObject()).getValue();
        double maxAge = static_cast<const RealPos&>(args[1].getVariable()->getRevObject()).getValue();

        if (std::isnan(minAge) or std::isnan(maxAge))
            throw RbException() << "Age bounds for taxon '" << dag_node->getValue().getName()
                                << "' must not be NaN.";
        if (minAge > maxAge)
            throw RbException() << "Invalid age range for taxon '" << dag_node->getValue().getName()
                                << "': minimum age (" << minAge << ") exceeds maximum age (" << maxAge << ").";

        // Validate both bounds before replacing the range, so failed updates leave it unchanged.
        dag_node->getValue().setAgeRange(RevBayesCore::TimeInterval(minAge, maxAge));
        return NULL;
    }

    return ModelObject<RevBayesCore::Taxon>::executeMethod( name, args, found );
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Taxon::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "taxon";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Taxon::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new ArgumentRule("taxonName"  , RlString::getClassTypeSpec(), "The name of the taxon.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("speciesName", RlString::getClassTypeSpec(), "The name of the species it belongs to.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("taxonName") ) );
        std::vector<std::string> ageLabels;
        ageLabels.push_back("age");
        ageLabels.push_back("min");
        memberRules.push_back( new ArgumentRule(ageLabels,     RealPos::getClassTypeSpec(), "The (minimum) age before the present when this taxon was sampled.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("unspecified") ) );
        memberRules.push_back( new ArgumentRule("max",         RealPos::getClassTypeSpec(), "The maximum age before the present when this taxon was sampled.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("age/min") ) );

        rules_set = true;
    }
    
    return memberRules;
}


/** Get Rev type of object */
const std::string& Taxon::getClassType(void)
{
    
    static std::string rev_type = "Taxon";
    
    return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Taxon::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
    return rev_type_spec; 
}


/** Get type spec */
const TypeSpec& Taxon::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/**
 * Initialize the member methods.
 */
void Taxon::initMethods( void )
{
    
    ArgumentRules* name_arg_rules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "getName", RlString::getClassTypeSpec(), name_arg_rules ) );
    
    ArgumentRules* species_name_arg_rules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "getSpeciesName", RlString::getClassTypeSpec(), species_name_arg_rules ) );
    
    ArgumentRules* age_arg_rules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "getAge", RlString::getClassTypeSpec(), age_arg_rules ) );
    
    age_arg_rules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "getMinAge", RlString::getClassTypeSpec(), age_arg_rules ) );

    age_arg_rules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "getMaxAge", RlString::getClassTypeSpec(), age_arg_rules ) );

    ArgumentRules* set_min_age_arg_rules = new ArgumentRules();
    set_min_age_arg_rules->push_back( new ArgumentRule( "a", RealPos::getClassTypeSpec(), "The min age of the clade.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "setMinAge", RlUtils::Void, set_min_age_arg_rules ) );

    ArgumentRules* set_max_age_arg_rules = new ArgumentRules();
    set_max_age_arg_rules->push_back( new ArgumentRule( "a", RealPos::getClassTypeSpec(), "The max age of the clade.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "setMaxAge", RlUtils::Void, set_max_age_arg_rules ) );

    ArgumentRules* set_age_range_arg_rules = new ArgumentRules();
    set_age_range_arg_rules->push_back(new ArgumentRule("minAge", RealPos::getClassTypeSpec(),
        "The minimum age of the taxon.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY));
    set_age_range_arg_rules->push_back(new ArgumentRule("maxAge", RealPos::getClassTypeSpec(),
        "The maximum age of the taxon.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY));
    methods.addFunction(new MemberProcedure("setAgeRange", RlUtils::Void, set_age_range_arg_rules));
}


/** Set a member variable */
void Taxon::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "taxonName") 
    {
        taxon = var ;
    } 
    else if ( name == "speciesName") 
    {
        species = var ;
    } 
    else if ( name == "age" || name == "min" || name == "age/min" )
    {
        minage = var ;
    }
    else if ( name == "max")
    {
        maxage = var ;
    } 
    else {
        RevObject::setConstParameter(name, var);
    }
    
}



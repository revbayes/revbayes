#include <cstddef>
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
    // we free the memory first
    if ( dag_node != NULL )
    {
        if ( dag_node->decrementReferenceCount() == 0 )
        {
            delete dag_node;
        }
    }
    
    // now allocate a new Taxon
    std::string taxonName = static_cast<const RlString &>( (taxon)->getRevObject() ).getValue() ;
    std::string taxonSpecies = static_cast<const RlString &>( (species)->getRevObject() ).getValue() ;
    taxonSpecies = ( taxonSpecies == "taxonName" ) ? taxonName : taxonSpecies;
    double minAge = static_cast<const RealPos &>( minage->getRevObject() ).getValue();
    double maxAge = minAge;

    if( maxage->getRevObject().isType( RealPos::getClassTypeSpec() ) )
    {
        maxAge = static_cast<const RealPos &>( maxage->getRevObject() ).getValue();
    }
    
    dag_node = new RevBayesCore::ConstantNode<RevBayesCore::Taxon>("", new RevBayesCore::Taxon( taxonName ) );
    
    dag_node->getValue().setSpeciesName( taxonSpecies );
    
    RevBayesCore::TimeInterval range(minAge, maxAge);
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
        
        double a = this->dag_node->getValue().getAge();
        return RevPtr<RevVariable>( new RevVariable( new RealPos( a ) ) );
    }
    if (name == "getMinAge")
    {
        found = true;

        double a = this->dag_node->getValue().getAgeRange().getMin();
        return RevPtr<RevVariable>( new RevVariable( new RealPos( a ) ) );
    }
    if (name == "getMaxAge")
    {
        found = true;

        double a = this->dag_node->getValue().getAgeRange().getMax();
        return RevPtr<RevVariable>( new RevVariable( new RealPos( a ) ) );
    }
    if (name == "setMinAge")
    {
        found = true;
        double a = static_cast<const RealPos&>( args[0].getVariable()->getRevObject() ).getValue();

        dag_node->getValue().setMinAge( a );
        return NULL;
    }
    if (name == "setMaxAge")
    {
        found = true;
        double a = static_cast<const RealPos&>( args[0].getVariable()->getRevObject() ).getValue();

        dag_node->getValue().setMaxAge( a );
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
        memberRules.push_back( new ArgumentRule(ageLabels,     RealPos::getClassTypeSpec(), "The (minimum) age before the present when this taxon was sampled.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(0.0) ) );
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



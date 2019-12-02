#include "ConstantNode.h"
#include "ModelVector.h"
#include "RlTaxon.h"
#include "RbUtil.h"
#include "RlString.h"
#include "RealPos.h"
#include "TypeSpec.h"

#include <sstream>

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


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Taxon::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Michael Landis" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Taxon::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string descriptions = "The taxon function creates a Taxon object.";
    
    return descriptions;
}


/**
 * Get the more detailed description of the function
 */
std::string Taxon::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "Each Taxon object records that taxon's name in addition to other information, such as age (which is non-zero for fossils). Character matrices and trees contain Taxon vectors (Taxon[]) that are used to match leaf nodes to data entries for phylogenetic analyses. For multispecies coalescent analyses, Taxon objects are also used to assign species memberships to individuals.";

    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Taxon::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# we can create a Taxon object\n";
    example += "taxon_gorilla = taxon(\"Gorilla_gorilla\")\n";
    example += "# we can create a dummy vector of Taxon objects for simulation\n";
    example += "for (i in 1:10) { taxa[i] = taxon(\"Taxon\"+i) }\n";
    example += "phy ~ dnBDP(lambda=1, mu=0, rootAge=1, taxa=taxa)\n";
    example += "# retrieve the taxon list for 'phy'\n";
    example += "phy.taxa()\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Taxon::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Taxon::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "readTaxonData" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Taxon::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Taxon object";
    
    return title;
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



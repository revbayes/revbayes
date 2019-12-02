#include "ArgumentRule.h"
#include "CoalaFunction.h"
#include "Func_coala.h"
#include "MatrixReal.h"
#include "ModelVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlMatrixReal.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"


using namespace RevLanguage;


/** Default constructor */
Func_coala::Func_coala( void ) : TypedFunction< ModelVector<Real> >()
{

}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_coala* Func_coala::clone( void ) const
{
    return new Func_coala( *this );
}


/** Execute function: Construct coala from RealPos values. */
RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >* Func_coala::createFunction( void ) const
{
    const RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> > *coordinates = static_cast< ModelVector<Real> & >( args[0].getVariable()->getRevObject() ).getDagNode();
    const RevBayesCore::MatrixReal &coa = static_cast< MatrixReal& >( args[1].getVariable()->getRevObject() ).getValue();
    const std::vector<double> &w = static_cast< ModelVector<RealPos>& >( args[2].getVariable()->getRevObject() ).getValue();


    RevBayesCore::CoalaFunction* func = new RevBayesCore::CoalaFunction( coordinates, coa, w );

    return func;
}


/** Get argument rules */
const ArgumentRules& Func_coala::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "coordinates",   ModelVector<Real>::getClassTypeSpec(),    "A vector of coordinates.",  ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "corAnalysis",   MatrixReal::getClassTypeSpec(),           "A correspondence analysis object.",  ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "weights",       ModelVector<RealPos>::getClassTypeSpec(), "A vector of weight for the coordinates.",  ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return argumentRules;
}


/** Get Rev type of object (static) */
const std::string& Func_coala::getClassType( void )
{
    static std::string rev_type = "Func_coala";

    return rev_type;
}


/** Get Rev type spec of object (static) */
const TypeSpec& Func_coala::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &Function::getClassTypeSpec() );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_coala::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnCoala";

    return f_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_coala::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Bastien Boussau" );

    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_coala::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description;

    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Func_coala::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details;

    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Func_coala::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";

    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Func_coala::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    RevBayesCore::RbHelpReference ref = RevBayesCore::RbHelpReference();
    ref.setCitation("A branch-heterogeneous model of protein evolution for efficient inference of ancestral sequences. Groussin M, Boussau B, Gouy M. Syst Biol. 2013 Jul;62(4):523-38.");
    ref.setDoi("10.1093/sysbio/syt016");
    ref.setUrl("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3676677/");

    references.push_back(ref);


    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_coala::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;


    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Func_coala::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "";

    return title;
}


/** Get Rev type spec of object from an instance */
const TypeSpec& Func_coala::getTypeSpec( void ) const
{
    return getClassTypeSpec();
}

#include "RlMultiValueEvent.h"
#include "RbUtil.h"
#include "TypeSpec.h"

#include <sstream>

using namespace RevLanguage;

/** Default constructor */
MultiValueEvent::MultiValueEvent(void) : ModelObject<RevBayesCore::MultiValueEvent>()
{
    
    initMethods();
    
}

/** Construct from core MultiValueEvent */
MultiValueEvent::MultiValueEvent(RevBayesCore::MultiValueEvent *c) : ModelObject<RevBayesCore::MultiValueEvent>( c )
{
    
    initMethods();
    
}

/** Construct from core MultiValueEvent */
MultiValueEvent::MultiValueEvent(const RevBayesCore::MultiValueEvent &t) : ModelObject<RevBayesCore::MultiValueEvent>( new RevBayesCore::MultiValueEvent( t ) )
{
    
    initMethods();
    
}

/** Construct from DAG node */
MultiValueEvent::MultiValueEvent(RevBayesCore::TypedDagNode<RevBayesCore::MultiValueEvent> *n) : ModelObject<RevBayesCore::MultiValueEvent>( n )
{
    
    initMethods();
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
MultiValueEvent* MultiValueEvent::clone(void) const
{
    
    return new MultiValueEvent(*this);
}


/* Map calls to member methods */
RevLanguage::RevPtr<RevLanguage::RevVariable> MultiValueEvent::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
//    if (name == "getName")
//    {
//        found = true;
//
//        std::string n = this->dag_node->getValue().getName();
//        return RevPtr<RevVariable>( new RevVariable( new RlString( n ) ) );
//    }
    
    return ModelObject<RevBayesCore::MultiValueEvent>::executeMethod( name, args, found );
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> MultiValueEvent::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string MultiValueEvent::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string descriptions = "";
    
    return descriptions;
}


/**
 * Get the more detailed description of the function
 */
std::string MultiValueEvent::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "";
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string MultiValueEvent::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# \n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> MultiValueEvent::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> MultiValueEvent::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
//    see_also.push_back( "readMultiValueEventData" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string MultiValueEvent::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "MultiValueEvent object";
    
    return title;
}


/** Get Rev type of object */
const std::string& MultiValueEvent::getClassType(void)
{
    
    static std::string rev_type = "MultiValueEvent";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& MultiValueEvent::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/** Get type spec */
const TypeSpec& MultiValueEvent::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/**
 * Initialize the member methods.
 */
void MultiValueEvent::initMethods( void )
{
    
//    ArgumentRules* name_arg_rules = new ArgumentRules();
//    methods.addFunction( new MemberProcedure( "getName", RlString::getClassTypeSpec(), name_arg_rules ) );
    
}



#include "RlMultiValueEvent.h"

#include "TypeSpec.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevVariable.h"

namespace RevLanguage { class Argument; }

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



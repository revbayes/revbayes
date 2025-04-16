#include <cstddef>
#include <sstream>
#include <string>

#include "ArgumentRule.h"
#include "RbException.h"
#include "SyntaxFormal.h"
#include "Workspace.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "SyntaxElement.h"

namespace RevLanguage { class Environment; }
namespace RevLanguage { class TypeSpec; }

using namespace RevLanguage;


/**
 * Constructor from formal argument specification with implicit type. Note that we
 * make the argument type dynamic by default (false to ArgumentRule constructor).
 * This is required if the function is to work in a dynamic expression as expected.
 * Of course, one can argue whether implicit-type formals are appropriate in such
 * functions, but it is at least conceivable that a user might use implicit-type
 * formals in a function they intend to use in a dynamic expression. If they
 * specify the return type of the function, there is also a reasonable chance that
 * their function will not generate a type error in a dynamic expression.
 */
SyntaxFormal::SyntaxFormal( const std::string& label, SyntaxElement* defaultVal ) :
    SyntaxElement(),
    argRule( NULL )
{
    // Find the type specification
    const TypeSpec& typeSpec = RevObject::getClassTypeSpec();
    
    // Make argument rule from the formal specification
    if ( defaultVal == NULL )
    {
        argRule = new ArgumentRule( label, typeSpec, "", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY );
    }
    else
    {
        RevObject* defaultObj = defaultVal->evaluateContent( Workspace::userWorkspacePtr() )->getRevObject().clone();
        argRule = new ArgumentRule( label, typeSpec, "", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, defaultObj );
    }
}


/**
 * Constructor from formal argument specification with explicit type. The only
 * difference from the constructor with implicit type is that we pay attention to
 * the type specification provided as part of the formal argument specification.
 *
 * Note that the type specification of formals may include a specifier prefix, such as
 * 'const', 'dynamic', 'mutable' or 'stochastic'. In that case, the type specification
 * will consist of two words, the first of which is the modifier, and the second of
 * which is the type.
 */
SyntaxFormal::SyntaxFormal( const std::string& type, const std::string& label, SyntaxElement* defaultVal ) :
    SyntaxElement(),
    argRule( NULL )
{
    // Divide up the type specification into modifier and type
    std::string modifier;
    std::string typeName;
    size_t breakpos = type.find_first_of( " " );
    if ( breakpos == std::string::npos )
    {
        modifier = "";
        typeName = type;
    }
    else
    {
        modifier = type.substr( 0, breakpos );
        typeName = type.substr( breakpos + 1, type.length() - breakpos );
    }
        
    // Check that we have a supported modifier
    if ( modifier != "" && modifier != "const" && modifier != "dynamic" )
        throw RbException() << "Formal type modifier '" << modifier << "' not supported (yet)" ; 
    
    // Find the type specification
    const TypeSpec& typeSpec = Workspace::userWorkspace().getClassTypeSpecOfType( typeName );

    // Generate the default object
    RevObject* defaultObj = NULL;
    if ( defaultVal != NULL )
        defaultObj = defaultVal->evaluateContent( Workspace::userWorkspacePtr() )->getRevObject().clone();
    
    // Now generate argument rule
    if ( modifier == "const" )
    {
        // Make const argument rule from element
        if ( defaultVal == NULL )
        {
            argRule = new ArgumentRule( label, typeSpec, "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY );
        }
        else
        {
            argRule = new ArgumentRule( label, typeSpec, "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, defaultObj );
        }
        
    }
    else // if ( modifier == "dynamic" || modifier == "" )
    {
        // Make dynamic argument rule from element
        if ( defaultVal == NULL )
        {
            argRule = new ArgumentRule( label, typeSpec, "", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY );
        }
        else
        {
            argRule = new ArgumentRule( label, typeSpec, "", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, defaultObj );
        }
        
    }
    
}


/** Deep copy constructor */
SyntaxFormal::SyntaxFormal( const SyntaxFormal& x ) :
    SyntaxElement(x) ,
    argRule( x.argRule->clone() )
{
}


/** Destructor deletes the argument rule, which we own */
SyntaxFormal::~SyntaxFormal()
{
    delete argRule;
}


/** Assignment operator */
SyntaxFormal& SyntaxFormal::operator=(const SyntaxFormal& x) {

    if (&x != this)
    {
        SyntaxElement::operator=(x);

        delete argRule;
        if ( x.argRule != NULL )
            argRule = x.argRule->clone();
        else
            argRule = NULL;
    }

    return (*this);
}


/** Type-safe clone of syntax element */
SyntaxFormal* SyntaxFormal::clone () const
{
    return new SyntaxFormal(*this);
}


/** Get the processed argument rule (non-const version) */
ArgumentRule* SyntaxFormal::getArgumentRule(void )
{
    return argRule;
}


/** Get the processed argument rule (const version) */
const ArgumentRule* SyntaxFormal::getArgumentRule(void ) const
{
    return argRule;
}


/** Get semantic value (not applicable so return NULL) */
RevPtr<RevVariable> SyntaxFormal::evaluateContent( const std::shared_ptr<Environment>& env, bool dynamic )
{
    return NULL;
}


/**
 * Set the protected flag of the formal argument specification.
 *
 * @todo Support the protected modifier of member variable slots.
 */
void SyntaxFormal::setIsProtected( bool prot )
{
//    argRule->setIsProtected( prot );
}


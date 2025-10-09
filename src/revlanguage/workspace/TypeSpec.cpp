#include <cstddef>
#include <ostream>
#include <string>

#include "TypeSpec.h"

using namespace RevLanguage;


/**
 * Constructor from type and parent type spec. This constructor
 * is appropriate for all classes except those that derive
 * immediately from ModelContainer or WorkspaceContainer.
 */
TypeSpec::TypeSpec( const std::string& objType, const TypeSpec* par ) :
    parent( par ),
    type( objType )
{
    
    element = getElementTypeSpecFromClass();
}


/**
 * Constructor from type, parent type spec, and element type spec.
 * This constructor is used by types that derive directly from ModelContainer
 * or WorkspaceContainer.
 */
TypeSpec::TypeSpec( const std::string& objType, const TypeSpec* par, const TypeSpec* elemType ) :
    element( elemType ),
    parent( par ),
    type( objType )
{
    
}


/**
 * Equals comparison. This assumes that the Rev type
 * names are unique, which should be safe. Rev does
 * not support overloaded types.
 */
bool TypeSpec::operator==( const TypeSpec& x ) const
{
    return type == x.type;
}


/**
 * Not equals comparison. Simply rely on negating
 * the equals comparison, as usual.
 */
bool TypeSpec::operator!=( const TypeSpec& x ) const
{
    return type != x.type;
}


/**
 * Get element type. Simply use the element type spec
 * that we computed in the constructor.
 */
const std::string& TypeSpec::getElementType( void ) const
{
    static std::string emptyString;
    
    if ( element == NULL)
        return emptyString;
    else
        return element->getType();
}


/**
 * Get element type spec. Simply return the element type spec
 * that we computed in the constructor.
 */
const TypeSpec* TypeSpec::getElementTypeSpec( void ) const
{
    return element;
}


/**
 * Get the element type by climbing the ladder to see whether we are
 * derived from ModelContainer or WorkspaceContainer. Then simply skip
 * the trailing square brackets to get the element type.
 */
const TypeSpec* TypeSpec::getElementTypeSpecFromClass(void) const
{
    // If we reached the end there is no element
    if ( parent == NULL )
        return NULL;
    
    // Get parent type
    std::string parentType;
    if ( parent != NULL )
        parentType = parent->getType();
    
    // Check if we now have the answer
    if ( parentType == "ModelContainer" || parentType == "WorkspaceContainer" )
        return this->element;
    
    // If not, continue climbing up our derivation tree
    return parent->getElementTypeSpecFromClass();
}


/**
 * Get parent type string. Simply return the type string of the
 * parent type spec but throw an error if the parent is NULL.
 */
const std::string& TypeSpec::getParentType( void ) const
{
    static std::string emptyString;
    
    if ( parent == NULL)
    {
        return emptyString;
    }
    else
    {
        return parent->getType();
    }
    
}


/**
 * Get parent type spec. Simply return the parent type spec
 * but throw an error if the parent is NULL.
 */
const TypeSpec* TypeSpec::getParentTypeSpec( void ) const
{
    return parent;
}


/**
 * Get the type string. Simply return the
 * type member variable.
 */
const std::string& TypeSpec::getType(void) const
{
    return type;
}


/**
 * Test whether we are the same type or derived from the same
 * type as the argument type spec. We do this by a preorder
 * traversal of the tree.
 */
bool TypeSpec::isDerivedOf( const TypeSpec& x ) const
{
    // First we test whether this type is the same as the argument type
    if ( operator==( x ) )
        return true;
    
    // We reached the end
    if ( parent == NULL )
        return false;

    // If not, continue climbing up our derivation tree
    return parent->isDerivedOf(x);
}


/**
 * Output the class vector (type spec) to a stream.
 */
std::ostream& RevLanguage::operator<<(std::ostream& o, const TypeSpec& x)
{
    o << x.getType();
    
    const TypeSpec* par = x.getParentTypeSpec();
    while ( par != NULL ) {
        
        o << ", " << par->getType();
        par = par->getParentTypeSpec();
    }

    return o;
}




#include "RevAbstractType.h"

#include "TypeSpec.h"
#include "RevObject.h"

using namespace RevLanguage;


/** Constructor from type specification */
RevAbstractType::RevAbstractType( const TypeSpec& t, RevObject* obj ) : WorkspaceObject(),
    myTypeSpec( t ),
    exampleObject( obj )
{
}


/** Copy constructor (really not needed) */
RevAbstractType::RevAbstractType(const RevAbstractType& x) : WorkspaceObject( x ),
    myTypeSpec( x.myTypeSpec )
{
}


/** Generate type-safe clone */
RevAbstractType* RevAbstractType::clone(void) const
{
    return new RevAbstractType(*this);
}


/** Get type spec */
const TypeSpec& RevAbstractType::getTypeSpec( void ) const
{
    
    return myTypeSpec;
}


/**
 * We represent an abstract type, so return true.
 */
bool RevAbstractType::isAbstract( void ) const
{
    return true;
}


/** Make a new instance of the example object of a non-abstract derived class */
RevObject* RevAbstractType::makeExampleObject( void ) const
{
    return exampleObject->clone();
}


/** Print the value (report an abstract type and then the type) */
void RevAbstractType::printValue(std::ostream &o, bool user) const
{
    
    o << "<abstract object of type '" << myTypeSpec.getType() << "'";
}



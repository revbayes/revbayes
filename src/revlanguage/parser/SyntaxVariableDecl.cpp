#include <iostream>
#include <list>
#include <cstddef>
#include <set>
#include <string>
#include <vector>

#include "Environment.h"
#include "RbException.h"
#include "SyntaxVariableDecl.h"
#include "Workspace.h"
#include "Natural.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "SyntaxElement.h"

using namespace RevLanguage;

/** Basic constructor. Construct from type name, lengths specification, and variable identifier. */
SyntaxVariableDecl::SyntaxVariableDecl( const std::string& elemType, std::list<SyntaxElement*>* lengths, const std::string& varName) :
    SyntaxElement(),
    elementTypeName( elemType ),
    lengthExpr( lengths ),
    variableName( varName )
{
}


/** Deep copy constructor */
SyntaxVariableDecl::SyntaxVariableDecl( const SyntaxVariableDecl& x ) :
    SyntaxElement( x ),
    elementTypeName( x.elementTypeName ),
    lengthExpr( NULL ),
    variableName( x.variableName )
{
    lengthExpr = new std::list<SyntaxElement*>();
    
    for ( std::list<SyntaxElement*>::const_iterator it = x.lengthExpr->begin(); it != x.lengthExpr->end(); ++it )
        lengthExpr->push_back( (*it)->clone() );
}


/** Destructor deletes operands */
SyntaxVariableDecl::~SyntaxVariableDecl( void )
{
    for (std::list<SyntaxElement*>::iterator it = lengthExpr->begin(); it != lengthExpr->end(); ++it )
        delete (*it);

    delete lengthExpr;
}


/** Assignment operator */
SyntaxVariableDecl& SyntaxVariableDecl::operator=(const SyntaxVariableDecl& x)
{
    if ( this != &x ) {
        
        SyntaxElement::operator=( x );
        
        for (std::list<SyntaxElement*>::iterator it = lengthExpr->begin(); it != lengthExpr->end(); ++it )
            delete (*it);
        
        lengthExpr->clear();

        elementTypeName = x.elementTypeName;
        variableName    = x.variableName;
        
        for ( std::list<SyntaxElement*>::const_iterator it = x.lengthExpr->begin(); it != x.lengthExpr->end(); ++it )
            lengthExpr->push_back( (*it)->clone() );
    }
    
    return (*this);
}



/** Type-safe clone of syntax element */
SyntaxVariableDecl* SyntaxVariableDecl::clone() const
{
    return new SyntaxVariableDecl( *this );
}


/**
 * @brief Get semantic value
 * 
 * Here we evaluate the length specification (statically) and create the
 * requested variable.
 *
 */
RevPtr<RevVariable> SyntaxVariableDecl::evaluateContent( Environment& env, bool dynamic )
{
    
    // Check if variable exists
    if ( env.existsVariable( variableName ) )
        throw RbException() << "Illegal attempt to redefine variable " << variableName ; 
    
    // Check if type exists
    if ( !Workspace::userWorkspace().existsType( elementTypeName ) )
        throw RbException() << "Type '" << elementTypeName << "' does not exist" ; 

    // Evaluate length specification
    std::vector<size_t> lengths;
    for ( std::list<SyntaxElement*>::iterator it = lengthExpr->begin(); it != lengthExpr->end(); ++it )
    {
        if ( (*it) == NULL )
        {
            lengths.push_back( 1 );
        }
        else
        {
            RevPtr<RevVariable> temp    = (*it)->evaluateContent( env, dynamic );
            const RevObject& value   = temp->getRevObject();
            
            size_t theLength;
            if ( value.isType( Natural::getClassTypeSpec() ) )
                theLength = size_t( static_cast<const Natural&>( value ).getValue() );
            else if ( value.isConvertibleTo( Natural::getClassTypeSpec(), true ) )
            {
                RevObject* convObj = value.convertTo( Natural::getClassTypeSpec() );
                theLength = size_t( static_cast<Natural*>( convObj )->getValue() );
                delete convObj;
            }
            else
                throw RbException( "Length specification does not evaluate to an object of type 'Natural'" );
            
            if ( theLength == 0 )
                throw RbException( "Invalid length specification (0)" );
            
            lengths.push_back( theLength );
        }
    }
    
    // We ask the user workspace for the new objects
    RevObject* newObject;
    if ( lengths.size() == 0 )
        newObject = Workspace::userWorkspace().makeNewDefaultObject( elementTypeName );
    else
    {
        throw RbException("This needs replacements!!!");
//        newObject = Workspace::userWorkspace().makeNewEmptyContainer( elementTypeName, lengths.size() );
//        static_cast<Container*>( newObject )->resize( lengths );
    }

    // Add the new RevVariable
    env.addVariable( variableName, new RevVariable( newObject ) );
    
    return NULL;
}


/**
 * Is the syntax element safe for use in a function (as
 * opposed to a procedure)? The element is safe but we
 * need to add the variable to the local variables
 */
bool SyntaxVariableDecl::isFunctionSafe( const Environment& env, std::set<std::string>& localVars ) const
{
    localVars.insert( variableName );
    return true;
}



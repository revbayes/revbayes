#include "RlCharacterHistory.h"

#include <stddef.h>
#include <iostream>

#include "ArgumentRule.h"
#include "TypeSpec.h"

using namespace RevLanguage;

CharacterHistory::CharacterHistory(void) : ModelObject<RevBayesCore::CharacterHistoryDiscrete>()
{
    initMethods();
}


CharacterHistory::CharacterHistory( const RevBayesCore::CharacterHistoryDiscrete &d) : ModelObject<RevBayesCore::CharacterHistoryDiscrete>( d.clone() )
{
    initMethods();
}


CharacterHistory::CharacterHistory( RevBayesCore::CharacterHistoryDiscrete *d) : ModelObject<RevBayesCore::CharacterHistoryDiscrete>( d )
{
    initMethods();
}


CharacterHistory::CharacterHistory( RevBayesCore::TypedDagNode<RevBayesCore::CharacterHistoryDiscrete> *d) : ModelObject<RevBayesCore::CharacterHistoryDiscrete>( d )
{
    initMethods();
}



CharacterHistory::~CharacterHistory()
{
}



CharacterHistory* CharacterHistory::clone() const
{
    return new CharacterHistory( *this );
}


/* Map calls to member methods */
RevPtr<RevVariable> CharacterHistory::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    

    return ModelObject<RevBayesCore::CharacterHistoryDiscrete>::executeMethod( name, args, found );
}


/* Get Rev type of object */
const std::string& CharacterHistory::getClassType(void)
{

    static std::string rev_type = "CharacterHistory";

    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& CharacterHistory::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( ModelObject<RevBayesCore::CharacterHistoryDiscrete>::getClassTypeSpec() ) );

    return rev_type_spec;
}


/** Get type spec */
const TypeSpec& CharacterHistory::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


void CharacterHistory::initMethods( void )
{

    // add the DAG node member methods
    // note that this is a safe case because all DAG nodes are member objects
    if ( dag_node != NULL )
    {
        const MethodTable &dagMethods = dynamic_cast<RevMemberObject*>( dag_node )->getMethods();
        methods.insertInheritedMethods( dagMethods );
    }
    
}

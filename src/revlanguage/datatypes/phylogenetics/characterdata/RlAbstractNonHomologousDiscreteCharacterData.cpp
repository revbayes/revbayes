#include "RlAbstractNonHomologousDiscreteCharacterData.h"

#include <cstddef>

#include "MemberProcedure.h"
#include "RlBoolean.h"
#include "RlString.h"
#include "ArgumentRules.h"
#include "MethodTable.h"
#include "RevMemberObject.h"
#include "RevVariable.h"
#include "TypeSpec.h"

namespace RevLanguage { class Argument; }


using namespace RevLanguage;

AbstractNonHomologousDiscreteCharacterData::AbstractNonHomologousDiscreteCharacterData(void) : ModelObject<RevBayesCore::AbstractNonHomologousDiscreteCharacterData>(  ),
    NonHomologousCharacterData( )
{
    
    initMethods();
    
}


AbstractNonHomologousDiscreteCharacterData::AbstractNonHomologousDiscreteCharacterData( const RevBayesCore::AbstractNonHomologousDiscreteCharacterData &d) : ModelObject<RevBayesCore::AbstractNonHomologousDiscreteCharacterData>( d.clone() ),
    NonHomologousCharacterData( )
{
    
    initMethods();
    
}


AbstractNonHomologousDiscreteCharacterData::AbstractNonHomologousDiscreteCharacterData( RevBayesCore::AbstractNonHomologousDiscreteCharacterData *d) : ModelObject<RevBayesCore::AbstractNonHomologousDiscreteCharacterData>( d ),
    NonHomologousCharacterData( )
{
    initMethods();
    
}


AbstractNonHomologousDiscreteCharacterData::AbstractNonHomologousDiscreteCharacterData( RevBayesCore::TypedDagNode<RevBayesCore::AbstractNonHomologousDiscreteCharacterData> *d) : ModelObject<RevBayesCore::AbstractNonHomologousDiscreteCharacterData>( d ),
    NonHomologousCharacterData( )
{
    
    initMethods();
    
}


AbstractNonHomologousDiscreteCharacterData::~AbstractNonHomologousDiscreteCharacterData()
{
    
}



AbstractNonHomologousDiscreteCharacterData* AbstractNonHomologousDiscreteCharacterData::clone() const
{
    return new AbstractNonHomologousDiscreteCharacterData( *this );
}


/* Map calls to member methods */
RevPtr<RevVariable> AbstractNonHomologousDiscreteCharacterData::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    RevPtr<RevVariable> retVal = NULL;
    
    if ( found == true )
    {
        return retVal;
    }
    
    retVal = executeCharacterDataMethod(name, args, found, &this->getValue());
    
    if ( found == true )
    {
        return retVal;
    }
    else if (name == "isHomologous")
    {
        found = true;
        
        bool ih = this->dag_node->getValue().isHomologyEstablished();
        
        return new RevVariable( new RlBoolean(ih) );
    }
    
    return ModelObject<RevBayesCore::AbstractNonHomologousDiscreteCharacterData>::executeMethod( name, args, found );
}


/* Get Rev type of object */
const std::string& AbstractNonHomologousDiscreteCharacterData::getClassType(void)
{
    
    static std::string rev_type = "AbstractNonHomologousDiscreteCharacterData";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& AbstractNonHomologousDiscreteCharacterData::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( ModelObject<RevBayesCore::AbstractNonHomologousDiscreteCharacterData>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/** Get the type spec of this class. We return a member variable because instances might have different element types. */
const TypeSpec& AbstractNonHomologousDiscreteCharacterData::getTypeSpec(void) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}


void AbstractNonHomologousDiscreteCharacterData::initMethods( void )
{
    
    // insert the character data specific methods
    MethodTable charDataMethods = getCharacterDataMethods();
    methods.insertInheritedMethods( charDataMethods );
    
    // add the DAG node member methods
    // note that this is a sage case because all DAG nodes are member objects
    const MethodTable &dagMethods = dynamic_cast<RevMemberObject*>( dag_node )->getMethods();
    methods.insertInheritedMethods( dagMethods );
    
    ArgumentRules* chartypeArgRules            = new ArgumentRules();
    ArgumentRules* ishomologousArgRules        = new ArgumentRules();
    
    
    methods.addFunction( new MemberProcedure( "chartype", RlString::getClassTypeSpec(),      chartypeArgRules           ) );
    methods.addFunction( new MemberProcedure( "isHomologous", RlBoolean::getClassTypeSpec(),     ishomologousArgRules       ) );

}


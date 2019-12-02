#include "ArgumentRule.h"
#include "ModelVector.h"
#include "MemberProcedure.h"
#include "Natural.h"
#include "OptionRule.h"
#include "RbVector.h"
#include "RealPos.h"
#include "RlRelativeNodeAgeWeightedConstraints.h"
#include "RlBoolean.h"
#include "RlString.h"
#include "RlSimplex.h"
#include "RevVariable.h"
#include <vector>

using namespace RevLanguage;

RlRelativeNodeAgeWeightedConstraints::RlRelativeNodeAgeWeightedConstraints(void) : ModelObject<RevBayesCore::RelativeNodeAgeWeightedConstraints>( )
{
    
    ArgumentRules* nameArgRules               = new ArgumentRules();
    
    methods.addFunction( new MemberProcedure( "filename", RlString::getClassTypeSpec(), nameArgRules           ) );
    
}


RlRelativeNodeAgeWeightedConstraints::RlRelativeNodeAgeWeightedConstraints( RevBayesCore::RelativeNodeAgeWeightedConstraints *v) : ModelObject<RevBayesCore::RelativeNodeAgeWeightedConstraints>( v )
{
    
    ArgumentRules* nameArgRules               = new ArgumentRules();
    
    methods.addFunction( new MemberProcedure( "filename", RlString::getClassTypeSpec(), nameArgRules           ) );
    
}


RlRelativeNodeAgeWeightedConstraints::RlRelativeNodeAgeWeightedConstraints( const RevBayesCore::RelativeNodeAgeWeightedConstraints &v) : ModelObject<RevBayesCore::RelativeNodeAgeWeightedConstraints>( v.clone() )
{
    
    ArgumentRules* nameArgRules               = new ArgumentRules();
    
    methods.addFunction( new MemberProcedure( "filename", RlString::getClassTypeSpec(), nameArgRules           ) );
    
}


RlRelativeNodeAgeWeightedConstraints::RlRelativeNodeAgeWeightedConstraints( RevBayesCore::TypedDagNode<RevBayesCore::RelativeNodeAgeWeightedConstraints> *m) : ModelObject<RevBayesCore::RelativeNodeAgeWeightedConstraints>( m )
{
    
    ArgumentRules* nameArgRules               = new ArgumentRules();
    
    methods.addFunction( new MemberProcedure( "filename", RlString::getClassTypeSpec(), nameArgRules           ) );
    
}


RlRelativeNodeAgeWeightedConstraints* RlRelativeNodeAgeWeightedConstraints::clone() const {
    return new RlRelativeNodeAgeWeightedConstraints( *this );
}


/* Map calls to member methods */
RevPtr<RevVariable> RlRelativeNodeAgeWeightedConstraints::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "filename")
    {
        found = true;
        
        std::string name = this->dag_node->getValue().getFilename();
        RlString *n = new RlString(name);
        return new RevVariable( n );
    }    
    
    return ModelObject<RevBayesCore::RelativeNodeAgeWeightedConstraints>::executeMethod( name, args, found );
}


/**
 * Size of the matrix.
 */
size_t RlRelativeNodeAgeWeightedConstraints::getNumberOfConstraints( void ) const
{
    return this->dag_node->getValue().getNumberOfConstraints();
}


/* Get Rev type of object */
const std::string& RlRelativeNodeAgeWeightedConstraints::getClassType(void) {
    
    static std::string rev_type = "RlRelativeNodeAgeWeightedConstraints";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& RlRelativeNodeAgeWeightedConstraints::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/** Get the type spec of this class. We return a member variable because instances might have different element types. */
const TypeSpec& RlRelativeNodeAgeWeightedConstraints::getTypeSpec(void) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}


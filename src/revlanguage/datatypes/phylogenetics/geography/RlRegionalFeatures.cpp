#include <stddef.h>
#include <vector>
#include <iosfwd>
#include <string>

#include "RbVector.h"
#include "RegionalFeatures.h"
#include "GeographicArea.h"
#include "ModelVector.h"
#include "MemberProcedure.h"
#include "Natural.h"
#include "OptionRule.h"
#include "RealPos.h"
#include "RlRegionalFeatures.h"
#include "RlString.h"
#include "RevVariable.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "RbVectorImpl.h"
#include "Real.h"
#include "RegionalFeatures.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RlConstantNode.h"
#include "RlUtils.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;


RlRegionalFeatures::RlRegionalFeatures(void) : ModelObject<RevBayesCore::RegionalFeatures>( )
{

    ArgumentRules* nLayers               = new ArgumentRules();
    
    methods.addFunction( new MemberProcedure( "nLayers", ModelVector<RlString>::getClassTypeSpec(), nLayers ) );

}


RlRegionalFeatures::RlRegionalFeatures( RevBayesCore::RegionalFeatures *v) : ModelObject<RevBayesCore::RegionalFeatures>( v ) {
   
    ArgumentRules* nLayers               = new ArgumentRules();
    
    methods.addFunction( new MemberProcedure( "nLayers", ModelVector<RlString>::getClassTypeSpec(), nLayers ) );
    
}


RlRegionalFeatures::RlRegionalFeatures( RevBayesCore::TypedDagNode<RevBayesCore::RegionalFeatures> *m) : ModelObject<RevBayesCore::RegionalFeatures>( m ) {
    
    ArgumentRules* nLayers               = new ArgumentRules();
    
    methods.addFunction( new MemberProcedure( "nLayers", ModelVector<RlString>::getClassTypeSpec(), nLayers ) );

}


RlRegionalFeatures* RlRegionalFeatures::clone() const
{
    return new RlRegionalFeatures( *this );
}


/* Map calls to member methods */
RevPtr<RevVariable> RlRegionalFeatures::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "nLayers") {
        found = true;

        return new RevVariable(new Natural((int)this->dag_node->getValue().getNumLayers())) ;
    }
    
    return ModelObject<RevBayesCore::RegionalFeatures>::executeMethod( name, args, found );
}


/* Get Rev type of object */
const std::string& RlRegionalFeatures::getClassType(void) {
    
    static std::string rev_type = "RlRegionalFeatures";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& RlRegionalFeatures::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/** Get the type spec of this class. We return a member variable because instances might have different element types. */
const TypeSpec& RlRegionalFeatures::getTypeSpec(void) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}


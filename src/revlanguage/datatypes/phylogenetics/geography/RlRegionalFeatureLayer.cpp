//
//  RlRegionalFeatureLayer.cpp
//  revbayes-tensorphylo-proj
//
//  Created by Michael Landis on 7/10/22.
//  Copyright © 2022 Michael Landis. All rights reserved.
//

#include "RlRegionalFeatureLayer.h"
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
#include "RlRegionalFeatureLayer.h"
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


RlRegionalFeatureLayer::RlRegionalFeatureLayer(void) : ModelObject<RevBayesCore::RegionalFeatureLayer>( )
{

    ArgumentRules* nLayers               = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nLayers", ModelVector<RlString>::getClassTypeSpec(), nLayers ) );
    
    ArgumentRules* get_args = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "get", ModelVector<ModelVector<Natural> >::getClassTypeSpec(), get_args ) );

}


RlRegionalFeatureLayer::RlRegionalFeatureLayer( RevBayesCore::RegionalFeatureLayer *v) : ModelObject<RevBayesCore::RegionalFeatureLayer>( v ) {
   
    ArgumentRules* nLayers               = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nLayers", ModelVector<RlString>::getClassTypeSpec(), nLayers ) );
    
    ArgumentRules* get_args = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "get", ModelVector<ModelVector<Natural> >::getClassTypeSpec(), get_args ) );

}


RlRegionalFeatureLayer::RlRegionalFeatureLayer( RevBayesCore::TypedDagNode<RevBayesCore::RegionalFeatureLayer> *m) : ModelObject<RevBayesCore::RegionalFeatureLayer>( m ) {
    
    ArgumentRules* nLayers               = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nLayers", ModelVector<RlString>::getClassTypeSpec(), nLayers ) );
    
    ArgumentRules* get_args = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "get", ModelVector<ModelVector<Natural> >::getClassTypeSpec(), get_args ) );

}


RlRegionalFeatureLayer* RlRegionalFeatureLayer::clone() const
{
    return new RlRegionalFeatureLayer( *this );
}


/* Map calls to member methods */
RevPtr<RevVariable> RlRegionalFeatureLayer::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "getFeatures") {
        found = true;
//        return new RevVariable(new ModelVector<ModelVector<Natural> >(this->dag_node->getValue().getLayers("within","categorical"))) ;
    }
    
    return ModelObject<RevBayesCore::RegionalFeatureLayer>::executeMethod( name, args, found );
}


/* Get Rev type of object */
const std::string& RlRegionalFeatureLayer::getClassType(void) {
    
    static std::string rev_type = "RlRegionalFeatureLayer";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& RlRegionalFeatureLayer::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/** Get the type spec of this class. We return a member variable because instances might have different element types. */
const TypeSpec& RlRegionalFeatureLayer::getTypeSpec(void) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}


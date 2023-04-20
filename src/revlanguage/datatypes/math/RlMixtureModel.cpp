//
//  RlMixtureModel.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 3/17/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlMemberFunction.h"
#include "RlMixtureModel.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MemberFunction.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "MixtureModel.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlDeterministicNode.h"
#include "RlTypedFunction.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevLanguage { class Argument; }

using namespace RevLanguage;

MixtureModel::MixtureModel(void) : ModelObject<RevBayesCore::MixtureModel>()
{
    initMethods();
}


MixtureModel::MixtureModel( const RevBayesCore::MixtureModel &v) : ModelObject<RevBayesCore::MixtureModel>( v.clone() )
{
    initMethods();
}


MixtureModel::MixtureModel( RevBayesCore::MixtureModel *v) : ModelObject<RevBayesCore::MixtureModel>( v )
{
    initMethods();
}


MixtureModel::MixtureModel( RevBayesCore::TypedDagNode<RevBayesCore::MixtureModel> *m) : ModelObject<RevBayesCore::MixtureModel>( m )
{
    initMethods();
}


MixtureModel* MixtureModel::clone() const
{
    return new MixtureModel( *this );
}


/* Map calls to member methods */
RevPtr<RevVariable> MixtureModel::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    ; // do nothing for now
    return ModelObject<RevBayesCore::MixtureModel>::executeMethod( name, args, found );
}


/* Get Rev type of object */
const std::string& MixtureModel::getClassType(void) {

    static std::string rev_type = "MixtureModel";

	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& MixtureModel::getClassTypeSpec(void) {

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );

	return rev_type_spec;
}


/** Get the type spec of this class. We return a member variable because instances might have different element types. */
const TypeSpec& MixtureModel::getTypeSpec(void) const {

    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}

void MixtureModel::initMethods(void) {
    
    
    // member functions
    ArgumentRules* mixtureModelArgRules = new ArgumentRules();
//    mixtureModelArgRules->push_back( new ArgumentRule( "rate", RealPos::getClassTypeSpec(), "The rate of the process (or duration of the process assuming rate=1).", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
//    mixtureModelArgRules->push_back( new ArgumentRule( "startAge", RealPos::getClassTypeSpec(), "The start age of the process.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(1.0) ) );
//    mixtureModelArgRules->push_back( new ArgumentRule( "endAge", RealPos::getClassTypeSpec(), "The end age of the process.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(0.0) ) );

    // what methods would we WANT to add?
    // - number of states?
    // - number of components?
    // - equilibrium frequencies for component m?
    // - etc.
    // methods.addFunction( new MemberFunction<MixtureModel, ModelVector<ModelVector<RealPos> > >( "getTransitionProbabilities", this, mixtureModelArgRules   ) );
    
}


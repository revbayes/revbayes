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
#include "RlSimplex.h"
#include "RlTree.h"
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
    if (name == "nComponents")
    {
        found = true;

        int s = dag_node->getValue().getNumberOfComponents();
        return new RevVariable( new Natural( s ) );
    }
    else if (name == "nStates")
    {
        found = true;

        int s = dag_node->getValue().getNumberOfStates();
        return new RevVariable( new Natural( s ) );
    }
    else if (name == "weights")
    {
        found = true;

        auto weights = dag_node->getValue().componentProbs();
        return new RevVariable( new Simplex( weights ) );
    }
    else if (name == "rate")
    {
        found = true;

        auto rate = dag_node->getValue().rate();
        if (not rate)
            throw RbException()<<"Cannot call .rate() on this mixture model: the rate is not defined.";

        return new RevVariable( new RealPos( *rate ) );
    }
    else if (name == "rootFrequencies")
    {
        found = true;

        const Natural& index = static_cast<const Natural&>( args[0].getVariable()->getRevObject() );
        int i = index.getValue() - 1;
        
        auto weights = dag_node->getValue().getRootFrequencies(i);
        return new RevVariable( new Simplex( weights ) );
    }

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

void MixtureModel::initMethods(void)
{
    // add method for call "nComponents" as a function
    ArgumentRules* nComponentsArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nComponents", Natural::getClassTypeSpec(), nComponentsArgRules) );

    // add method for call "nStates" as a function
    ArgumentRules* nStatesArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nStates", Natural::getClassTypeSpec(), nStatesArgRules) );

    // add method for call "weights" as a function
    ArgumentRules* nWeightsArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "weights", Simplex::getClassTypeSpec(), nWeightsArgRules) );

    // add method for call "weights" as a function
    ArgumentRules* rateArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "rate", RealPos::getClassTypeSpec(), rateArgRules) );

    // add method for call "rootFrequencies" as a function
    ArgumentRules* rootFrequenciesArgRules = new ArgumentRules();
    rootFrequenciesArgRules->push_back( new ArgumentRule( "index", Natural::getClassTypeSpec(), "The mixture component index.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "rootFrequencies", Simplex::getClassTypeSpec(), rootFrequenciesArgRules) );

    // add method for call "getTransitionProbabilities" as a function
    ArgumentRules* getTransitionProbArgRules = new ArgumentRules();
    getTransitionProbArgRules->push_back( new ArgumentRule( "tree", Tree::getClassTypeSpec(), "The mixture component index.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    getTransitionProbArgRules->push_back( new ArgumentRule( "node", Natural::getClassTypeSpec(), "The mixture component index.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    getTransitionProbArgRules->push_back( new ArgumentRule( "rate", RealPos::getClassTypeSpec(), "The rate of the process (or duration of the process assuming rate=1).", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberFunction<MixtureModel, ModelVector<ModelVector<ModelVector<RealPos> > > >( "getTransitionProbabilities", this, getTransitionProbArgRules   ) );
}


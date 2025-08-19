//
//  RlRegionalFeatureLayer.cpp
//  revbayes-tensorphylo-proj
//
//  Created by Michael Landis on 7/10/22.
//  Copyright © 2022 Michael Landis. All rights reserved.
//

#include "RlRegionalFeatureLayer.h"
#include <cstddef>
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
#include "RbVector.h"
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
#include "OptionRule.h"
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


namespace RevLanguage { class Argument; }

using namespace RevLanguage;

RlRegionalFeatureLayer::RlRegionalFeatureLayer(void) : ModelObject<RevBayesCore::RegionalFeatureLayer>( )
{
    initMethods();
}

/** Construct from core reference */
RlRegionalFeatureLayer::RlRegionalFeatureLayer(const RevBayesCore::RegionalFeatureLayer &m) : ModelObject<RevBayesCore::RegionalFeatureLayer>( new RevBayesCore::RegionalFeatureLayer( m ) )
{
    initMethods();
}

RlRegionalFeatureLayer::RlRegionalFeatureLayer( RevBayesCore::RegionalFeatureLayer *v) : ModelObject<RevBayesCore::RegionalFeatureLayer>( v ) {
    initMethods();
}


RlRegionalFeatureLayer::RlRegionalFeatureLayer( RevBayesCore::TypedDagNode<RevBayesCore::RegionalFeatureLayer> *m) : ModelObject<RevBayesCore::RegionalFeatureLayer>( m ) {
    initMethods();
}


RlRegionalFeatureLayer* RlRegionalFeatureLayer::clone() const
{
    return new RlRegionalFeatureLayer( *this );
}


/* Map calls to member methods */
RevPtr<RevVariable> RlRegionalFeatureLayer::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "get")
    {
        found = true;
        
        const RevBayesCore::RegionalFeatureLayer& layer = this->getDagNode()->getValue();
        std::string relationship = layer.getFeatureRelationship();
        std::string type = layer.getFeatureType();
        
        if (relationship == "within" && type == "categorical") {
            std::vector<std::vector<double> > val = layer.getFeatureValues();
            RevBayesCore::RbVector<RevBayesCore::RbVector<std::int64_t> > z(1);
            for (size_t i = 0; i < val[0].size(); i++) {
                z[0].push_back((std::int64_t)val[0][i]);
            }
            return new RevVariable( new ModelVector<ModelVector<Natural> >( z ) );
            
        } else if (relationship == "within" && type == "quantitative") {
            std::vector<std::vector<double> > val = layer.getFeatureValues();
            RevBayesCore::RbVector<RevBayesCore::RbVector<double> > z(1);
            for (size_t i = 0; i < val[0].size(); i++) {
                z[0].push_back((double)val[0][i]);
            }
            return new RevVariable( new ModelVector<ModelVector<Real> >( z ) );
        } else if (relationship == "between" && type == "categorical") {
            std::vector<std::vector<double> > val = layer.getFeatureValues();
            RevBayesCore::RbVector<RevBayesCore::RbVector<std::int64_t> > z;
            for (size_t i = 0; i < val.size(); i++) {
                z.push_back(RevBayesCore::RbVector<std::int64_t>());
                for (size_t j = 0; j < val[i].size(); j++) {
                    z[i].push_back((std::int64_t)val[i][j]);
                }
            }
            const ModelVector<ModelVector<Natural> >& x = static_cast<const ModelVector<ModelVector<Natural> >&>(z);
            return new RevVariable( new ModelVector<ModelVector<Natural> >( x ) );
        } else if (relationship == "between" && type == "quantitative") {
            std::vector<std::vector<double> > val = layer.getFeatureValues();
            RevBayesCore::RbVector<RevBayesCore::RbVector<double> > z;
            for (size_t i = 0; i < val.size(); i++) {
                z.push_back(val[i]);
            }
            const ModelVector<ModelVector<Real> >& x = static_cast<const ModelVector<ModelVector<Real> >&>(z);
            return new RevVariable( new ModelVector<ModelVector<Real> >( x ) );
        }
    }
    
    return ModelObject<RevBayesCore::RegionalFeatureLayer>::executeMethod( name, args, found );
}

void RlRegionalFeatureLayer::initMethods(void) {
    ArgumentRules* nLayers = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nLayers", ModelVector<RlString>::getClassTypeSpec(), nLayers ) );
    
    ArgumentRules* get_args = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "get", ModelVector<ModelVector<Natural> >::getClassTypeSpec(), get_args ) );

    return;
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


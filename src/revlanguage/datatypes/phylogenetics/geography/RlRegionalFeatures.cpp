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
#include "RbVector.h"
#include "RealPos.h"
#include "RlRegionalFeatures.h"
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


RlRegionalFeatures::RlRegionalFeatures(void) : ModelObject<RevBayesCore::RegionalFeatures>( )
{
    initMethods();
}

RlRegionalFeatures::RlRegionalFeatures( RevBayesCore::RegionalFeatures *v) : ModelObject<RevBayesCore::RegionalFeatures>( v ) {
    initMethods();
}

/** Construct from core reference */
RlRegionalFeatures::RlRegionalFeatures(const RevBayesCore::RegionalFeatures &m) : ModelObject<RevBayesCore::RegionalFeatures>( new RevBayesCore::RegionalFeatures( m ) )
{
    initMethods();
}

RlRegionalFeatures::RlRegionalFeatures( RevBayesCore::TypedDagNode<RevBayesCore::RegionalFeatures> *m) : ModelObject<RevBayesCore::RegionalFeatures>( m ) {
    initMethods();
}

void RlRegionalFeatures::initMethods(void) {
   
    std::vector<std::string> relationshipOptions;
    relationshipOptions.push_back( "within" );
    relationshipOptions.push_back( "between" );
    std::vector<std::string> typeOptions;
    typeOptions.push_back( "categorical" );
    typeOptions.push_back( "quantitative" );
    
    ArgumentRules* numTimeArgRules   = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "numTimeslices", Natural::getClassTypeSpec(), numTimeArgRules ) );
    
    ArgumentRules* numLayersArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "numLayers", ModelVector<Natural>::getClassTypeSpec(), numLayersArgRules ) );
    
    ArgumentRules* normalizeArgRules = new ArgumentRules();
    normalizeArgRules->push_back( new OptionRule( "relationship", new RlString("within"), relationshipOptions, "" ) );
    methods.addFunction( new MemberProcedure( "normalize", RlUtils::Void, normalizeArgRules ) );
    
    ArgumentRules* getArgRules = new ArgumentRules();
    getArgRules->push_back( new OptionRule( "relationship", new RlString("within"), relationshipOptions, "" ) );
    getArgRules->push_back( new OptionRule( "type", new RlString("categorical"), typeOptions, "" ) );
    getArgRules->push_back( new ArgumentRule( "timeIndex", Natural::getClassTypeSpec(), "Index of time slice.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    
    methods.addFunction( new MemberProcedure( "get", ModelVector<ModelVector<Natural> >::getClassTypeSpec(), getArgRules ) );

    return;
}

RlRegionalFeatures* RlRegionalFeatures::clone() const
{
    return new RlRegionalFeatures( *this );
}


/* Map calls to member methods */
RevPtr<RevVariable> RlRegionalFeatures::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "numLayers") {
        found = true;
        std::map<std::string, std::map<std::string, size_t> > val = this->dag_node->getValue().getNumLayers();
        std::vector<std::int64_t> x;
        x.push_back( val["within"]["categorical"] );
        x.push_back( val["within"]["quantitative"] );
        x.push_back( val["between"]["categorical"] );
        x.push_back( val["between"]["quantitative"] );
        return new RevVariable( new ModelVector<Natural>(x) );
    }
    if (name == "numTimeslices") {
        found = true;
        int val = (int)this->dag_node->getValue().getNumTimeslices();
        return new RevVariable( new Natural(val) );
    }
    if (name == "normalize") {
        found = true;
        
        std::string relationship = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();

        if (relationship == "within") {
            this->dag_node->getValue().normalizeWithinQuantitative();
        } else if (relationship == "between") {
            this->dag_node->getValue().normalizeBetweenQuantitative();
        }

        return NULL;
    }
    if (name == "get")
    {
        found = true;
        std::string relationship = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();
        std::string type = static_cast<const RlString &>( args[1].getVariable()->getRevObject() ).getValue();
        size_t time_index = static_cast<const Natural &>( args[2].getVariable()->getRevObject() ).getValue() - 1;

        // get relevant layer
        const std::vector<RevBayesCore::RegionalFeatureLayer>& y = this->dag_node->getValue().getLayers(relationship, type, time_index);
        return new RevVariable( new ModelVector<RlRegionalFeatureLayer>( y ) );
        
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


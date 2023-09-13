#include <cstddef>
#include <vector>
#include <iosfwd>
#include <string>

#include "RbVector.h"
#include "TimeAtlas.h"
#include "GeographicArea.h"
#include "ModelVector.h"
#include "MemberProcedure.h"
#include "Natural.h"
#include "OptionRule.h"
#include "RealPos.h"
#include "RlAtlas.h"
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
#include "RevObject.h"
#include "RevPtr.h"
#include "RlConstantNode.h"
#include "RlUtils.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

RlAtlas::RlAtlas(void) : ModelObject<RevBayesCore::TimeAtlas>( )
{

    ArgumentRules* nAreasArgRules               = new ArgumentRules();
    ArgumentRules* nEpochsArgRules              = new ArgumentRules();
    ArgumentRules* namesArgRules                = new ArgumentRules();
    ArgumentRules* epochsArgRules                = new ArgumentRules();
    
    methods.addFunction( new MemberProcedure( "names", ModelVector<RlString>::getClassTypeSpec(), namesArgRules ) );
    methods.addFunction( new MemberProcedure( "nAreas", Natural::getClassTypeSpec(), nAreasArgRules ) );
    methods.addFunction( new MemberProcedure( "nEpochs", Natural::getClassTypeSpec(), nEpochsArgRules ) );
    methods.addFunction( new MemberProcedure( "epochTimes", ModelVector<RealPos>::getClassTypeSpec(), epochsArgRules ) );
    
    ArgumentRules* adjacentArgRules             = new ArgumentRules();
    std::vector<std::string> optionsValue;
    optionsValue.push_back( "dispersal" );
    optionsValue.push_back( "dispersal-upper" );
    optionsValue.push_back( "extinction" );
    optionsValue.push_back( "latlon" );
    optionsValue.push_back( "altitude" );
    optionsValue.push_back( "size" );
    adjacentArgRules->push_back( new OptionRule( "value", new RlString( "dispersal" ), optionsValue, "" ) );
    //    adjacentArgRules->push_back( new OptionRule( "elements", new RlString( "off-diagonal" ), optionsElements ) );
    //    adjacentArgRules->push_back( new ArgumentRule("symmetric", RlBoolean::getClassTypeSpec(), ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
    methods.addFunction( new MemberProcedure( "getValues", RlUtils::Void, adjacentArgRules) );

    // Add method for call "size" as a function
    ArgumentRules* sizeArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "size", Natural::getClassTypeSpec(), sizeArgRules) );

}


RlAtlas::RlAtlas( RevBayesCore::TimeAtlas *v) : ModelObject<RevBayesCore::TimeAtlas>( v ),
atlas(v)
{
    ArgumentRules* nAreasArgRules               = new ArgumentRules();
    ArgumentRules* nEpochsArgRules              = new ArgumentRules();
    ArgumentRules* namesArgRules                = new ArgumentRules();
    ArgumentRules* epochsArgRules                = new ArgumentRules();
    
    methods.addFunction( new MemberProcedure( "names", ModelVector<RlString>::getClassTypeSpec(), namesArgRules ) );
    methods.addFunction( new MemberProcedure( "nAreas", Natural::getClassTypeSpec(), nAreasArgRules ) );
    methods.addFunction( new MemberProcedure( "nEpochs", Natural::getClassTypeSpec(), nEpochsArgRules ) );
    methods.addFunction( new MemberProcedure( "epochTimes", ModelVector<RealPos>::getClassTypeSpec(), epochsArgRules ) );
    
    ArgumentRules* adjacentArgRules             = new ArgumentRules();
    std::vector<std::string> optionsValue;
    optionsValue.push_back( "dispersal" );
    optionsValue.push_back( "dispersal-upper" );
    optionsValue.push_back( "extinction" );
    optionsValue.push_back( "latlon" );
    optionsValue.push_back( "altitude" );
    optionsValue.push_back( "size" );
    adjacentArgRules->push_back( new OptionRule( "value", new RlString( "dispersal" ), optionsValue, "" ) );
    //    adjacentArgRules->push_back( new OptionRule( "elements", new RlString( "off-diagonal" ), optionsElements ) );
    //    adjacentArgRules->push_back( new ArgumentRule("symmetric", RlBoolean::getClassTypeSpec(), ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
    methods.addFunction( new MemberProcedure( "getValues", RlUtils::Void, adjacentArgRules) );
    
    // Add method for call "size" as a function
    ArgumentRules* sizeArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "size", Natural::getClassTypeSpec(), sizeArgRules) );
    
}


RlAtlas::RlAtlas( RevBayesCore::TypedDagNode<RevBayesCore::TimeAtlas> *m) : ModelObject<RevBayesCore::TimeAtlas>( m ),
atlas(&m->getValue())
{
    ArgumentRules* nAreasArgRules               = new ArgumentRules();
    ArgumentRules* nEpochsArgRules              = new ArgumentRules();
    ArgumentRules* namesArgRules                = new ArgumentRules();
    ArgumentRules* epochsArgRules                = new ArgumentRules();
    
    methods.addFunction( new MemberProcedure( "names", ModelVector<RlString>::getClassTypeSpec(), namesArgRules ) );
    methods.addFunction( new MemberProcedure( "nAreas", Natural::getClassTypeSpec(), nAreasArgRules ) );
    methods.addFunction( new MemberProcedure( "nEpochs", Natural::getClassTypeSpec(), nEpochsArgRules ) );
    methods.addFunction( new MemberProcedure( "epochTimes", ModelVector<RealPos>::getClassTypeSpec(), epochsArgRules ) );
    
    ArgumentRules* adjacentArgRules             = new ArgumentRules();
    std::vector<std::string> optionsValue;
    optionsValue.push_back( "dispersal" );
    optionsValue.push_back( "dispersal-upper" );
    optionsValue.push_back( "extinction" );
    optionsValue.push_back( "latlon" );
    optionsValue.push_back( "altitude" );
    optionsValue.push_back( "size" );
    adjacentArgRules->push_back( new OptionRule( "value", new RlString( "dispersal" ), optionsValue, "" ) );
    //    adjacentArgRules->push_back( new OptionRule( "elements", new RlString( "off-diagonal" ), optionsElements ) );
    //    adjacentArgRules->push_back( new ArgumentRule("symmetric", RlBoolean::getClassTypeSpec(), ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
    methods.addFunction( new MemberProcedure( "getValues", RlUtils::Void, adjacentArgRules) );
    
    // Add method for call "size" as a function
    ArgumentRules* sizeArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "size", Natural::getClassTypeSpec(), sizeArgRules) );
    
}


RlAtlas* RlAtlas::clone() const
{
    return new RlAtlas( *this );
}


/* Map calls to member methods */
RevPtr<RevVariable> RlAtlas::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "nAreas")
    {
        found = true;
        
        return new RevVariable(new Natural((int)this->dag_node->getValue().getNumAreas())) ;
    }
    else if (name == "nEpochs")
    {
        found = true;
        
        return new RevVariable(new Natural((int)this->dag_node->getValue().getNumEpochs())) ;
    }
    else if (name == "names")
    {
        found = true;
        
        ModelVector<RlString> *n = new ModelVector<RlString>();
        const std::vector<std::vector<RevBayesCore::GeographicArea*> >& areas = this->dag_node->getValue().getAreas();
        for (size_t i = 0; i < areas[0].size(); ++i)
        {
            std::string name = areas[0][i]->getName();
            n->push_back( name );
        }
        return new RevVariable( n );
    }
    else if (name == "epochTimes")
    {
        found = true;
        ModelVector<RealPos> *n = new ModelVector<RealPos>( this->dag_node->getValue().getEpochs() );
        return new RevVariable( n );
    }
    else if (name == "getValues")
    {
        found = true;
        
        // get the member with give index
        std::string value = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();
        std::vector<std::vector<RevBayesCore::GeographicArea*> > areas = this->dag_node->getValue().getAreas();        
        ModelVector<ModelVector<ModelVector<Real > > > *f = new ModelVector<ModelVector<ModelVector<Real > > >();
        for (size_t i = 0; i < areas.size(); i++)
        {
            RevBayesCore::RbVector<RevBayesCore::RbVector<double> > v;
            for (size_t j = 0; j < areas[i].size(); j++)
            {
                if (value == "dispersal")
                    v.push_back(areas[i][j]->getDispersalValues());
                else if (value == "dispersal-upper")
                    v.push_back(areas[i][j]->getDispersalValues());
                else if (value == "extinction")
                    v.push_back(areas[i][j]->getExtinctionValues());
                else if (value == "latlon")
                    v.push_back(areas[i][j]->getLatlon());
                else if (value == "altitude")
                    v.push_back( RevBayesCore::RbVector<double>(1, areas[i][j]->getAltitude()) );
                else if (value == "size")
                    v.push_back( RevBayesCore::RbVector<double>(1, areas[i][j]->getSize()) );
            }
            f->push_back( ModelVector<ModelVector<Real> >(v) );
        }
        return new RevVariable(f);
    }
    
    return ModelObject<RevBayesCore::TimeAtlas>::executeMethod( name, args, found );
}


/* Get Rev type of object */
const std::string& RlAtlas::getClassType(void) {
    
    static std::string rev_type = "RlAtlas";
    
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& RlAtlas::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/** Get the type spec of this class. We return a member variable because instances might have different element types. */
const TypeSpec& RlAtlas::getTypeSpec(void) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}


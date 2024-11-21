#include "RlDiscretizedContinuousCharacterData.h"

#include <cstddef>

#include "ModelVector.h"
#include "RlDiscretizedContinuousState.h"
#include "DiscretizedContinuousCharacterData.h"
#include "ConstantNode.h"
#include "RlMemberFunction.h"
#include "Natural.h"
#include "Real.h"
#include "RealPos.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlContinuousTaxonData.h"
#include "RlDistanceMatrix.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ContinuousTaxonData.h"
#include "DistanceMatrix.h"
#include "MemberFunction.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "RbException.h"
#include "RevMemberObject.h"
#include "RevVariable.h"
#include "RlAbstractDiscreteTaxonData.h"
#include "RlDeterministicNode.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"

using namespace RevLanguage;

/** Default constructor */
DiscretizedContinuousCharacterData::DiscretizedContinuousCharacterData(void) : ModelObject<RevBayesCore::DiscretizedContinuousCharacterData>(),
    HomologousCharacterData( )
{
    
    initMethods();

}

/** Construct from core data type */
DiscretizedContinuousCharacterData::DiscretizedContinuousCharacterData(const RevBayesCore::DiscretizedContinuousCharacterData &d) : ModelObject<RevBayesCore::DiscretizedContinuousCharacterData>( d.clone() ),
    HomologousCharacterData( )
{
    
    
    initMethods();
    
}

/** Construct from core data type */
DiscretizedContinuousCharacterData::DiscretizedContinuousCharacterData(RevBayesCore::DiscretizedContinuousCharacterData *d) : ModelObject<RevBayesCore::DiscretizedContinuousCharacterData>( d ),
    HomologousCharacterData( )
{
    
    initMethods();
}


DiscretizedContinuousCharacterData::DiscretizedContinuousCharacterData( RevBayesCore::TypedDagNode<RevBayesCore::DiscretizedContinuousCharacterData> *d) : ModelObject<RevBayesCore::DiscretizedContinuousCharacterData>( d ),
    HomologousCharacterData( )
{
    
    initMethods();

}



DiscretizedContinuousCharacterData::~DiscretizedContinuousCharacterData()
{
}


void DiscretizedContinuousCharacterData::concatenate(const RevObject &d, std::string type) const
{
    const DiscretizedContinuousCharacterData* tmp = dynamic_cast<const DiscretizedContinuousCharacterData*>( &d );
    if ( tmp != NULL )
    {
        concatenate( *tmp, type );
    }
    else
    {
        throw RbException() << "Cannot add an object of type '" <<  d.getType() << "' to a DiscretizedContinuousCharacterData object.";
    }
}



void DiscretizedContinuousCharacterData::concatenate(const DiscretizedContinuousCharacterData &d, std::string type) const
{
    
    // we need to make this a constant DAG node so that we can actually modify the value
    // otherwise the value might be overwritten again, e.g., if this is a deterministic node.
    //    clone_obj->makeConstantValue();

    // now concatenate
    getDagNode()->getValue().concatenate( d.getValue(), type );
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
DiscretizedContinuousCharacterData* DiscretizedContinuousCharacterData::clone(void) const
{
    
	return new DiscretizedContinuousCharacterData(*this);
}


/* Map calls to member methods */
RevPtr<RevVariable> DiscretizedContinuousCharacterData::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    RevPtr<RevVariable> retVal = dynamic_cast<RevMemberObject *>( dag_node )->executeMethod(name, args, found);
    
    if ( found == true )
    {
        return retVal;
    }
    
    retVal = executeCharacterDataMethod(name, args, found, &this->getValue());
    
    if ( found == true )
    {
        return retVal;
    }
    else if (name == "[]")
    {
        found = true;

		// get the member with given index
        const Natural& index = static_cast<const Natural&>( args[0].getVariable()->getRevObject() );

        if (this->dag_node->getValue().getNumberOfTaxa() < (size_t)(index.getValue()) )
        {
            throw RbException("Index out of bounds in []");
        }

        const RevBayesCore::AbstractDiscreteTaxonData& element = dag_node->getValue().getTaxonData(size_t(index.getValue()) - 1);

        return new RevVariable( new AbstractDiscreteTaxonData( element.clone() ) );
    }
    else if (name == "getDx")
    {
        found = true;

        if( args[0].getVariable()->getRevObject() == RevNullObject::getInstance() )
        {
            // if the argument is null, return the whole vector
        	std::vector<double> dx = dag_node->getValue().getDeltaX();
        	return new RevVariable( new ModelVector<RealPos>(dx) );
        }
        else
        {
    		// get the member with given index
            size_t index = static_cast<const Natural&>( args[0].getVariable()->getRevObject() ).getValue() - 1;
            double dx = dag_node->getValue().getDeltaXForCharacter(index);
            return new RevVariable( new RealPos(dx) );
        }

    }
    else if (name == "getPoints")
    {
        found = true;

		// get the member with given index
        size_t index = static_cast<const Natural&>( args[0].getVariable()->getRevObject() ).getValue() - 1;
        std::vector<double> pts = dag_node->getValue().getPointsForCharacter(index);
        return new RevVariable( new ModelVector<Real>(pts) );
    }


    return ModelObject<RevBayesCore::DiscretizedContinuousCharacterData >::executeMethod( name, args, found );
    
}


/** Get Rev type of object */
const std::string& DiscretizedContinuousCharacterData::getClassType(void)
{
    
    static std::string rev_type = "DiscretizedContinuousCharacterData";
    
	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& DiscretizedContinuousCharacterData::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( ModelObject<RevBayesCore::AbstractHomologousDiscreteCharacterData>::getClassTypeSpec() ) );
    
	return rev_type_spec;
}



/** Get type spec */
const TypeSpec& DiscretizedContinuousCharacterData::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}


void DiscretizedContinuousCharacterData::initMethods( void )
{
    
    // add the DAG node member methods
    // note that this is a safe case because all DAG nodes are member objects
    if ( dag_node != NULL )
    {
        const MethodTable &dagMethods = dynamic_cast<RevMemberObject*>( dag_node )->getMethods();
        methods.insertInheritedMethods( dagMethods );
    }
    
    // insert the character data specific methods
    MethodTable charDataMethods = getCharacterDataMethods();
    methods.insertInheritedMethods( charDataMethods );
    
    ArgumentRules* squareBracketArgRules = new ArgumentRules();
    squareBracketArgRules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the taxon.", ArgumentRule::BY_VALUE, ArgumentRule::ANY  ) );
    methods.addFunction( new MemberProcedure( "[]", AbstractDiscreteTaxonData::getClassTypeSpec(), squareBracketArgRules) );

    ArgumentRules* getDxArgRules = new ArgumentRules();
    getDxArgRules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
    methods.addFunction( new MemberProcedure( "getDx", RealPos::getClassTypeSpec(), getDxArgRules) );

    ArgumentRules* getPointsArgRules = new ArgumentRules();
    getPointsArgRules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "getPoints", ModelVector<Real>::getClassTypeSpec(), getPointsArgRules) );

}




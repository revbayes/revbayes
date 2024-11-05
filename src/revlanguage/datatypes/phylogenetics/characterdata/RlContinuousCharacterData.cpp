#include "RlContinuousCharacterData.h"

#include <cstddef>

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
#include "RlDeterministicNode.h"
#include "RlDiscretizedContinuousCharacterData.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"

using namespace RevLanguage;

/** Default constructor */
ContinuousCharacterData::ContinuousCharacterData(void) : ModelObject<RevBayesCore::ContinuousCharacterData>(),
    HomologousCharacterData( )
{
    
    initMethods();

}

/** Construct from core data type */
ContinuousCharacterData::ContinuousCharacterData(const RevBayesCore::ContinuousCharacterData &d) : ModelObject<RevBayesCore::ContinuousCharacterData>( d.clone() ),
    HomologousCharacterData( )
{
    
    
    initMethods();
    
}

/** Construct from core data type */
ContinuousCharacterData::ContinuousCharacterData(RevBayesCore::ContinuousCharacterData *d) : ModelObject<RevBayesCore::ContinuousCharacterData>( d ),
    HomologousCharacterData( )
{
    
    initMethods();
}


ContinuousCharacterData::ContinuousCharacterData( RevBayesCore::TypedDagNode<RevBayesCore::ContinuousCharacterData> *d) : ModelObject<RevBayesCore::ContinuousCharacterData>( d ),
    HomologousCharacterData( )
{
    
    initMethods();

}



ContinuousCharacterData::~ContinuousCharacterData()
{
}


void ContinuousCharacterData::concatenate(const RevObject &d, std::string type) const
{
    const ContinuousCharacterData* tmp = dynamic_cast<const ContinuousCharacterData*>( &d );
    if ( tmp != NULL )
    {
        concatenate( *tmp, type );
    }
    else
    {
        throw RbException("Cannot add an object of type '" + d.getType() + "' to a ContinuousCharacterData object.");
    }
}



void ContinuousCharacterData::concatenate(const ContinuousCharacterData &d, std::string type) const
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
ContinuousCharacterData* ContinuousCharacterData::clone(void) const
{
    
	return new ContinuousCharacterData(*this);
}


/* Map calls to member methods */
RevPtr<RevVariable> ContinuousCharacterData::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
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
    else if (name == "[]" || name == "getTaxon")
    {
        found = true;
        
        if ( args[0].getVariable()->getRevObject().isType( Natural::getClassTypeSpec() ) )
        {
            // get the member with give index
            const Natural& index = static_cast<const Natural&>( args[0].getVariable()->getRevObject() );
        
            if (this->dag_node->getValue().getNumberOfTaxa() < (size_t)(index.getValue()) )
            {
                throw RbException("Index out of bounds in []");
            }
        
            const RevBayesCore::ContinuousTaxonData& element = static_cast< RevBayesCore::ContinuousCharacterData& >( this->dag_node->getValue() ).getTaxonData(size_t(index.getValue()) - 1);
        
            return new RevVariable( new ContinuousTaxonData( new RevBayesCore::ContinuousTaxonData( element ) ) );
        }
        else
        {
            // get the argument for the index
            const std::string &index = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();
            
            const RevBayesCore::ContinuousTaxonData& element = static_cast< RevBayesCore::ContinuousCharacterData& >( this->dag_node->getValue() ).getTaxonData(index);            
            return new RevVariable( new ContinuousTaxonData( new RevBayesCore::ContinuousTaxonData( element ) ) );
        }
    }
    else if ( name == "maxDifference" )
    {
        found = true;
        
        // get the argument for the index
        std::int64_t index = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue() - 1;
        double max = this->dag_node->getValue().getMaxDifference( index );
        
        return new RevVariable( new Real(max) );
    }
    else if ( name == "maxSpeciesDifference" )
    {
        found = true;
        
        // get the argument for the index
        std::int64_t index = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue() - 1;
        double max = this->dag_node->getValue().getMaxSpeciesDifference( index );
        
        return new RevVariable( new Real(max) );
    }
    else if ( name == "meanDifference" )
    {
        found = true;
        
        // get the argument for the index
        std::int64_t index = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue() - 1;
        double mean = this->dag_node->getValue().getMeanDifference( index );
        
        return new RevVariable( new Real(mean) );
    }
    else if ( name == "meanSpeciesDifference" )
    {
        found = true;
        
        // get the argument for the index
        std::int64_t index = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue() - 1;
        double mean = this->dag_node->getValue().getMeanSpeciesDifference( index );
        
        return new RevVariable( new Real(mean) );
    }
    else if ( name == "meanValue" )
    {
        found = true;
        
        // get the argument for the index
        std::int64_t index = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue() - 1;
        double mean = this->dag_node->getValue().getMeanValue( index );
        
        return new RevVariable( new Real(mean) );
    }
    else if ( name == "minDifference" )
    {
        found = true;
        
        // get the argument for the index
        std::int64_t index = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue() - 1;
        double min = this->dag_node->getValue().getMinSpeciesDifference( index );
        
        return new RevVariable( new Real(min) );
    }
    else if ( name == "minSpeciesDifference" )
    {
        found = true;
        
        // get the argument for the index
        std::int64_t index = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue() - 1;
        double min = this->dag_node->getValue().getMinDifference( index );
        
        return new RevVariable( new Real(min) );
    }
    else if ( name == "getPairwiseSpeciesDifference" )
    {
        found = true;
        
        // get the argument for the index
        std::int64_t index = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue() - 1;
        RevBayesCore::DistanceMatrix pd = this->dag_node->getValue().getPairwiseSpeciesDifference( index );
        
        return new RevVariable( new DistanceMatrix(pd) );
    }
    else if ( name == "varDifference" )
    {
        found = true;
        
        // get the argument for the index
        std::int64_t index = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue() - 1;
        double var = this->dag_node->getValue().getVarSpeciesDifference( index );
        
        return new RevVariable( new Real(var) );
    }
    else if ( name == "varSpeciesDifference" )
    {
        found = true;
        
        // get the argument for the index
        std::int64_t index = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue() - 1;
        double var = this->dag_node->getValue().getVarDifference( index );
        
        return new RevVariable( new Real(var) );
    }
    else if ( name == "var" )
    {
        found = true;
        
        // get the argument for the index
        std::int64_t index = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue() - 1;
        double var = this->dag_node->getValue().getVarValue( index );
        
        return new RevVariable( new RealPos(var) );
    }
    else if ( name == "withinSpeciesVar" )
    {
        found = true;
        
        // get the argument for the index
        std::int64_t index_species  = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue() - 1;
        std::int64_t index_site     = static_cast<const Natural &>( args[1].getVariable()->getRevObject() ).getValue() - 1;
        double var = this->dag_node->getValue().getWithinSpeciesVariance( index_species, index_site );
        
        return new RevVariable( new RealPos(var) );
    }
    else if ( name == "discretizeCharacter" )
    {
    	found = true;

    	// get the arguments
        double error    = static_cast<const RealPos &>( args[0].getVariable()->getRevObject() ).getValue();
        std::int64_t num_bins   = static_cast<const Natural &>( args[1].getVariable()->getRevObject() ).getValue();
        double cushion  = static_cast<const RealPos &>( args[2].getVariable()->getRevObject() ).getValue();

        // make the discretized data
        RevBayesCore::DiscretizedContinuousCharacterData *disc_data = this->dag_node->getValue().discretizeCharacter(error, num_bins, cushion);

    	return new RevVariable( new DiscretizedContinuousCharacterData(disc_data) );
    }
    
    
    return ModelObject<RevBayesCore::ContinuousCharacterData>::executeMethod( name, args, found );
    
}


/** Get Rev type of object */
const std::string& ContinuousCharacterData::getClassType(void)
{
    
    static std::string rev_type = "ContinuousCharacterData";
    
	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& ContinuousCharacterData::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( ModelObject<RevBayesCore::ContinuousCharacterData>::getClassTypeSpec() ) );
    
	return rev_type_spec;
}



/** Get type spec */
const TypeSpec& ContinuousCharacterData::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}


void ContinuousCharacterData::initMethods( void )
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
    
    // Add method for call "x[]" as a function
    ArgumentRules* squareBracketArgRules = new ArgumentRules();
    squareBracketArgRules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the taxon.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "[]", ContinuousTaxonData::getClassTypeSpec(), squareBracketArgRules) );
    
    ArgumentRules* squareBracketArgRules2 = new ArgumentRules();
    squareBracketArgRules2->push_back( new ArgumentRule( "name" , RlString::getClassTypeSpec(), "The name of the taxon.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "getTaxon", ContinuousTaxonData::getClassTypeSpec(), squareBracketArgRules2) );
    
    ArgumentRules* pairwise_difference_arg_rules = new ArgumentRules();
    pairwise_difference_arg_rules->push_back( new ArgumentRule( "charIndex" , Natural::getClassTypeSpec(), "The index of the character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "getPairwiseSpeciesDifference", DistanceMatrix::getClassTypeSpec(), pairwise_difference_arg_rules) );

    ArgumentRules* max_diff_arg_rules = new ArgumentRules();
    max_diff_arg_rules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "maxDifference", Real::getClassTypeSpec(), max_diff_arg_rules ) );

    ArgumentRules* max_species_diff_arg_rules = new ArgumentRules();
    max_species_diff_arg_rules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "maxSpeciesDifference", Real::getClassTypeSpec(), max_species_diff_arg_rules ) );

    ArgumentRules* mean_value_arg_rules = new ArgumentRules();
    mean_value_arg_rules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "meanValue", Real::getClassTypeSpec(), mean_value_arg_rules ) );

    ArgumentRules* mean_diff_arg_rules = new ArgumentRules();
    mean_diff_arg_rules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "meanDifference", Real::getClassTypeSpec(), mean_diff_arg_rules ) );
    
    ArgumentRules* mean_species_diff_arg_rules = new ArgumentRules();
    mean_species_diff_arg_rules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "meanSpeciesDifference", Real::getClassTypeSpec(), mean_species_diff_arg_rules ) );

    ArgumentRules* min_diff_arg_rules = new ArgumentRules();
    min_diff_arg_rules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "minDifference", Real::getClassTypeSpec(), min_diff_arg_rules ) );

    ArgumentRules* min_species_diff_arg_rules = new ArgumentRules();
    min_species_diff_arg_rules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "minSpeciesDifference", Real::getClassTypeSpec(), min_species_diff_arg_rules ) );

    ArgumentRules* var_diff_arg_rules = new ArgumentRules();
    var_diff_arg_rules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "varDifference", Real::getClassTypeSpec(), var_diff_arg_rules ) );

    ArgumentRules* var_species_diff_arg_rules = new ArgumentRules();
    var_species_diff_arg_rules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "varSpeciesDifference", Real::getClassTypeSpec(), var_species_diff_arg_rules ) );
    
    ArgumentRules* var_arg_rules = new ArgumentRules();
    var_arg_rules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "var", RealPos::getClassTypeSpec(), var_arg_rules ) );
    
    ArgumentRules* within_species_var_arg_rules = new ArgumentRules();
    within_species_var_arg_rules->push_back( new ArgumentRule( "species" , Natural::getClassTypeSpec(), "The index of the species.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    within_species_var_arg_rules->push_back( new ArgumentRule( "site" , Natural::getClassTypeSpec(), "The index of the character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "withinSpeciesVar", RealPos::getClassTypeSpec(), within_species_var_arg_rules ) );

    ArgumentRules* discretize_arg_rules = new ArgumentRules();
    discretize_arg_rules->push_back(new ArgumentRule( "error", RealPos::getClassTypeSpec(), "The amount of error in the measured tip states (currently a global value). NULL means all the probability will be placed in the bin closes to the observed value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(0.0) ) );
    discretize_arg_rules->push_back(new ArgumentRule( "nBins", Natural::getClassTypeSpec(), "The number of discrete bins.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    discretize_arg_rules->push_back(new ArgumentRule( "cushionWidth", RealPos::getClassTypeSpec(), "The amount of padding around the range of observed values.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(5.0) ) );
    methods.addFunction( new MemberProcedure( "discretizeCharacter", AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), discretize_arg_rules) );

    
    // member functions
    ArgumentRules* get_element_arg_rules = new ArgumentRules();
    get_element_arg_rules->push_back( new ArgumentRule( "i", Natural::getClassTypeSpec(), "The index of the taxon.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    get_element_arg_rules->push_back( new ArgumentRule( "j", Natural::getClassTypeSpec(), "The index of the character.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberFunction<ContinuousCharacterData, Real>( "get", this, get_element_arg_rules   ) );

}




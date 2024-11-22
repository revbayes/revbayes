#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "RealPos.h"
#include "RlBirthDeathProcess.h"
#include "RlClade.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "RlTimeTree.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "StochasticNode.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

using namespace RevLanguage;


/**
 * Default constructor.
 * 
 * The default constructor does nothing except allocating the object.
 */
BirthDeathProcess::BirthDeathProcess() : TypedDistribution<TimeTree>() 
{
    
}


/**
 * Get Rev type of object 
 *
 * \return The class' name.
 */
const std::string& BirthDeathProcess::getClassType(void) 
{ 
    
    static std::string rev_type = "BirthDeathProcess";
    
	return rev_type; 
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& BirthDeathProcess::getClassTypeSpec(void) 
{ 
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<TimeTree>::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}



/** 
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the homogeneous birth-death process are:
 * (1) time of the process since the origin.
 * (2) time of the process since the rootAge.
 * (3) the sampling probability.
 * (4) the sampling strategy.
 * (5) the condition.
 * (6) the number of taxa.
 * (7) the taxon names.
 * (8) the clade constraints.
 *
 * \return The member rules.
 */
const MemberRules& BirthDeathProcess::getParameterRules(void) const 
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {
        memberRules.push_back( new ArgumentRule( "rootAge", RealPos::getClassTypeSpec()    , "The time of the process starting at the root, if applicable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        std::vector<std::string> optionsStrategy;
        optionsStrategy.push_back( "uniform" );
        optionsStrategy.push_back( "diversified" );
        memberRules.push_back( new OptionRule( "samplingStrategy", new RlString("uniform"), optionsStrategy, "The sampling strategy of including taxa at the present." ) );

        std::vector<std::string> optionsCondition;
        optionsCondition.push_back( "time" );
        optionsCondition.push_back( "survival" );
        optionsCondition.push_back( "nTaxa" );
        memberRules.push_back( new OptionRule( "condition", new RlString("time"), optionsCondition, "The condition of the process." ) );
        memberRules.push_back( new ArgumentRule( "taxa"  , ModelVector<Taxon>::getClassTypeSpec(), "The taxa used for initialization.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "incompleteClades"  , ModelVector<Clade>::getClassTypeSpec(), "Vector of incompletely sampled clades with number of missing species.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        
        rules_set = true;
    }
    
    return memberRules;
}


/** 
 * Set a member variable.
 * 
 * Sets a member variable with the given name and store the pointer to the variable.
 * The value of the variable might still change but this function needs to be called again if the pointer to
 * the variable changes. The current values will be used to create the distribution object.
 *
 * \param[in]    name     Name of the member variable.
 * \param[in]    var      Pointer to the variable.
 */
void BirthDeathProcess::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) 
{
    
    if ( name == "rootAge" )
    {
        start_age = var;
    }
    else if ( name == "samplingStrategy" ) 
    {
        sampling_strategy = var;
    }
    else if ( name == "taxa" )
    {
        taxa = var;
    }
    else if ( name == "condition" ) 
    {
        condition = var;
    }
    else if ( name == "incompleteClades" )
    {
        incomplete_clades = var;
    }
    else 
    {
        TypedDistribution<TimeTree>::setConstParameter(name, var);
    }
    
}

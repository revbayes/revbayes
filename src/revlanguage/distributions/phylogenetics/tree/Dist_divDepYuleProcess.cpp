#include <math.h>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_divDepYuleProcess.h"
#include "DiversityDependentPureBirthProcess.h"
#include "ModelVector.h"
#include "Natural.h"
#include "OptionRule.h"
#include "RealPos.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "RlTimeTree.h"
#include "StochasticNode.h"
#include "Taxon.h"
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
Dist_divDepYuleProcess::Dist_divDepYuleProcess() : TypedDistribution<TimeTree>() 
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process. 
 */
Dist_divDepYuleProcess* Dist_divDepYuleProcess::clone( void ) const 
{
    return new Dist_divDepYuleProcess(*this);
}


/**
 * Create a new internal distribution object.
 *
 * This function simply dynamically allocates a new internal distribution object that can be 
 * associated with the variable. The internal distribution object is created by calling its
 * constructor and passing the distribution-parameters (other DAG nodes) as arguments of the 
 * constructor. The distribution constructor takes care of the proper hook-ups.
 *
 * \return A new internal distribution object.
 */
RevBayesCore::DiversityDependentPureBirthProcess* Dist_divDepYuleProcess::createDistribution( void ) const 
{
    
    // get the parameters
    
    // the root age
    RevBayesCore::TypedDagNode<double>* ra = static_cast<const RealPos &>( rootAge->getRevObject() ).getDagNode();
    // speciation rate
    RevBayesCore::TypedDagNode<double>* s       = static_cast<const RealPos &>( initialLambda->getRevObject() ).getDagNode();
    // extinction rate
    RevBayesCore::TypedDagNode<std::int64_t>* k          = static_cast<const Natural &>( capacity->getRevObject() ).getDagNode();
    // condition
    const std::string& cond                     = static_cast<const RlString &>( condition->getRevObject() ).getValue();
    // taxon names
    const std::vector<RevBayesCore::Taxon> t = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getDagNode()->getValue();
    
    // create the internal distribution object
    RevBayesCore::DiversityDependentPureBirthProcess*   d = new RevBayesCore::DiversityDependentPureBirthProcess(ra, s, k, cond, t);
    
    return d;
}



/**
 * Get Rev type of object 
 *
 * \return The class' name.
 */
const std::string& Dist_divDepYuleProcess::getClassType( void ) 
{ 
    
    static std::string rev_type = "Dist_divDepYuleProcess";
    
	return rev_type; 
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_divDepYuleProcess::getClassTypeSpec( void ) 
{ 
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<TimeTree>::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_divDepYuleProcess::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "DiversityDependentYule";
    
    return d_name;
}


/** 
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the diversity-dependent pure-birth process are:
 * (1) the initial speciation rate lambda which must be a positive real.
 * (2) the carrying capacity that must be a natural number.
 * (3) all member rules specified by BirthDeathProcess.
 *
 * \return The member rules.
 */
const MemberRules& Dist_divDepYuleProcess::getParameterRules(void) const 
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {
        memberRules.push_back( new ArgumentRule( "lambda"  , RealPos::getClassTypeSpec(), "The initial speciation rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "capacity", Natural::getClassTypeSpec(), "The carrying capacity.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "rootAge" , RealPos::getClassTypeSpec(), "The time of the process since the root, if applicable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        std::vector<std::string> optionsCondition;
        optionsCondition.push_back( "time" );
        optionsCondition.push_back( "survival" );
        optionsCondition.push_back( "nTaxa" );
        memberRules.push_back( new OptionRule( "condition"    , new RlString("survival"), optionsCondition, "The condition of the process." ) );
        memberRules.push_back( new ArgumentRule( "taxa"      , ModelVector<Taxon>::getClassTypeSpec(), "The names of the taxa used for simulation.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        // add the rules from the base class
        const MemberRules &parentRules = TypedDistribution<TimeTree>::getParameterRules();
        memberRules.insert(memberRules.end(), parentRules.begin(), parentRules.end());
        
        rules_set = true;
    }
    
    return memberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_divDepYuleProcess::getTypeSpec( void ) const 
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
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
void Dist_divDepYuleProcess::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) 
{
    
    if ( name == "lambda" ) 
    {
        initialLambda = var;
    }
    else if ( name == "capacity" ) 
    {
        capacity = var;
    }
    else if ( name == "rootAge" )
    {
        rootAge = var;
    }
    else if ( name == "taxa" )
    {
        taxa = var;
    }
    else if ( name == "condition" ) 
    {
        condition = var;
    }
    else {
        TypedDistribution<TimeTree>::setConstParameter(name, var);
    }
    
}

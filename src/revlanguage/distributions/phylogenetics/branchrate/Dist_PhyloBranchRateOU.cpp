#include <stddef.h>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "PhyloBrownianProcess.h"
#include "Dist_PhyloBranchRateOU.h"
#include "Real.h"
#include "RlTimeTree.h"
#include "ArgumentRules.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDistribution.h"
#include "TypeSpec.h"

using namespace RevLanguage;


Dist_PhyloBranchRateOU* Dist_PhyloBranchRateOU::clone( void ) const
{
    return new Dist_PhyloBranchRateOU(*this);
}


RevBayesCore::PhyloBranchRatesOU* Dist_PhyloBranchRateOU::createDistribution( void ) const
{
    // get the parameters

    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    
    RevBayesCore::TypedDagNode<double>* ro  = static_cast<const RealPos&>( root_state->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* si  = static_cast<const RealPos&>( sigma->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* al  = static_cast<const RealPos&>( alpha->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* th  = static_cast<const Real&>(    theta->getRevObject() ).getDagNode();

    RevBayesCore::PhyloBranchRatesOU* d    = new RevBayesCore::PhyloBranchRatesOU( tau, ro, si, al, th );
    
    return d;

}



/* Get Rev type of object */
const std::string& Dist_PhyloBranchRateOU::getClassType(void) {
    
    static std::string rev_type = "Dist_PhyloBranchRateOU";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_PhyloBranchRateOU::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< ModelVector<RealPos> >::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_PhyloBranchRateOU::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
//    a_names.push_back( "PhyloOU" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_PhyloBranchRateOU::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhyloBranchRateOU";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_PhyloBranchRateOU::getParameterRules(void) const
{
    
    static MemberRules dist;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        dist.push_back( new ArgumentRule( "tree"     , TimeTree::getClassTypeSpec(), "The tree along which the rates evolves.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist.push_back( new ArgumentRule( "rootState", RealPos::getClassTypeSpec() , "The value of the root.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist.push_back( new ArgumentRule( "sigma"    , RealPos::getClassTypeSpec() , "The rate of drift/evolution of the Ornstein-Uhlenbeck process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist.push_back( new ArgumentRule( "alpha"    , RealPos::getClassTypeSpec() , "The rate of attraction/selection Ornstein-Uhlenbeck process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Real(0.0) ) );
        dist.push_back( new ArgumentRule( "theta"    , Real::getClassTypeSpec()    , "The optimum value of Ornstein-Uhlenbeck process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Real(0.0) ) );

        rules_set = true;
    }
    
    return dist;
}


const TypeSpec& Dist_PhyloBranchRateOU::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}




/** Print value for user */

 void Dist_PhyloBranchRateOU::printValue(std::ostream& o) const
{
    
    o << " BranchRateOU(";
    
    o << "tau=";
    if ( tree != NULL )
    {
        o << tree->getName();
    }
    else
    {
        o << "?";
    }

    o << ",";
     
    o << "sigma=";
    if ( sigma != NULL )
    {
        o << sigma->getName();
    }
    else
    {
        o << "?";
    }

    o << ",";
     
    o << "alpha=";
    if ( alpha != NULL )
    {
        o << alpha->getName();
    }
    else
    {
        o << "?";
    }
    
    o << "theta=";
    if ( theta != NULL )
    {
        o << theta->getName();
    }
    else
    {
        o << "?";
    }

     o << ")";
}



/** Set a member variable */
void Dist_PhyloBranchRateOU::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "rootState" )
    {
        root_state = var;
    }
    else if ( name == "sigma" )
    {
        sigma = var;
    }
    else if ( name == "alpha" )
    {
        alpha = var;
    }
    else if ( name == "theta" )
    {
        theta = var;
    }
    else
    {
        TypedDistribution< ModelVector<RealPos> >::setConstParameter(name, var);
    }
}



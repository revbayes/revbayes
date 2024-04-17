#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "PhyloBrownianProcess.h"
#include "Dist_PhyloNodeStateOU.h"
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


Dist_PhyloNodeStateOU* Dist_PhyloNodeStateOU::clone( void ) const
{
    return new Dist_PhyloNodeStateOU(*this);
}


RevBayesCore::PhyloNodeStateOU* Dist_PhyloNodeStateOU::createDistribution( void ) const
{
    // get the parameters

    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    
    RevBayesCore::TypedDagNode<double>* r   = static_cast<const RealPos&>( root_state->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* s   = static_cast<const RealPos&>( sigma->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* a   = static_cast<const RealPos&>( alpha->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* th  = static_cast<const Real&>( theta->getRevObject() ).getDagNode();

    RevBayesCore::PhyloNodeStateOU* d    = new RevBayesCore::PhyloNodeStateOU( tau, r, s, a, th );
    
    return d;

}



/* Get Rev type of object */
const std::string& Dist_PhyloNodeStateOU::getClassType(void) {
    
    static std::string rev_type = "Dist_PhyloNodeStateOU";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_PhyloNodeStateOU::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< ModelVector<Real> >::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_PhyloNodeStateOU::getDistributionFunctionAliases( void ) const
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
std::string Dist_PhyloNodeStateOU::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhyloNodeStateOU";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_PhyloNodeStateOU::getParameterRules(void) const
{
    
    static MemberRules dist;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        dist.push_back( new ArgumentRule( "tree"     , TimeTree::getClassTypeSpec(), "The tree along which the continuous character evolves.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist.push_back( new ArgumentRule( "rootState", Real::getClassTypeSpec()    , "The value of the root.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist.push_back( new ArgumentRule( "sigma"    , RealPos::getClassTypeSpec() , "The branch-length multiplier to scale the variance of the Ornstein-Uhlenbeck process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist.push_back( new ArgumentRule( "alpha"    , RealPos::getClassTypeSpec() , "The selection/atraction parameter of the Ornstein-Uhlenbeck process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist.push_back( new ArgumentRule( "theta"    , Real::getClassTypeSpec()    , "The optimum parameter of the Ornstein-Uhlenbeck process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return dist;
}


const TypeSpec& Dist_PhyloNodeStateOU::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}




/** Print value for user */

 void Dist_PhyloNodeStateOU::printValue(std::ostream& o) const
{
    
    o << " dnPhyloNodeStateOU(";
    
    o << "tau=";
    if ( tree != NULL ) {
        o << tree->getName();
    } else {
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
    
    o << ",";
    
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
void Dist_PhyloNodeStateOU::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
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
        TypedDistribution< ModelVector<Real> >::setConstParameter(name, var);
    }
}



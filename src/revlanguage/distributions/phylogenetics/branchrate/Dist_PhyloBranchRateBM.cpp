#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "PhyloBrownianProcess.h"
#include "Dist_PhyloBranchRateBM.h"
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


Dist_PhyloBranchRateBM* Dist_PhyloBranchRateBM::clone( void ) const
{
    return new Dist_PhyloBranchRateBM(*this);
}


RevBayesCore::PhyloBranchRatesBM* Dist_PhyloBranchRateBM::createDistribution( void ) const
{
    // get the parameters

    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    
    RevBayesCore::TypedDagNode<double>* r  = static_cast<const RealPos&>( root_state->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* s  = static_cast<const RealPos&>( sigma->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* dr  = static_cast<const Real&>( drift->getRevObject() ).getDagNode();
    
    RevBayesCore::PhyloBranchRatesBM* d    = new RevBayesCore::PhyloBranchRatesBM( tau, r, s, dr );
    
    return d;

}



/* Get Rev type of object */
const std::string& Dist_PhyloBranchRateBM::getClassType(void) {
    
    static std::string rev_type = "Dist_PhyloBranchRateBM";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_PhyloBranchRateBM::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< ModelVector<RealPos> >::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_PhyloBranchRateBM::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
//    a_names.push_back( "PhyloBM" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_PhyloBranchRateBM::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhyloBranchRateBM";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_PhyloBranchRateBM::getParameterRules(void) const
{
    
    static MemberRules dist;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        dist.push_back( new ArgumentRule( "tree"     , TimeTree::getClassTypeSpec(), "The tree along which the continuous character evolves.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist.push_back( new ArgumentRule( "rootState", RealPos::getClassTypeSpec() , "The value of the root.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist.push_back( new ArgumentRule( "sigma"    , RealPos::getClassTypeSpec() , "The branch-length multiplier to scale the variance of the Brownian motion.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist.push_back( new ArgumentRule( "drift"    , Real::getClassTypeSpec()    , "The drift parameter of the Brownian motion.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Real(0.0) ) );
        rules_set = true;
    }
    
    return dist;
}


const TypeSpec& Dist_PhyloBranchRateBM::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}




/** Print value for user */

 void Dist_PhyloBranchRateBM::printValue(std::ostream& o) const
{
    
    o << " brownian(";
    
    o << "tau=";
    if ( tree != NULL ) {
        o << tree->getName();
    } else {
        o << "?";
    }

     o << ",";
     
     o << "sigma=";
     if ( sigma != NULL ) {
         o << sigma->getName();
     } else {
         o << "?";
     }

     o << ",";
     
     o << "drift=";
     if ( drift != NULL ) {
         o << drift->getName();
     } else {
         o << "?";
     }

     o << ")";
}



/** Set a member variable */
void Dist_PhyloBranchRateBM::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
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
    else if ( name == "drift" )
    {
        drift = var;
    }
    else
    {
        TypedDistribution< ModelVector<RealPos> >::setConstParameter(name, var);
    }
}



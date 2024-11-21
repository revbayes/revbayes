#include <math.h>
#include <cstddef>
#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "Dist_wishart.h"
#include "Natural.h"
#include "RealPos.h"
#include "RlMatrixRealSymmetric.h"
#include "StochasticNode.h"
#include "WishartDistribution.h"
#include "ArgumentRules.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "MatrixReal.h"
#include "RbException.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistribution.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

using namespace RevLanguage;

Dist_wishart::Dist_wishart() : TypedDistribution<MatrixRealSymmetric>()
{
    
}


Dist_wishart::~Dist_wishart()
{
    
}



Dist_wishart* Dist_wishart::clone( void ) const
{
    return new Dist_wishart(*this);
}


RevBayesCore::WishartDistribution* Dist_wishart::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal>* om = NULL;
    RevBayesCore::TypedDagNode<double>* ka = NULL;
    /*
    if (omega != NULL)  {
        om = static_cast<const MatrixRealSymmetric &>( omega->getRevObject() ).getDagNode();
    }
     */
    if (kappa != NULL)
    {
        ka = static_cast<const RealPos&>( kappa->getRevObject() ).getDagNode();
    }
    
    RevBayesCore::TypedDagNode<std::int64_t>* deg = static_cast<const Natural &>( df->getRevObject()).getDagNode();

    RevBayesCore::TypedDagNode<std::int64_t>* dm = NULL;
//    int dm = -1;
    if (dim != NULL)
    {
        dm = static_cast<const Natural &>( dim->getRevObject()).getDagNode();
//        dm = static_cast<const Natural &>( dim->getValue()).getValue();
    }
    RevBayesCore::WishartDistribution* w    =  NULL;
    
    if ( om != NULL )
    {
            w = new RevBayesCore::WishartDistribution( om, deg );
    }
    else
    {
        if (dm == NULL || ka == NULL)
        {
            throw RbException("error in Dist_wishart: should specify arguments");
        }
        w = new RevBayesCore::WishartDistribution( dm, ka, deg );
    }
    return w;
}



/* Get Rev type of object */
const std::string& Dist_wishart::getClassType(void)
{
    
    static std::string rev_type = "Dist_wishart";
    
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_wishart::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_wishart::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "Wishart";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_wishart::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
//        distExpMemberRules.push_back( new ArgumentRule( "omega", true, MatrixRealSymmetric::getClassTypeSpec() ) );
        dist_member_rules.push_back( new ArgumentRule( "df"   , Natural::getClassTypeSpec(), "The degrees of dreedom.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "kappa", RealPos::getClassTypeSpec(), "The scaling parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "dim"  , Natural::getClassTypeSpec(), "The dimension of the distribution.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_wishart::getTypeSpec( void ) const {
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_wishart::printValue(std::ostream& o) const {
    
    o << " Wishart(omega=";
/*
    if ( omega != NULL ) {
        o << omega->getName();
    } else {
*/
 if (kappa != NULL)  {
            if (dim == NULL)    {
                throw RbException("error in Wishart distribution: kappa and dim should both be non null");
            }
            o << kappa->getName() << ".I_" << dim->getName();
        }
        else{
            o << "?";
        }
   // }
    o << ")";
}


/** Set a member variable */
void Dist_wishart::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "omega" )
    {
//        omega = var;
    }
    else if ( name == "kappa" )
    {
        kappa = var;
    }
    else if ( name == "df" )
    {
        df = var;
    }
    else if ( name == "dim" )
    {
        dim = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}

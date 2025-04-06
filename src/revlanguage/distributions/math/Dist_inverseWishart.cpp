#include <cmath>
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "Dist_inverseWishart.h"
#include "Natural.h"
#include "RealPos.h"
#include "ModelVector.h"
#include "RlMatrixRealSymmetric.h"
#include "StochasticNode.h"
#include "InverseWishartDistribution.h"
#include "ArgumentRules.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "MatrixReal.h"
#include "ModelObject.h"
#include "RbException.h"
#include "RevNullObject.h"
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

namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevLanguage;

Dist_inverseWishart::Dist_inverseWishart() : TypedDistribution<MatrixRealSymmetric>()
{
    
}


Dist_inverseWishart::~Dist_inverseWishart()
{
    
}



Dist_inverseWishart* Dist_inverseWishart::clone( void ) const
{
    return new Dist_inverseWishart(*this);
}


RevBayesCore::InverseWishartDistribution* Dist_inverseWishart::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal>* sg = NULL;
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* dv = NULL;
    RevBayesCore::TypedDagNode<double>* ka = NULL;
    RevBayesCore::TypedDagNode<std::int64_t>* deg = NULL;
    RevBayesCore::TypedDagNode<std::int64_t>* dm = NULL;
    
    if ( sigma->getRevObject() != RevNullObject::getInstance() )
    {
        sg = static_cast<const MatrixRealSymmetric &>( sigma->getRevObject() ).getDagNode();
    }
    
    if ( diagonal->getRevObject() != RevNullObject::getInstance() )
    {
        dv = static_cast<const ModelVector<RealPos> &>( diagonal->getRevObject() ).getDagNode();
    }
    
    if ( kappa->getRevObject() != RevNullObject::getInstance() )
    {
        ka = static_cast<const RealPos&>( kappa->getRevObject() ).getDagNode();
    }
    
    if ( df->getRevObject() != RevNullObject::getInstance() )
    {
        deg = static_cast<const Natural &>( df->getRevObject()).getDagNode();
    }

    if ( dim->getRevObject() != RevNullObject::getInstance() )
    {
        dm = static_cast<const Natural &>( dim->getRevObject()).getDagNode();
    }
    
    RevBayesCore::InverseWishartDistribution* w =  NULL;

    if ( sg != NULL && sg->getValue().getDim() != 0 )
    {
        // parameter is sigma
        w = new RevBayesCore::InverseWishartDistribution( sg, deg );
    }
    else if (dm == NULL || dm->getValue() == 0)
    {
        // parameter is Diagonal(kappaVector))
        w = new RevBayesCore::InverseWishartDistribution( dv, deg );
    }
    else
    {
        // parameter is kappa * Id
        w = new RevBayesCore::InverseWishartDistribution( dm, ka, deg );
    }
    
    return w;
}



/* Get class name of object */
const std::string& Dist_inverseWishart::getClassType(void) {
    
    static std::string revClassType = "Dist_inverseWishart";
    
	return revClassType;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_inverseWishart::getClassTypeSpec(void)
{
    
    static TypeSpec revClassTypeSpec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
    
	return revClassTypeSpec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_inverseWishart::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "invWishart" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_inverseWishart::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "InverseWishart";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_inverseWishart::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        dist_member_rules.push_back( new ArgumentRule( "sigma"   , MatrixRealSymmetric::getClassTypeSpec() , "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL  ) );
        dist_member_rules.push_back( new ArgumentRule( "diagonal", ModelVector<RealPos>::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL  ) );
        dist_member_rules.push_back( new ArgumentRule( "df"      , Natural::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "kappa"   , RealPos::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "dim"     , Natural::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_inverseWishart::getTypeSpec( void ) const {
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_inverseWishart::printValue(std::ostream& o) const {
    
    o << " InverseWishart(sigma=";
/*
    if ( sigma != NULL ) {
        o << sigma->getName();
    } else {
*/
    if (kappa != NULL)  {
        if (dim == NULL) {
            throw RbException("error in Wishart distribution: kappa and dim should both be non null");
        }
        o << kappa->getName() << ".I_" << dim->getName();
    } else {
        o << "?";
    }

    o << ")";
}


/** Set a member variable */
void Dist_inverseWishart::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "sigma" ) {
        sigma = var;
    }
    else if ( name == "diagonal" ) {
        diagonal = var;
    }
    else if ( name == "kappa" ) {
        kappa = var;
    }
    else if ( name == "df" ) {
        df = var;
    }
    else if ( name == "dim" ) {
        dim = var;
    }
    else {
        Distribution::setConstParameter(name, var);
    }
}

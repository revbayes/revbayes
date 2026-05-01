#include <cmath>
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "Dist_decomposedInverseWishart.h"
#include "Natural.h"
#include "RealPos.h"
#include "ModelVector.h"
#include "RlMatrixRealSymmetric.h"
#include "StochasticNode.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagMemberFunction.h"
#include "DagNode.h"
#include "DecomposedInverseWishartDistribution.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MatrixReal.h"
#include "ModelObject.h"
#include "RbException.h"
#include "RbVector.h"
#include "Real.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistribution.h"
#include "RlMatrixReal.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

Dist_decomposedInverseWishart::Dist_decomposedInverseWishart(void) : TypedDistribution<MatrixReal>()
{
    
}


Dist_decomposedInverseWishart::~Dist_decomposedInverseWishart(void)
{
    
}


Dist_decomposedInverseWishart* Dist_decomposedInverseWishart::clone(void) const
{

    return new Dist_decomposedInverseWishart(*this);
}


RevBayesCore::DecomposedInverseWishartDistribution* Dist_decomposedInverseWishart::createDistribution(void) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal>* sg = static_cast<const MatrixRealSymmetric &>( sigma->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* dv = static_cast<const ModelVector<RealPos> &>( diagonal->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<double>* ka = static_cast<const RealPos&>( kappa->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<long>* deg = static_cast<const Natural &>( df->getRevObject()).getDagNode();
    RevBayesCore::TypedDagNode<long>* dm = static_cast<const Natural &>( dim->getRevObject()).getDagNode();
    
    RevBayesCore::DecomposedInverseWishartDistribution* w =  0;
    if ( sg->getValue().getDim() == 0 )
    {
        // parameter is sigma
        w = new RevBayesCore::DecomposedInverseWishartDistribution( sg, deg );
    }
    else if (dm->getValue() == 0)
    {
        // parameter is Diagonal(kappaVector))
        w = new RevBayesCore::DecomposedInverseWishartDistribution( dv, deg );
    }
    else
    {
        // parameter is kappa * Id
        w = new RevBayesCore::DecomposedInverseWishartDistribution( dm, ka, deg );
    }
    return w;
}


/* Get class name of object */
const std::string& Dist_decomposedInverseWishart::getClassType(void) {
    
    static std::string revClassType = "Dist_decomposedInverseWishart";
	return revClassType;
}


/* Get class type spec describing type of object */
const TypeSpec& Dist_decomposedInverseWishart::getClassTypeSpec(void)
{
    
    static TypeSpec revClassTypeSpec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
	return revClassTypeSpec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_decomposedInverseWishart::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "DecomposedInvWishart";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_decomposedInverseWishart::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    if ( !rules_set )
    {
        dist_member_rules.push_back( new ArgumentRule( "sigma"   , MatrixRealSymmetric::getClassTypeSpec() , "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new MatrixRealSymmetric()  ) );
        dist_member_rules.push_back( new ArgumentRule( "diagonal", ModelVector<RealPos>::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new ModelVector<RealPos>()  ) );
        dist_member_rules.push_back( new ArgumentRule( "df"      , Natural::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Natural(0L) ) );
        dist_member_rules.push_back( new ArgumentRule( "kappa"   , RealPos::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Real(0.0) ) );
        dist_member_rules.push_back( new ArgumentRule( "dim"     , Natural::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Natural(0L) ) );
        rules_set = true;
    }
    return dist_member_rules;
}


const TypeSpec& Dist_decomposedInverseWishart::getTypeSpec(void) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    return ts;
}


/** Print value for user */
void Dist_decomposedInverseWishart::printValue(std::ostream& o) const {
    
    o << " DecomposedInverseWishart(sigma=";
/*
    if ( sigma != NULL ) {
        o << sigma->getName();
    } else {
*/
    if (kappa != NULL)  {
        if (dim == NULL) {
            throw RbException("error in Decomposed Inversev Wishart distribution: kappa and dim should both be non null");
        }
        o << kappa->getName() << ".I_" << dim->getName();
    } else {
        o << "?";
    }

    o << ")";
}


/** Set a member variable */
void Dist_decomposedInverseWishart::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "sigma" )
        {
        sigma = var;
        }
    else if ( name == "diagonal" )
        {
        diagonal = var;
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

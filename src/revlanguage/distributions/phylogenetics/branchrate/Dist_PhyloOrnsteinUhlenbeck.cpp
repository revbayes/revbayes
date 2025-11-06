/*
 * File:   Dist_branchHeterogeneousOrnsteinUhlenbeck.cpp
 * Author: Bastien Boussau
 *
 * Created on 05/10/14.
 */

#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "Dist_PhyloOrnsteinUhlenbeck.h"
#include "ModelVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RlTimeTree.h"
#include "Tree.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "PhyloOrnsteinUhlenbeckProcess.h"
#include "RbException.h"
#include "RbVector.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlDistribution.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;


Dist_PhyloOrnsteinUhlenbeck* Dist_PhyloOrnsteinUhlenbeck::clone( void ) const {
    return new Dist_PhyloOrnsteinUhlenbeck(*this);
}


RevBayesCore::PhyloOrnsteinUhlenbeckProcess* Dist_PhyloOrnsteinUhlenbeck::createDistribution( void ) const {
    // get the parameters
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();

    size_t nNodes = tau->getValue().getNumberOfNodes();
    RevBayesCore::PhyloOrnsteinUhlenbeckProcess* dist = new RevBayesCore::PhyloOrnsteinUhlenbeckProcess( tau );

    //Sigma:
    if ( sigma->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* s = static_cast<const ModelVector<RealPos> &>( sigma->getRevObject() ).getDagNode();

        // sanity check
        if ( (nNodes-1) != s->getValue().size() )
        {
            throw RbException( "The number of sigma values does not match the number of branches" );
        }
        
        dist->setSigma( s );
    }
    else
    {
        RevBayesCore::TypedDagNode< double >* s = static_cast<const RealPos &>( sigma->getRevObject() ).getDagNode();
        dist->setSigma( s );
    }
    
    //Mean:
    if ( mean->getRevObject().isType( ModelVector<Real>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* m = static_cast<const ModelVector<Real> &>( mean->getRevObject() ).getDagNode();
        
        // sanity check
        if ( (nNodes-1) != m->getValue().size() )
        {
            throw RbException( "The number of mean values does not match the number of branches" );
        }
        
        dist->setMean( m );
    }
    else
    {
        RevBayesCore::TypedDagNode< double >* m = static_cast<const Real &>( mean->getRevObject() ).getDagNode();
        dist->setMean( m );
    }

    //Phi:
    if ( phi->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* p = static_cast<const ModelVector<RealPos> &>( phi->getRevObject() ).getDagNode();
        
        // sanity check
        if ( (nNodes-1) != p->getValue().size() )
        {
            throw RbException( "The number of phi values does not match the number of branches" );
        }
        
        dist->setPhi( p );
    }
    else
    {
        RevBayesCore::TypedDagNode< double >* p = static_cast<const RealPos &>( phi->getRevObject() ).getDagNode();
        dist->setPhi( p );
    }
    
    return dist;
    
}



/* Get Rev type of object */
const std::string& Dist_PhyloOrnsteinUhlenbeck::getClassType(void)
{
    
    static std::string rev_type = "Dist_PhyloOrnsteinUhlenbeck";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_PhyloOrnsteinUhlenbeck::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_PhyloOrnsteinUhlenbeck::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "PhyloOU" );
    
    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_PhyloOrnsteinUhlenbeck::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhyloOrnsteinUhlenbeck";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_PhyloOrnsteinUhlenbeck::getMemberRules(void) const
{
    
    static MemberRules dist;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        dist.push_back( new ArgumentRule( "tree" , TimeTree::getClassTypeSpec(), "The tree along which the character evolves.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        std::vector<TypeSpec> posTypes;
        posTypes.push_back( RealPos::getClassTypeSpec() );
        posTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        dist.push_back( new ArgumentRule( "sigma", posTypes, "The variance parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        std::vector<TypeSpec> types;
        types.push_back( RealPos::getClassTypeSpec() );
        types.push_back( ModelVector<Real>::getClassTypeSpec() );
        dist.push_back( new ArgumentRule( "mean" , types, "The location/mean of the process", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist.push_back( new ArgumentRule( "phi"  , posTypes, "The attraction (speed) to the mean.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return dist;
}


const TypeSpec& Dist_PhyloOrnsteinUhlenbeck::getTypeSpec( void ) const {
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}




/** Print value for user */

void Dist_PhyloOrnsteinUhlenbeck::printValue(std::ostream& o) const {
    
    o << " OrnsteinUhlenbeck(";
    
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
    
    o << "mean=";
    if ( mean != NULL ) {
        o << mean->getName();
    } else {
        o << "?";
    }
    
    o << ",";
    
    o << "phi=";
    if ( phi != NULL ) {
        o << phi->getName();
    } else {
        o << "?";
    }
    
    o << ")";
}



/** Set a member variable */
void Dist_PhyloOrnsteinUhlenbeck::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "sigma" )
    {
        sigma = var;
    }
    else if ( name == "mean" )
    {
        mean = var;
    }
    else if ( name == "phi" )
    {
        phi = var;
    }
    else if ( name == "root" )
    {
        root = var;
    }
    else {
        Distribution::setConstParameter(name, var);
    }
}


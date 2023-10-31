#include "Dist_PhyloMultiSampleOrnsteinUhlenbeck.h"

#include <cstddef>
#include <ostream>

#include "PhyloMultiSampleOrnsteinUhlenbeckProcess.h"
#include "RlTaxon.h"
#include "RlTree.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RbVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDistribution.h"
#include "StringUtilities.h"
#include "Tree.h"
#include "TypeSpec.h"

namespace RevBayesCore { class Taxon; }

using namespace RevLanguage;


Dist_PhyloMultiSampleOrnsteinUhlenbeck::Dist_PhyloMultiSampleOrnsteinUhlenbeck() : TypedDistribution< ContinuousCharacterData >()
{
    
}


Dist_PhyloMultiSampleOrnsteinUhlenbeck::~Dist_PhyloMultiSampleOrnsteinUhlenbeck()
{
    
}



Dist_PhyloMultiSampleOrnsteinUhlenbeck* Dist_PhyloMultiSampleOrnsteinUhlenbeck::clone( void ) const
{
    
    return new Dist_PhyloMultiSampleOrnsteinUhlenbeck(*this);
}


RevBayesCore::TypedDistribution< RevBayesCore::ContinuousCharacterData >* Dist_PhyloMultiSampleOrnsteinUhlenbeck::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const Tree &>( tree->getRevObject() ).getDagNode();
    size_t n = size_t( static_cast<const Natural &>( n_sites->getRevObject() ).getValue() );
    size_t n_nodes = tau->getValue().getNumberOfNodes();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* var = static_cast<const ModelVector<RealPos> &>( within_species_variances->getRevObject() ).getDagNode();
    const std::vector<RevBayesCore::Taxon>      &t  = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();

    RevBayesCore::PhyloMultiSampleOrnsteinUhlenbeckProcess *dist = new RevBayesCore::PhyloMultiSampleOrnsteinUhlenbeckProcess(tau, var, t, n);
    
    // set alpha
    if ( alpha->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* a = static_cast<const ModelVector<RealPos> &>( alpha->getRevObject() ).getDagNode();
        dist->setAlpha( a );
    }
    else
    {
        RevBayesCore::TypedDagNode< double >* a = static_cast<const RealPos &>( alpha->getRevObject() ).getDagNode();
        dist->setAlpha( a );
    }
    
    // set theta
    if ( theta->getRevObject().isType( ModelVector<Real>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* t = static_cast<const ModelVector<Real> &>( theta->getRevObject() ).getDagNode();
        dist->setTheta( t );
    }
    else
    {
        RevBayesCore::TypedDagNode< double >* t = static_cast<const Real &>( theta->getRevObject() ).getDagNode();
        dist->setTheta( t );
    }
    
    // set sigma
    if ( sigma->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* s = static_cast<const ModelVector<RealPos> &>( sigma->getRevObject() ).getDagNode();
        dist->setSigma( s );
    }
    else
    {
        RevBayesCore::TypedDagNode< double >* s = static_cast<const RealPos &>( sigma->getRevObject() ).getDagNode();
        dist->setSigma( s );
    }
    
    
    // set the root states
    //    if ( rootStates->getRevObject().isType( ModelVector<Real>::getClassTypeSpec() ) )
    //    {
    //        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* rs = static_cast<const ModelVector<Real> &>( rootStates->getRevObject() ).getDagNode();
    //        dist->setRootState( rs );
    //    }
    //    else
    //    {
    RevBayesCore::TypedDagNode< double >* rs = static_cast<const Real &>( root_states->getRevObject() ).getDagNode();
    dist->setRootState( rs );
    //    }
    
    return dist;
}



/* Get Rev type of object */
const std::string& Dist_PhyloMultiSampleOrnsteinUhlenbeck::getClassType(void)
{
    
    static std::string rev_type = "Dist_PhyloMultiSampleOrnsteinUhlenbeck";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_PhyloMultiSampleOrnsteinUhlenbeck::getClassTypeSpec(void)
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
std::string Dist_PhyloMultiSampleOrnsteinUhlenbeck::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhyloMultiSampleOrnsteinUhlenbeck";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_PhyloMultiSampleOrnsteinUhlenbeck::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        dist_member_rules.push_back( new ArgumentRule( "tree" , Tree::getClassTypeSpec(), "The tree along which the character evolves.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        std::vector<TypeSpec> alphaTypes;
        alphaTypes.push_back( RealPos::getClassTypeSpec() );
        alphaTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "alpha" , alphaTypes, "The rate of attraction/selection (per branch).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        
        std::vector<TypeSpec> thetaTypes;
        thetaTypes.push_back( Real::getClassTypeSpec() );
        thetaTypes.push_back( ModelVector<Real>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "theta" , thetaTypes, "The optimum value (per branch).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        
        std::vector<TypeSpec> sigmaTypes;
        sigmaTypes.push_back( RealPos::getClassTypeSpec() );
        sigmaTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "sigma" , sigmaTypes, "The rate of random drift (per branch).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        
        dist_member_rules.push_back( new ArgumentRule( "withinSpeciesVariances" , ModelVector<RealPos>::getClassTypeSpec(), "The per species within-species variances.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        std::vector<TypeSpec> rootStateTypes;
        rootStateTypes.push_back( Real::getClassTypeSpec() );
        rootStateTypes.push_back( ModelVector<Real>::getClassTypeSpec() );
        Real *defaultRootStates = new Real(0.0);
        dist_member_rules.push_back( new ArgumentRule( "rootStates" , rootStateTypes, "The vector of root states.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, defaultRootStates ) );
        
        dist_member_rules.push_back( new ArgumentRule( "nSites"         ,  Natural::getClassTypeSpec(), "The number of sites which is used for the initialized (random draw) from this distribution.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(10) ) );
        dist_member_rules.push_back( new ArgumentRule( "taxa"  , ModelVector<Taxon>::getClassTypeSpec(), "The vector of taxa which have species and individual names.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_PhyloMultiSampleOrnsteinUhlenbeck::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_PhyloMultiSampleOrnsteinUhlenbeck::printValue(std::ostream& o) const
{
    
    o << "PhyloMultiSampleOrnsteinUhlenbeck(tree=";
    if ( tree != NULL )
    {
        o << tree->getName();
    }
    else
    {
        o << "?";
    }
    o << ", sigma=";
    if ( sigma != NULL )
    {
        o << sigma->getName();
    }
    else
    {
        o << "?";
    }
    o << ", theta=";
    if ( theta != NULL )
    {
        o << theta->getName();
    }
    else
    {
        o << "?";
    }
    o << ", alpha=";
    if ( alpha != NULL )
    {
        o << alpha->getName();
    }
    else
    {
        o << "?";
    }
    o << ", nSites=";
    if ( n_sites != NULL )
    {
        o << n_sites->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
    
}


/** Set a member variable */
void Dist_PhyloMultiSampleOrnsteinUhlenbeck::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "alpha" )
    {
        alpha = var;
    }
    else if ( name == "theta" )
    {
        theta = var;
    }
    else if ( name == "sigma" )
    {
        sigma = var;
    }
    else if ( name == "withinSpeciesVariances" )
    {
        within_species_variances = var;
    }
    else if ( name == "rootStates" )
    {
        root_states = var;
    }
    else if ( name == "nSites" )
    {
        n_sites = var;
    }
    else if ( name == "taxa" )
    {
        taxa = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
    
}


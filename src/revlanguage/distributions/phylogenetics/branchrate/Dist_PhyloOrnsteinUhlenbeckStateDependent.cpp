#include "Dist_PhyloOrnsteinUhlenbeckStateDependent.h"

#include <stddef.h>
#include <ostream>

#include "PhyloOrnsteinUhlenbeckStateDependent.h"
#include "RlTree.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "Natural.h"
#include "Real.h"
#include "RealPos.h"
#include "RlCharacterHistory.h"
#include "RlDistribution.h"
#include "StringUtilities.h"
#include "Tree.h"
#include "TypeSpec.h"

namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevLanguage;


Dist_PhyloOrnsteinUhlenbeckStateDependent::Dist_PhyloOrnsteinUhlenbeckStateDependent() : TypedDistribution< ContinuousCharacterData >()
{
    
}


Dist_PhyloOrnsteinUhlenbeckStateDependent::~Dist_PhyloOrnsteinUhlenbeckStateDependent()
{
    
}



Dist_PhyloOrnsteinUhlenbeckStateDependent* Dist_PhyloOrnsteinUhlenbeckStateDependent::clone( void ) const
{
    
    return new Dist_PhyloOrnsteinUhlenbeckStateDependent(*this);
}


RevBayesCore::TypedDistribution< RevBayesCore::ContinuousCharacterData >* Dist_PhyloOrnsteinUhlenbeckStateDependent::createDistribution( void ) const
{
    
    // get the parameters
    size_t n = size_t( static_cast<const Natural &>( n_sites->getRevObject() ).getValue() );
    
    const CharacterHistory& rl_char_hist = static_cast<const RevLanguage::CharacterHistory&>( character_history->getRevObject() );
    RevBayesCore::TypedDagNode<RevBayesCore::CharacterHistoryDiscrete>* char_hist   =  rl_char_hist.getDagNode();
    
    RevBayesCore::PhyloOrnsteinUhlenbeckStateDependent *dist = new RevBayesCore::PhyloOrnsteinUhlenbeckStateDependent(char_hist, n);
    
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
    
     // set observation variance
    RevBayesCore::TypedDagNode< double >* ov = static_cast<const RealPos &>( observation_variance->getRevObject() ).getDagNode();
    dist->setObservationVariance( ov );
    
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
const std::string& Dist_PhyloOrnsteinUhlenbeckStateDependent::getClassType(void)
{
    
    static std::string rev_type = "Dist_PhyloOrnsteinUhlenbeckStateDependent";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_PhyloOrnsteinUhlenbeckStateDependent::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_PhyloOrnsteinUhlenbeckStateDependent::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "PhyloOUSD" );
    a_names.push_back( "PhOUSD" );

    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_PhyloOrnsteinUhlenbeckStateDependent::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhyloOrnsteinUhlenbeckStateDependent";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_PhyloOrnsteinUhlenbeckStateDependent::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        dist_member_rules.push_back( new ArgumentRule("characterHistory", CharacterHistory::getClassTypeSpec(), "The character history object from which we obtain the state indices.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        std::vector<TypeSpec> alphaTypes;
        alphaTypes.push_back( RealPos::getClassTypeSpec() );
        alphaTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "alpha" , alphaTypes, "The rate of attraction/selection (per state).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
        
        std::vector<TypeSpec> thetaTypes;
        thetaTypes.push_back( Real::getClassTypeSpec() );
        thetaTypes.push_back( ModelVector<Real>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "theta" , thetaTypes, "The optimum value (per state).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        
        std::vector<TypeSpec> sigmaTypes;
        sigmaTypes.push_back( RealPos::getClassTypeSpec() );
        sigmaTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "sigma" , sigmaTypes, "The rate of random drift (per state).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        
        std::vector<TypeSpec> rootStateTypes;
        rootStateTypes.push_back( Real::getClassTypeSpec() );
        rootStateTypes.push_back( ModelVector<Real>::getClassTypeSpec() );
        Real *defaultRootStates = new Real(0.0);
        dist_member_rules.push_back( new ArgumentRule( "rootStates" , rootStateTypes, "The vector of root states.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, defaultRootStates ) );
         
        std::vector<TypeSpec> observation_varianceTypes;
        observation_varianceTypes.push_back( RealPos::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "observation_variance" , observation_varianceTypes, "The observation variance of the trait (in units of trait^2).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        

        dist_member_rules.push_back( new ArgumentRule( "nSites"         ,  Natural::getClassTypeSpec(), "The number of sites which is used for the initialized (random draw) from this distribution.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(10) ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_PhyloOrnsteinUhlenbeckStateDependent::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_PhyloOrnsteinUhlenbeckStateDependent::printValue(std::ostream& o) const
{
    
    o << "PhyloOrnsteinUhlenbeckProcess(tree=";
    if ( character_history != NULL )
    {
        o << character_history->getName();
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
void Dist_PhyloOrnsteinUhlenbeckStateDependent::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "characterHistory" )
    {
        character_history = var;
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
    else if ( name == "observation_variance" )
    {
        observation_variance = var;
    }
    else if ( name == "rootStates" )
    {
        root_states = var;
    }
    else if ( name == "nSites" )
    {
        n_sites = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
    
}


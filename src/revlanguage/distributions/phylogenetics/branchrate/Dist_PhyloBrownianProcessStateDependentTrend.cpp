#include "Dist_PhyloBrownianProcessStateDependentTrend.h"

#include <stddef.h>
#include <ostream>

#include "PhyloBrownianProcessStateDependentTrend.h"
#include "RlTree.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "Natural.h"
#include "Real.h"
#include "OptionRule.h"
#include "RealPos.h"
#include "RlCharacterHistory.h"
#include "RlDistribution.h"
#include "RlString.h"
#include "StringUtilities.h"
#include "Tree.h"
#include "TypeSpec.h"

namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevLanguage;


Dist_PhyloBrownianProcessStateDependentTrend::Dist_PhyloBrownianProcessStateDependentTrend() : TypedDistribution< ContinuousCharacterData >()
{
    
}


Dist_PhyloBrownianProcessStateDependentTrend::~Dist_PhyloBrownianProcessStateDependentTrend()
{
    
}



Dist_PhyloBrownianProcessStateDependentTrend* Dist_PhyloBrownianProcessStateDependentTrend::clone( void ) const
{
    
    return new Dist_PhyloBrownianProcessStateDependentTrend(*this);
}


RevBayesCore::TypedDistribution< RevBayesCore::ContinuousCharacterData >* Dist_PhyloBrownianProcessStateDependentTrend::createDistribution( void ) const
{
    
    // get the parameters
    
    const CharacterHistory& rl_char_hist = static_cast<const RevLanguage::CharacterHistory&>( character_history->getRevObject() );
    RevBayesCore::TypedDagNode<RevBayesCore::CharacterHistoryDiscrete>* char_hist   =  rl_char_hist.getDagNode();

    RevBayesCore::PhyloBrownianProcessStateDependentTrend *dist = new RevBayesCore::PhyloBrownianProcessStateDependentTrend(char_hist);
    
    // set tau
    if ( tau->getRevObject().isType( ModelVector<Real>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* t = static_cast<const ModelVector<Real> &>( tau->getRevObject() ).getDagNode();
        dist->setTau( t );
    }
    else
    {
        RevBayesCore::TypedDagNode< double >* t = static_cast<const Real &>( tau->getRevObject() ).getDagNode();
        dist->setTau( t );
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

    RevBayesCore::TypedDagNode< double >* rs = static_cast<const Real &>( root_state->getRevObject() ).getDagNode();
    dist->setRootState( rs );
    
    return dist;
}



/* Get Rev type of object */
const std::string& Dist_PhyloBrownianProcessStateDependentTrend::getClassType(void)
{
    
    static std::string rev_type = "Dist_PhyloBrownianProcessStateDependentTrend";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_PhyloBrownianProcessStateDependentTrend::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_PhyloBrownianProcessStateDependentTrend::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "PhyloBMSDT" );
    a_names.push_back( "PhBMSDT" );

    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_PhyloBrownianProcessStateDependentTrend::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhyloBrownianProcessStateDependentTrend";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_PhyloBrownianProcessStateDependentTrend::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        dist_member_rules.push_back( new ArgumentRule("characterHistory", CharacterHistory::getClassTypeSpec(), "The character history object from which we obtain the state indices.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        std::vector<TypeSpec> tauTypes;
        tauTypes.push_back( Real::getClassTypeSpec() );
        tauTypes.push_back( ModelVector<Real>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "tau" , tauTypes, "The trend (per state).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Real(0.0) ) );
        
        std::vector<TypeSpec> sigmaTypes;
        sigmaTypes.push_back( RealPos::getClassTypeSpec() );
        sigmaTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "sigma" , sigmaTypes, "The rate of random drift (per state).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        
        std::vector<TypeSpec> rootStateTypes;
        rootStateTypes.push_back( Real::getClassTypeSpec() );
        Real *defaultrootState = new Real(0.0);
        dist_member_rules.push_back( new ArgumentRule( "rootState" , rootStateTypes, "The root state.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, defaultrootState ) );
       
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_PhyloBrownianProcessStateDependentTrend::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_PhyloBrownianProcessStateDependentTrend::printValue(std::ostream& o) const
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
    o << ", tau=";
    if ( tau != NULL )
    {
        o << tau->getName();
    }
    else
    {
        o << "?";
    }
    o << ", root_state=";
    if ( root_state != NULL )
    {
        o << root_state->getName();
    }
    else
    {
        o << "?";
    }
    
    o << ")";
    
}


/** Set a member variable */
void Dist_PhyloBrownianProcessStateDependentTrend::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "characterHistory" )
    {
        character_history = var;
    }
    else if ( name == "tau" )
    {
        tau = var;
    }
    else if ( name == "sigma" )
    {
        sigma = var;
    }
    else if ( name == "rootState" )
    {
        root_state = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
    
}


#include "Dist_PhyloBrownianProcessStateDependent.h"

#include <stddef.h>
#include <ostream>

#include "PhyloBrownianProcessStateDependent.h"
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


Dist_PhyloBrownianProcessStateDependent::Dist_PhyloBrownianProcessStateDependent() : TypedDistribution< ContinuousCharacterData >()
{
    
}


Dist_PhyloBrownianProcessStateDependent::~Dist_PhyloBrownianProcessStateDependent()
{
    
}



Dist_PhyloBrownianProcessStateDependent* Dist_PhyloBrownianProcessStateDependent::clone( void ) const
{
    
    return new Dist_PhyloBrownianProcessStateDependent(*this);
}


RevBayesCore::TypedDistribution< RevBayesCore::ContinuousCharacterData >* Dist_PhyloBrownianProcessStateDependent::createDistribution( void ) const
{
    
    // get the parameters
    size_t n = size_t( static_cast<const Natural &>( n_sites->getRevObject() ).getValue() );
    
    const CharacterHistory& rl_char_hist = static_cast<const RevLanguage::CharacterHistory&>( character_history->getRevObject() );
    RevBayesCore::TypedDagNode<RevBayesCore::CharacterHistoryDiscrete>* char_hist   =  rl_char_hist.getDagNode();
    
    RevBayesCore::PhyloBrownianProcessStateDependent *dist = new RevBayesCore::PhyloBrownianProcessStateDependent(char_hist, n);
    
    
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
    
    return dist;
}



/* Get Rev type of object */
const std::string& Dist_PhyloBrownianProcessStateDependent::getClassType(void)
{
    
    static std::string rev_type = "Dist_PhyloBrownianProcessStateDependent";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_PhyloBrownianProcessStateDependent::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_PhyloBrownianProcessStateDependent::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "PhyloBMSD" );
    a_names.push_back( "PhBMSD" );

    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_PhyloBrownianProcessStateDependent::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhyloBrownianProcessStateDependent";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_PhyloBrownianProcessStateDependent::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        dist_member_rules.push_back( new ArgumentRule("characterHistory", CharacterHistory::getClassTypeSpec(), "The character history object from which we obtain the state indices.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        std::vector<TypeSpec> sigmaTypes;
        sigmaTypes.push_back( RealPos::getClassTypeSpec() );
        sigmaTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "sigma" , sigmaTypes, "The rate of random drift (per state).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        
        dist_member_rules.push_back( new ArgumentRule( "nSites"         ,  Natural::getClassTypeSpec(), "The number of sites which is used for the initialized (random draw) from this distribution.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(10) ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_PhyloBrownianProcessStateDependent::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_PhyloBrownianProcessStateDependent::printValue(std::ostream& o) const
{
    
    o << "PhyloBrownianProcessStateDependent(tree=";
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
void Dist_PhyloBrownianProcessStateDependent::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "characterHistory" )
    {
        character_history = var;
    }
    else if ( name == "sigma" )
    {
        sigma = var;
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


#include "Dist_PhyloMultiSampleOrnsteinUhlenbeckStateDependent.h"

#include <stddef.h>
#include <ostream>

#include "PhyloMultiSampleOrnsteinUhlenbeckStateDependent.h"
#include "RlTree.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "Natural.h"
#include "Real.h"
#include "OptionRule.h"
#include "RealPos.h"
#include "RbBoolean.h"
#include "RlBoolean.h"
#include "RlCharacterHistory.h"
#include "RlDistribution.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "StringUtilities.h"
#include "Tree.h"
#include "TypeSpec.h"

namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevLanguage;


Dist_PhyloMultiSampleOrnsteinUhlenbeckStateDependent::Dist_PhyloMultiSampleOrnsteinUhlenbeckStateDependent() : TypedDistribution< ContinuousCharacterData >()
{

}


Dist_PhyloMultiSampleOrnsteinUhlenbeckStateDependent::~Dist_PhyloMultiSampleOrnsteinUhlenbeckStateDependent()
{

}



Dist_PhyloMultiSampleOrnsteinUhlenbeckStateDependent* Dist_PhyloMultiSampleOrnsteinUhlenbeckStateDependent::clone( void ) const
{

    return new Dist_PhyloMultiSampleOrnsteinUhlenbeckStateDependent(*this);
}


RevBayesCore::TypedDistribution< RevBayesCore::ContinuousCharacterData >* Dist_PhyloMultiSampleOrnsteinUhlenbeckStateDependent::createDistribution( void ) const
{

    // get the parameters
    size_t n = size_t( static_cast<const Natural &>( n_sites->getRevObject() ).getValue() );

    const CharacterHistory& rl_char_hist = static_cast<const RevLanguage::CharacterHistory&>( character_history->getRevObject() );
    RevBayesCore::TypedDagNode<RevBayesCore::CharacterHistoryDiscrete>* char_hist   =  rl_char_hist.getDagNode();
    size_t number_states = char_hist->getValue().getNumberOfStates();

   //    set the root treatment
    const std::string& rt = static_cast<const RlString &>( root_treatment->getRevObject() ).getValue();
    RevBayesCore::PhyloMultiSampleOrnsteinUhlenbeckStateDependent::ROOT_TREATMENT rtr;
    if (rt == "optimum")
    {
        rtr = RevBayesCore::PhyloMultiSampleOrnsteinUhlenbeckStateDependent::ROOT_TREATMENT::OPTIMUM;
    }
    else if (rt == "equilibrium")
    {
        rtr = RevBayesCore::PhyloMultiSampleOrnsteinUhlenbeckStateDependent::ROOT_TREATMENT::EQUILIBRIUM;
    }
    else if (rt == "parameter")
    {
        rtr = RevBayesCore::PhyloMultiSampleOrnsteinUhlenbeckStateDependent::ROOT_TREATMENT::PARAMETER;
    }
    else
    {
        throw RbException("argument rootTreatment must be one of \"optimum\", \"equilibrium\" or \"parameter\"");
    }

    const std::vector<RevBayesCore::Taxon> &ta = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();

    // set the within-species variances, or the variances of the within-species variances
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* var;
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* var_var;
    bool use_emp_var = static_cast<const RlBoolean &>( useEmpiricalSpeciesVariances->getRevObject() ).getDagNode()->getValue();
    if ( use_emp_var == true)
    {
         if ( within_species_variances->getRevObject() != RevNullObject::getInstance() )
         {
             throw RbException("To use empirical species variances, do not specify \"withinSpeciesVariances\"");
         }
         else if ( variances_of_within_species_variances->getRevObject() == RevNullObject::getInstance())
         {
             throw RbException("To use empirical species variances, you have to specify \"variancesOfWithinSpeciesVariances\"");
         }
         else
         {
             RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* var = static_cast<const ModelVector<RealPos> &>( within_species_variances->getRevObject() ).getDagNode();
         }
    }
    else        // useEmpiricalSpeciesVariances == FALSE
    {
        if ( within_species_variances->getRevObject() == RevNullObject::getInstance() )
        {
            throw RbException("If \"useEmpiricalSpeciesVariances\" is set to false, you have to specify \"withinSpeciesVariances\"");
        }
        else if ( variances_of_within_species_variances->getRevObject() != RevNullObject::getInstance())
        {
            throw RbException("If \"useEmpiricalSpeciesVariances\" is set to false, do not specify \"variancesOfWithinSpeciesVariances\"");
        }
        else
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* var_var = static_cast<const ModelVector<RealPos> &>( variances_of_within_species_variances->getRevObject() ).getDagNode();
        }
    }



    RevBayesCore::PhyloMultiSampleOrnsteinUhlenbeckStateDependent *dist = new RevBayesCore::PhyloMultiSampleOrnsteinUhlenbeckStateDependent(char_hist, n, rtr, var, var_var, ta);

    // set alpha
    if ( alpha->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* a = static_cast<const ModelVector<RealPos> &>( alpha->getRevObject() ).getDagNode();
        if ( a->getValue().size() == number_states )
        {
            dist->setAlpha( a );
        }
        else
        {
            throw RbException() << "The number of states (" << number_states << ") in the character history doesn't match the number of alpha parameters (" << a->getValue().size() << ")";
        }
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
        if ( t->getValue().size() == number_states )
        {
            dist->setTheta( t );
        }
        else
        {
            throw RbException() << "The number of states (" << number_states << ") in the character history doesn't match the number of theta parameters (" << t->getValue().size() << ")";
        }
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
        if ( s->getValue().size() == number_states )
        {
            dist->setSigma( s );
        }
        else
        {
            throw RbException() << "The number of states (" << number_states << ") in the character history doesn't match the number of sigma parameters (" << s->getValue().size() << ")";
        }
    }
    else
    {
        RevBayesCore::TypedDagNode< double >* s = static_cast<const RealPos &>( sigma->getRevObject() ).getDagNode();
        dist->setSigma( s );
    }

    // set the root value
    if ( rt == "optimum" || rt == "equilibrium" )
    {
         if ( root_value->getRevObject() != RevNullObject::getInstance() )
         {
             throw RbException("To use the root treatment \"optimum\" or \"equilibrium\", you should not specify the argument rootValue ");
         }
    }
    else if ( rt == "parameter" )
    {
        if ( root_value->getRevObject() != RevNullObject::getInstance() )
        {
            RevBayesCore::TypedDagNode< double >* rvl = static_cast<const Real &>( root_value->getRevObject() ).getDagNode();
            dist->setRootValue( rvl );
        }
        else
        {
            throw RbException("To use the root treatment \"parameter\", you need to specify the argument rootValue ");
        }
    }


    return dist;
}



/* Get Rev type of object */
const std::string& Dist_PhyloMultiSampleOrnsteinUhlenbeckStateDependent::getClassType(void)
{

    static std::string rev_type = "Dist_PhyloMultiSampleOrnsteinUhlenbeckStateDependent";

    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_PhyloMultiSampleOrnsteinUhlenbeckStateDependent::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_PhyloMultiSampleOrnsteinUhlenbeckStateDependent::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "PhyloMultiSampleOUSD" );
    a_names.push_back( "PhMultiSampleOUSD" );

    return a_names;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_PhyloMultiSampleOrnsteinUhlenbeckStateDependent::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhyloMultiSampleOrnsteinUhlenbeckStateDependent";

    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_PhyloMultiSampleOrnsteinUhlenbeckStateDependent::getParameterRules(void) const
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

        std::vector<TypeSpec> rootValueTypes;
        rootValueTypes.push_back( Real::getClassTypeSpec() );
        //Real *defaultRootValue = new Real(0.0);
        dist_member_rules.push_back( new ArgumentRule( "rootValue" , rootValueTypes, "The value of the continuous trait at root.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

        std::vector<std::string> rootTreatmentTypes;
        rootTreatmentTypes.push_back( "optimum" );
        rootTreatmentTypes.push_back( "equilibrium" );
        rootTreatmentTypes.push_back( "parameter" );
        dist_member_rules.push_back( new OptionRule ("rootTreatment", new RlString("optimum"), rootTreatmentTypes, "Whether the root value should be assumed to be equal to the optimum at the root (the default), assumed to be a random variable distributed according to the equilibrium state of the OU process, or whether to estimate the ancestral value as an independent parameter.") );
        // setting the default
        //RevBayesCore::PhyloMultiSampleOrnsteinUhlenbeckStateDependent::ROOT_TREATMENT rtr = RevBayesCore::PhyloMultiSampleOrnsteinUhlenbeckStateDependent::ROOT_TREATMENT::OPTIMUM;

        dist_member_rules.push_back( new ArgumentRule( "useEmpiricalSpeciesMeans",     RlBoolean::getClassTypeSpec(), "Should the species means assumed to be equal to the empirical species mean or should we estimate them?",                        ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        dist_member_rules.push_back( new ArgumentRule( "useEmpiricalSpeciesVariances", RlBoolean::getClassTypeSpec(), "If \"TRUE\", the within-species variance estimates are informed by the empirical within-species variances.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );

        dist_member_rules.push_back( new ArgumentRule( "withinSpeciesVariances" , ModelVector<RealPos>::getClassTypeSpec(), "The within-species variance for each species.",                         ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "variancesOfWithinSpeciesVariances" , ModelVector<RealPos>::getClassTypeSpec(), "The variances of within-species variance for each species.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "taxa"  , ModelVector<Taxon>::getClassTypeSpec(), "The vector of taxa which have species and individual names.",                              ArgumentRule::BY_VALUE,              ArgumentRule::ANY ) );

        dist_member_rules.push_back( new ArgumentRule( "nSites",  Natural::getClassTypeSpec(), "The number of sites which is used for the initialized (random draw) from this distribution.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(1) ) );

        rules_set = true;
    }

    return dist_member_rules;
}


const TypeSpec& Dist_PhyloMultiSampleOrnsteinUhlenbeckStateDependent::getTypeSpec( void ) const
{

    static TypeSpec ts = getClassTypeSpec();

    return ts;
}


/** Print value for user */
void Dist_PhyloMultiSampleOrnsteinUhlenbeckStateDependent::printValue(std::ostream& o) const
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
    o << ", alpha=";
    if ( alpha != NULL )
    {
        o << alpha->getName();
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
    o << ", rootValue=";
    if ( root_value != NULL )
    {
        o << root_value->getName();
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
void Dist_PhyloMultiSampleOrnsteinUhlenbeckStateDependent::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
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
    else if ( name == "rootValue" )
    {
        root_value = var;
    }
    else if ( name == "nSites" )
    {
        n_sites = var;
    }
    else if ( name == "rootTreatment" )
    {
        root_treatment = var;
    }
    else if ( name == "withinSpeciesVariances" )
    {
        within_species_variances = var;
    }
    else if ( name == "variancesOfWithinSpeciesVariances" )
    {
        variances_of_within_species_variances = var;
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

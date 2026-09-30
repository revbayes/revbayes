#include "Dist_PhyloMultiSampleOrnsteinUhlenbeckStateDependent.h"

#include <cstddef>
#include <stddef.h>
#include <ostream>

#include "PhyloMultiSampleOrnsteinUhlenbeckStateDependent.h"
#include "RlTree.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "Natural.h"
#include "MatrixReal.h"
#include "Real.h"
#include "OptionRule.h"
#include "RealPos.h"
#include "RlCharacterHistory.h"
#include "RlDistribution.h"
#include "RlMatrixReal.h"
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
    if ( n != 1 )
    {
         throw RbException("The state-dependent OU model only supports univariate continuous character. To infer multiple univariate continuous characters under the same character history, please create multiple OUSD distributions.");
    }
    
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

    // set the treatment for variance of species means for species with one sample only
    // const std::string& sst = static_cast<const RlString &>( single_sample_treatment->getRevObject() ).getValue();
    // RevBayesCore::PhyloMultiSampleOrnsteinUhlenbeckStateDependent::SINGLE_SAMPLE_TREATMENT sstr;
    // if (sst == "mean")
    // {
    //     sstr = RevBayesCore::PhyloMultiSampleOrnsteinUhlenbeckStateDependent::SINGLE_SAMPLE_TREATMENT::MEAN;
    // }
    // else if (sst == "median")
    // {
    //     sstr = RevBayesCore::PhyloMultiSampleOrnsteinUhlenbeckStateDependent::SINGLE_SAMPLE_TREATMENT::MEDIAN;
    // }
    // else if (sst == "as_is")
    // {
    //     sstr = RevBayesCore::PhyloMultiSampleOrnsteinUhlenbeckStateDependent::SINGLE_SAMPLE_TREATMENT::AS_IS;
    // }
    // else
    // {
    //     throw RbException("argument singleSampleTreatment must be one of \"mean\", \"median\" or \"as_is\"");
    // }

    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* wsv = static_cast<const ModelVector<RealPos> &>( species_var->getRevObject() ).getDagNode();

    const std::vector<RevBayesCore::Taxon> &ta  = static_cast<const ModelVector<Taxon> &>( taxon_map->getRevObject() ).getValue();


    // RevBayesCore::PhyloMultiSampleOrnsteinUhlenbeckStateDependent *dist = new RevBayesCore::PhyloMultiSampleOrnsteinUhlenbeckStateDependent(char_hist, n, rtr, sstr, ta, wsv);
    RevBayesCore::PhyloMultiSampleOrnsteinUhlenbeckStateDependent *dist = new RevBayesCore::PhyloMultiSampleOrnsteinUhlenbeckStateDependent(char_hist, n, rtr, ta, wsv);

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

    // if ( species_var->getRevObject() != RevNullObject::getInstance() )
    // {
//         if ( num_samples_per_species->getRevObject() == RevNullObject::getInstance() )
//         {
//             throw RbException() << "Please also provide the number of samples per species if you want to include uncertainty at the tips.";
// 
//         }

        // RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double>>* sp_var  = static_cast<const ModelVector<RealPos>&>( species_var->getRevObject() ).getDagNode();
// 
//         if (sp_var->getValue().size() != n)
//         {
//             throw RbException()<< "The number of sites (" << n << ") specified doesn't match the size of the within-species variance matrix (" << sp_var->getValue().size() << ")";
//         }
//         else
//         {
        // dist->setWithinSpeciesVariance( sp_var );
        // }

    // }

//     if ( num_samples_per_species->getRevObject() != RevNullObject::getInstance() )
//     {
//         if ( species_var->getRevObject() == RevNullObject::getInstance() )
//         {
//             throw RbException() << "Please also provide the number of samples per species if you want to include uncertainty at the tips.";
// 
//         }
//         else
//         {
//             RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double>>* n_samples  = static_cast<const ModelVector<RealPos>&>( num_samples_per_species->getRevObject() ).getDagNode();

            // if (n_samples->getValue().size() != n)
            // {
            //     throw RbException()<< "The number of sites (" << n << ") specified doesn't match the size of the number-of-samples-per-species matrix (" << n_samples->getValue().size() << ")";
            // }
            // else
            // {
//             dist->setNumberOfSamplesPerSpecies( n_samples );
//             }
//         }
// 
//     }

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
    a_names.push_back( "PhyloMSOUSD" );
    a_names.push_back( "PhyloMSSDOU" );
    a_names.push_back( "PhMSOUSD" );

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


        dist_member_rules.push_back( new ArgumentRule( "withinSpeciesVariance" , ModelVector<RealPos>::getClassTypeSpec(), "The within-species variance for each species.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

        dist_member_rules.push_back( new ArgumentRule( "numberOfSamplesPerSpecies" , ModelVector<RealPos>::getClassTypeSpec(), "The number of samples for each species.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

        // std::vector<std::string> singleSampleTreatmentTypes;
        // singleSampleTreatmentTypes.push_back( "mean" );
        // singleSampleTreatmentTypes.push_back( "median" );
        // singleSampleTreatmentTypes.push_back( "as_is" );
        // dist_member_rules.push_back( new OptionRule ("singleSampleTreatment", new RlString("mean"), singleSampleTreatmentTypes, "What to be input as the variance of species mean at the tip is the species contains one sample only. Options \"mean\" and \"median\" calculate the mean/median of the variance of species mean for species with multiple sample. Option \"as_is\" uses the value provided in the vector of \"withinSpeciesVariance\" directly.") );

        dist_member_rules.push_back( new ArgumentRule( "nSites",  Natural::getClassTypeSpec(), "The number of continuous character.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(1) ) );

        dist_member_rules.push_back( new ArgumentRule( "taxa"  , ModelVector<Taxon>::getClassTypeSpec(), "The vector of taxa which have species and individual names.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
                
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
    if ( species_var != NULL )
    {
        o << species_var->getName();
    }
    else
    {
        o << "?";
    }
    // if ( num_samples_per_species != NULL )
    // {
    //     o << num_samples_per_species->getName();
    // }
    // else
    // {
    //     o << "?";
    // }
    // if ( single_sample_treatment != NULL )
    // {
    //     o << single_sample_treatment->getName();
    // }
    // else
    // {
    //     o << "?";
    // }
    if ( taxon_map != NULL )
    {
        o << taxon_map->getName();
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
    else if ( name == "withinSpeciesVariance" )
    {
        species_var = var;
    }
    //else if ( name == "numberOfSamplesPerSpecies" )
    //{
    //    num_samples_per_species = var;
    //}
    // else if ( name == "singleSampleTreatment" )
    // {
    //     single_sample_treatment = var;
    // }
    else if ( name == "taxa" )
    {
        taxon_map = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }

}

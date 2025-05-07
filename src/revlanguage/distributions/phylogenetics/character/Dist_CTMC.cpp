#include "Dist_CTMC.h"

#include <cstddef>
#include <ostream>

#include "CTMCProcess.h"

#include "AminoAcidState.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "BinaryState.h"
#include "CodonState.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DnaState.h"
#include "MatrixReal.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "Natural.h"
#include "NaturalNumbersState.h"
#include "OptionRule.h"
#include "PoMoState.h"
#include "RbBoolean.h"
#include "RbException.h"
#include "RbVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RevVariable.h"
#include "RlDistributionMemberFunction.h"
#include "RlMatrixReal.h"
#include "RlRateGenerator.h"
#include "RlString.h"
#include "RlSimplex.h"
#include "RnaState.h"
#include "StandardState.h"

using namespace RevLanguage;

Dist_CTMC::Dist_CTMC() : TypedDistribution< AbstractDiscreteTaxonData >()
{

}


Dist_CTMC::~Dist_CTMC()
{

}



Dist_CTMC* Dist_CTMC::clone( void ) const
{

    return new Dist_CTMC(*this);
}


RevBayesCore::TypedDistribution< RevBayesCore::AbstractDiscreteTaxonData >* Dist_CTMC::createDistribution( void ) const
{

    // get the parameters
    size_t this_num_sites = size_t( static_cast<const Natural &>( nSites->getRevObject() ).getValue() );
    const std::string& dt = static_cast<const RlString &>( type->getRevObject() ).getValue();

    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* site_rates_node = NULL;
    if ( site_rates != NULL && site_rates->getRevObject() != RevNullObject::getInstance() )
    {
        site_rates_node = static_cast<const ModelVector<RealPos> &>( site_rates->getRevObject() ).getDagNode();
    }
    
    const RevBayesCore::TypedDagNode< RevBayesCore::Simplex > *site_rates_probs_node = NULL;
    if ( site_rates_probs->getRevObject() != RevNullObject::getInstance() )
    {
        site_rates_probs_node = static_cast<const Simplex &>( site_rates_probs->getRevObject() ).getDagNode();
    }
    if (site_rates_probs_node != NULL)
    {
        if (site_rates_node == NULL)
        {
            throw RbException( "Provided site rates probs but not using site rates." );
        }
        else if (site_rates_probs_node->getValue().size() != site_rates_node->getValue().size())
        {
            throw RbException( "The number of site rates probs does not match the number of site rates." );
        }
    }

    RevBayesCore::TypedDistribution< RevBayesCore::AbstractDiscreteTaxonData > *d = NULL;
    const RevBayesCore::TypedDagNode< RevBayesCore::Simplex > *rf = NULL;
    if ( root_frequencies->getRevObject() != RevNullObject::getInstance() )
    {
        rf = static_cast<const Simplex &>( root_frequencies->getRevObject() ).getDagNode();
    }

    const RevBayesCore::TypedDagNode< RevBayesCore::Simplex > *sp = NULL;

    if ( site_matrices->getRevObject().isType( Simplex::getClassTypeSpec() ) )
    {

        sp = static_cast<const Simplex &>( site_matrices->getRevObject() ).getDagNode();

        if ( q->getRevObject().isType( ModelVector<RateGenerator>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( q->getRevObject() ).getDagNode();
            // sanity check
            if ( sp->getValue().size() != rm->getValue().size())
            {
                throw RbException( "The number of substitution matrices does not match the number of matrix probabilities" );
            }
        }
        else
        {
            throw RbException( "The number of substitution matrices does not match the number of matrix mixture probabilities" );
        }
    }

    if ( dt == "DNA" )
    {
        RevBayesCore::CTMCProcess<RevBayesCore::DnaState> *dist = new RevBayesCore::CTMCProcess<RevBayesCore::DnaState>(4,this_num_sites);

        // set the root frequencies (by default these are NULL so this is OK)
        dist->setRootFrequencies( rf );

        dist->setSiteMatricesProbs(sp);

        // set the rate matrix
        if ( q->getRevObject().isType( ModelVector<RateGenerator>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( q->getRevObject() ).getDagNode();
            
            dist->setRateMatrix( rm );
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rm = static_cast<const RateGenerator &>( q->getRevObject() ).getDagNode();
            dist->setRateMatrix( rm );
        }

        if ( site_rates_node != NULL && site_rates_node->getValue().size() > 0 )
        {
            dist->setSiteRates( site_rates_node );
        }
        
        if ( site_rates_probs_node != NULL && site_rates_probs_node->getValue().size() > 0 )
        {
            dist->setSiteRatesProbs( site_rates_probs_node );
        }
        

        d = dist;
    }
    else if ( dt == "RNA" )
    {
        RevBayesCore::CTMCProcess<RevBayesCore::RnaState> *dist = new RevBayesCore::CTMCProcess<RevBayesCore::RnaState>(4,this_num_sites);

        // set the root frequencies (by default these are NULL so this is OK)
        dist->setRootFrequencies( rf );

        dist->setSiteMatricesProbs(sp);

        // set the rate matrix
        if ( q->getRevObject().isType( ModelVector<RateGenerator>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( q->getRevObject() ).getDagNode();
            
            dist->setRateMatrix( rm );
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rm = static_cast<const RateGenerator &>( q->getRevObject() ).getDagNode();
            dist->setRateMatrix( rm );
        }

        if ( site_rates_node != NULL && site_rates_node->getValue().size() > 0 )
        {
            dist->setSiteRates( site_rates_node );
        }
        
        if ( site_rates_probs_node != NULL && site_rates_probs_node->getValue().size() > 0 )
        {
            dist->setSiteRatesProbs( site_rates_probs_node );
        }
        

        d = dist;
    }
    else if ( dt == "AA" || dt == "Protein" )
    {
        RevBayesCore::CTMCProcess<RevBayesCore::AminoAcidState> *dist = new RevBayesCore::CTMCProcess<RevBayesCore::AminoAcidState>(20,this_num_sites);

        // set the root frequencies (by default these are NULL so this is OK)
        dist->setRootFrequencies( rf );

        dist->setSiteMatricesProbs(sp);

        // set the rate matrix
        if ( q->getRevObject().isType( ModelVector<RateGenerator>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( q->getRevObject() ).getDagNode();
            
            dist->setRateMatrix( rm );
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rm = static_cast<const RateGenerator &>( q->getRevObject() ).getDagNode();
            dist->setRateMatrix( rm );
        }

        if ( site_rates_node != NULL && site_rates_node->getValue().size() > 0 )
        {
            dist->setSiteRates( site_rates_node );
        }
        
        if ( site_rates_probs_node != NULL && site_rates_probs_node->getValue().size() > 0 )
        {
            dist->setSiteRatesProbs( site_rates_probs_node );
        }
        

        d = dist;
    }
    else if ( dt == "Codon" )
    {
        
        RevBayesCore::CTMCProcess<RevBayesCore::CodonState> *dist = new RevBayesCore::CTMCProcess<RevBayesCore::CodonState>(61,this_num_sites);

        // set the root frequencies (by default these are NULL so this is OK)
        dist->setRootFrequencies( rf );

        dist->setSiteMatricesProbs(sp);

        // set the rate matrix
        if ( q->getRevObject().isType( ModelVector<RateGenerator>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( q->getRevObject() ).getDagNode();
            
            dist->setRateMatrix( rm );
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rm = static_cast<const RateGenerator &>( q->getRevObject() ).getDagNode();
            dist->setRateMatrix( rm );
        }

        if ( site_rates_node != NULL && site_rates_node->getValue().size() > 0 )
        {
            dist->setSiteRates( site_rates_node );
        }
        
        if ( site_rates_probs_node != NULL && site_rates_probs_node->getValue().size() > 0 )
        {
            dist->setSiteRatesProbs( site_rates_probs_node );
        }
        

        d = dist;
    }
    else if ( dt == "PoMo" )
    {

        // we get the number of states from the rate matrix (we don't know, because PoMo is flexible about its rates)
        // set the rate matrix
        size_t nChars = 1;
        if ( q->getRevObject().isType( ModelVector<RateGenerator>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( q->getRevObject() ).getDagNode();
            nChars = rm->getValue()[0].getNumberOfStates();
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rm = static_cast<const RateGenerator &>( q->getRevObject() ).getDagNode();
            nChars = rm->getValue().getNumberOfStates();
        }
        
        RevBayesCore::CTMCProcess<RevBayesCore::PoMoState> *dist = new RevBayesCore::CTMCProcess<RevBayesCore::PoMoState>(nChars,this_num_sites);

        // set the root frequencies (by default these are NULL so this is OK)
        dist->setRootFrequencies( rf );

        dist->setSiteMatricesProbs(sp);

        // set the rate matrix
        if ( q->getRevObject().isType( ModelVector<RateGenerator>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( q->getRevObject() ).getDagNode();
            
            dist->setRateMatrix( rm );
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rm = static_cast<const RateGenerator &>( q->getRevObject() ).getDagNode();
            dist->setRateMatrix( rm );
        }

        if ( site_rates_node != NULL && site_rates_node->getValue().size() > 0 )
        {
            dist->setSiteRates( site_rates_node );
        }
        
        if ( site_rates_probs_node != NULL && site_rates_probs_node->getValue().size() > 0 )
        {
            dist->setSiteRatesProbs( site_rates_probs_node );
        }
        

        d = dist;
    }
    else if ( dt == "Standard" )
    {
        // we get the number of states from the rates matrix
        // set the rate matrix
        size_t nChars = 1;
        if ( q->getRevObject().isType( ModelVector<RateGenerator>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( q->getRevObject() ).getDagNode();
            nChars = rm->getValue()[0].getNumberOfStates();
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rm = static_cast<const RateGenerator &>( q->getRevObject() ).getDagNode();
            nChars = rm->getValue().getNumberOfStates();
        }

        
        RevBayesCore::CTMCProcess<RevBayesCore::StandardState> *dist = new RevBayesCore::CTMCProcess<RevBayesCore::StandardState>(nChars,this_num_sites);

        // set the root frequencies (by default these are NULL so this is OK)
        dist->setRootFrequencies( rf );

        dist->setSiteMatricesProbs(sp);

        // set the rate matrix
        if ( q->getRevObject().isType( ModelVector<RateGenerator>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( q->getRevObject() ).getDagNode();
            
            dist->setRateMatrix( rm );
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rm = static_cast<const RateGenerator &>( q->getRevObject() ).getDagNode();
            dist->setRateMatrix( rm );
        }

        if ( site_rates_node != NULL && site_rates_node->getValue().size() > 0 )
        {
            dist->setSiteRates( site_rates_node );
        }
        
        if ( site_rates_probs_node != NULL && site_rates_probs_node->getValue().size() > 0 )
        {
            dist->setSiteRatesProbs( site_rates_probs_node );
        }
        

        d = dist;
    }
    else if ( dt == "NaturalNumbers" )
    {
        // we get the number of states from the rates matrix
        size_t n_chars = 1;
        if ( q->getRevObject().isType( ModelVector<RateGenerator>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( q->getRevObject() ).getDagNode();
            n_chars = rm->getValue()[0].getNumberOfStates();
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rm = static_cast<const RateGenerator &>( q->getRevObject() ).getDagNode();
            n_chars = rm->getValue().getNumberOfStates();
        }

        
        RevBayesCore::CTMCProcess<RevBayesCore::NaturalNumbersState> *dist = new RevBayesCore::CTMCProcess<RevBayesCore::NaturalNumbersState>(n_chars,this_num_sites);

        // set the root frequencies (by default these are NULL so this is OK)
        dist->setRootFrequencies( rf );

        dist->setSiteMatricesProbs(sp);

        // set the rate matrix
        if ( q->getRevObject().isType( ModelVector<RateGenerator>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( q->getRevObject() ).getDagNode();
            
            dist->setRateMatrix( rm );
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rm = static_cast<const RateGenerator &>( q->getRevObject() ).getDagNode();
            dist->setRateMatrix( rm );
        }

        if ( site_rates_node != NULL && site_rates_node->getValue().size() > 0 )
        {
            dist->setSiteRates( site_rates_node );
        }
        
        if ( site_rates_probs_node != NULL && site_rates_probs_node->getValue().size() > 0 )
        {
            dist->setSiteRatesProbs( site_rates_probs_node );
        }
        

        d = dist;
    }
    else if ( dt == "Binary" || dt == "Restriction" )
    {
        // we get the number of states from the rates matrix
        // set the rate matrix
        size_t nChars = 1;

        if ( q->getRevObject().isType( ModelVector<RateGenerator>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( q->getRevObject() ).getDagNode();
            nChars = rm->getValue()[0].getNumberOfStates();
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rm = static_cast<const RateGenerator &>( q->getRevObject() ).getDagNode();
            nChars = rm->getValue().getNumberOfStates();
        }

        // sanity check
        if ( nChars != 2 )
        {
            throw RbException( "Only binary characters allowed for type=Binary/Restriction" );
        }

        
        RevBayesCore::CTMCProcess<RevBayesCore::BinaryState> *dist = new RevBayesCore::CTMCProcess<RevBayesCore::BinaryState>(2,this_num_sites);

        // set the root frequencies (by default these are NULL so this is OK)
        dist->setRootFrequencies( rf );

        dist->setSiteMatricesProbs(sp);

        // set the rate matrix
        if ( q->getRevObject().isType( ModelVector<RateGenerator>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( q->getRevObject() ).getDagNode();
            
            dist->setRateMatrix( rm );
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rm = static_cast<const RateGenerator &>( q->getRevObject() ).getDagNode();
            dist->setRateMatrix( rm );
        }

        if ( site_rates_node != NULL && site_rates_node->getValue().size() > 0 )
        {
            dist->setSiteRates( site_rates_node );
        }
        
        if ( site_rates_probs_node != NULL && site_rates_probs_node->getValue().size() > 0 )
        {
            dist->setSiteRatesProbs( site_rates_probs_node );
        }
        

        d = dist;
    }


    return d;
}



/* Get Rev type of object */
const std::string& Dist_CTMC::getClassType(void)
{

    static std::string rev_type = "Dist_CTMC";

    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_CTMC::getClassTypeSpec(void)
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
std::string Dist_CTMC::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "CTMC";
    
    return d_name;
}


MethodTable Dist_CTMC::getDistributionMethods( void ) const
{
    
    MethodTable methods = TypedDistribution<AbstractDiscreteTaxonData>::getDistributionMethods();
    
    // member functions
    ArgumentRules* siteLikelihoodsArgRules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_CTMC, ModelVector<Real> >( "siteLikelihoods", variable, siteLikelihoodsArgRules, true ) );
    
    // f(Xh,rk|theta); returning a matrixreal where each row corresponds a site and each column corresponds an array (length = (num_site_rates + pInv != 0)) of rate categories [site_rate_index + pInv != 0] with the first element under inv if pInv > 0
    ArgumentRules* siteRateLikelihoodsArgRules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_CTMC, MatrixReal >( "siteRateLikelihoods", variable, siteRateLikelihoodsArgRules, true ) );
    
    // f(Xh,mk|theta); returning a matrixreal where each row corresponds a site and each column corresponds an array (length = ((num_site_rates + pInv != 0) * num_site_matrices)) of mixture categories [(site_rate_index + pInv != 0) * num_site_matrices + site_matrix_index] with the first num_site_matrices elements under inv if pInv > 0
    ArgumentRules* siteMixtureLikelihoodsArgRules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_CTMC, MatrixReal >( "siteMixtureLikelihoods", variable, siteMixtureLikelihoodsArgRules, true ) );
    
    ArgumentRules* siteRatesArgRules = new ArgumentRules();
    
    std::vector<std::string> optionsMethod;
    optionsMethod.push_back( "sampling" );
    optionsMethod.push_back( "weightedAverage" );
    siteRatesArgRules->push_back( new OptionRule( "estimateMethod", new RlString("sampling"), optionsMethod, "The method used to estimate the site specific rate." ) );
    
    methods.addFunction( new DistributionMemberFunction<Dist_CTMC, ModelVector<RealPos> >( "siteRates", variable, siteRatesArgRules, true ) );
    
    return methods;
}


/** Return member rules (no members) */
const MemberRules& Dist_CTMC::getParameterRules(void) const
{

    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {

        std::vector<TypeSpec> rateMatrixTypes;
        rateMatrixTypes.push_back( RateGenerator::getClassTypeSpec() );
        rateMatrixTypes.push_back( ModelVector<RateGenerator>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "Q", rateMatrixTypes, "The global or site-mixture rate matrices.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        // optional argument for the root frequencies
        dist_member_rules.push_back( new ArgumentRule( "rootFrequencies", Simplex::getClassTypeSpec(), "The root specific frequencies of the characters, if applicable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

        ModelVector<RealPos> *defaultSiteRates = new ModelVector<RealPos>();
        std::vector<TypeSpec> matrix_probs_types;
        matrix_probs_types.push_back(Simplex::getClassTypeSpec());
        matrix_probs_types.push_back(RlBoolean::getClassTypeSpec());

        dist_member_rules.push_back( new ArgumentRule( "siteMatrices", matrix_probs_types, "Simplex of site matrix mixture probabilities. Treats Q as vector of site mixture categories instead of branch-specific matrices.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "siteRates", ModelVector<RealPos>::getClassTypeSpec(), "The rate categories for the sites.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, defaultSiteRates ) );
        dist_member_rules.push_back( new ArgumentRule( "siteRatesProbs", Simplex::getClassTypeSpec(), "The probability weights of rate categories for the sites.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

        dist_member_rules.push_back( new ArgumentRule( "nSites", Natural::getClassTypeSpec(), "The number of sites, used for simulation.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural() ) );

        std::vector<std::string> options;
        options.push_back( "DNA" );
        options.push_back( "RNA" );
        options.push_back( "AA" );
        options.push_back( "Codon" );
        options.push_back( "PoMo" );
        options.push_back( "Protein" );
        options.push_back( "Standard" );
        options.push_back( "NaturalNumbers" );
        options.push_back( "Binary" );
        options.push_back( "Restriction" );
        dist_member_rules.push_back( new OptionRule( "type", new RlString("DNA"), options, "The data type, used for simulation and initialization." ) );

        rules_set = true;
    }

    return dist_member_rules;
}


const TypeSpec& Dist_CTMC::getTypeSpec( void ) const
{

    static TypeSpec ts = getClassTypeSpec();

    return ts;
}


/** Print value for user */
void Dist_CTMC::printValue(std::ostream& o) const
{

    o << "CTMC(Q=";
    if ( q != NULL )
    {
        o << q->getName();
    }
    else
    {
        o << "?";
    }
    o << ", site_rates=";
    if ( site_rates != NULL )
    {
        o << site_rates->getName();
    }
    else
    {
        o << "?";
    }
    o << ", site_rates_probs=";
    if ( site_rates_probs != NULL )
    {
        o << site_rates_probs->getName();
    }
    else
    {
        o << "?";
    }
    o << ", site_matrices=";
    if ( site_matrices != NULL )
    {
        o << site_matrices->getName();
    }
    else
    {
        o << "?";
    }
    o << ", nSites=";
    if ( nSites != NULL )
    {
        o << nSites->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";

}


/** Set a member variable */
void Dist_CTMC::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "Q" )
    {
        q = var;
    }
    else if ( name == "rootFrequencies" )
    {
        root_frequencies = var;
    }
    else if ( name == "siteRates" )
    {
        site_rates = var;
    }
    else if ( name == "siteRatesProbs" )
    {
        site_rates_probs = var;
    }
    else if ( name == "siteMatrices" )
    {
        site_matrices = var;
    }
    else if ( name == "nSites" )
    {
        nSites = var;
    }
    else if ( name == "type" )
    {
        type = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }

}


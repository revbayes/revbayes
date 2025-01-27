#include "Dist_phyloCTMC.h"

#include <cstddef>
#include <ostream>

#include "RlDistributionMemberFunction.h"
#include "PhyloCTMCSiteHomogeneous.h"
#include "PhyloCTMCSiteHomogeneousBinary.h"
#include "PhyloCTMCSiteHomogeneousNucleotide.h"
#include "OptionRule.h"
#include "Probability.h"
#include "RevNullObject.h"
#include "RlBoolean.h"
#include "RlMatrixReal.h"
#include "RlRateGenerator.h"
#include "RlSiteMixtureModel.h"
#include "RlString.h"
#include "RlTree.h"
#include "StandardState.h"
#include "RlSimplex.h"
#include "PoMoState.h"
#include "NaturalNumbersState.h"
#include "AbstractPhyloCTMCSiteHomogeneous.h"
#include "AminoAcidState.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "CodonState.h"
#include "DoubletState.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DiscreteTaxonData.h"
#include "DistributionMemberFunction.h"
#include "DnaState.h"
#include "HomologousDiscreteCharacterData.h"
#include "IndirectReferenceFunction.h"
#include "MatrixReal.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "Natural.h"
#include "PhyloCTMCSiteHomogeneousConditional.h"
#include "RateGenerator.h"
#include "RbBoolean.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Real.h"
#include "RealPos.h"
#include "RlConstantNode.h"
#include "RlDistribution.h"
#include "RnaState.h"
#include "Simplex.h"
#include "StringUtilities.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

Dist_phyloCTMC::Dist_phyloCTMC() : TypedDistribution< AbstractHomologousDiscreteCharacterData >()
{

}


Dist_phyloCTMC::~Dist_phyloCTMC()
{

}



Dist_phyloCTMC* Dist_phyloCTMC::clone( void ) const
{

    return new Dist_phyloCTMC(*this);
}

int Dist_phyloCTMC::computeNumberOfStates() const
{
    // we get the number of states from the rate matrix (we don't know, because PoMo is flexible about its rates)
    if ( q->getRevObject().isType( ModelVector<RateGenerator>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( q->getRevObject() ).getDagNode();
        return rm->getValue()[0].getNumberOfStates();
    }
    else if (q->getRevObject().isType( SiteMixtureModel::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::SiteMixtureModel >* mm = static_cast<const SiteMixtureModel &>( q->getRevObject() ).getDagNode();
        return mm->getValue().getNumberOfStates();
    }
    else
    {
        RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rm = static_cast<const RateGenerator &>( q->getRevObject() ).getDagNode();
        return rm->getValue().getNumberOfStates();
    }

    return 1;
}

template <typename Dist>
RevBayesCore::TypedDistribution< RevBayesCore::AbstractHomologousDiscreteCharacterData >* Dist_phyloCTMC::setDistParameters(Dist* dist) const
{
    const std::string& dt = static_cast<const RlString &>( type->getRevObject() ).getValue();

    // get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const Tree &>( tree->getRevObject() ).getDagNode();
    size_t nNodes = tau->getValue().getNumberOfNodes();

    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* site_ratesNode = NULL;
    if ( site_rates != NULL && site_rates->getRevObject() != RevNullObject::getInstance() )
    {
        site_ratesNode = static_cast<const ModelVector<RealPos> &>( site_rates->getRevObject() ).getDagNode();
    }
    
    const RevBayesCore::TypedDagNode< RevBayesCore::Simplex > *site_rates_probsNode = NULL;
    if ( site_rates_probs->getRevObject() != RevNullObject::getInstance() )
    {
        site_rates_probsNode = static_cast<const Simplex &>( site_rates_probs->getRevObject() ).getDagNode();
    }
    if (site_rates_probsNode != NULL)
    {
        if (site_ratesNode == NULL)
        {
            throw RbException( "Provided site rates probs but not using site rates." );
        }
        else if (site_rates_probsNode->getValue().size() != site_ratesNode->getValue().size())
        {
            throw RbException( "The number of site rates probs does not match the number of site rates." );
        }
    }
    
    RevBayesCore::TypedDagNode< double >* p_invNode = NULL;
    if ( p_inv != NULL && p_inv->getRevObject() != RevNullObject::getInstance() )
    {
        p_invNode = static_cast<const Probability &>( p_inv->getRevObject() ).getDagNode();
    }

    RevBayesCore::TypedDistribution< RevBayesCore::AbstractHomologousDiscreteCharacterData > *d = NULL;
    const RevBayesCore::TypedDagNode< RevBayesCore::Simplex > *rf = NULL;
    if ( root_frequencies->getRevObject() != RevNullObject::getInstance() )
    {
        rf = static_cast<const Simplex &>( root_frequencies->getRevObject() ).getDagNode();
    }

    const RevBayesCore::TypedDagNode< RevBayesCore::Simplex > *sp = NULL;

    bool use_site_matrices = false;

    if ( site_matrices->getRevObject().isType( Simplex::getClassTypeSpec() ) )
    {
        use_site_matrices = true;

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
    else if ( site_matrices->getRevObject().isType( RlBoolean::getClassTypeSpec() ) )
    {
        use_site_matrices = static_cast<const RlBoolean &>( site_matrices->getRevObject() ).getDagNode()->getValue();
    }

    // set the root frequencies (by default these are NULL so this is OK)
    dist->setRootFrequencies( rf );

    // set the probability for invariant site (by default this p_inv=0.0)
    dist->setPInv( p_invNode );

    // set the clock rates
    if ( rate->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    {
	RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* clockRates = static_cast<const ModelVector<RealPos> &>( rate->getRevObject() ).getDagNode();

	// sanity check
	if ( (nNodes-1) != clockRates->getValue().size() )
	{
	    throw RbException( "The number of clock rates does not match the number of branches" );
	}

	dist->setClockRate( clockRates );
    }
    else
    {
	RevBayesCore::TypedDagNode<double>* clockRate = static_cast<const RealPos &>( rate->getRevObject() ).getDagNode();
	dist->setClockRate( clockRate );
    }
    dist->setUseSiteMatrices(use_site_matrices, sp);

    // set the rate matrix
    if ( q->getRevObject().isType( ModelVector<RateGenerator>::getClassTypeSpec() ) )
    {
	RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( q->getRevObject() ).getDagNode();
	if (use_site_matrices == false)
	{
	    // sanity check
	    if ( (nNodes-1) != rm->getValue().size())
	    {
		throw RbException( "The number of substitution matrices does not match the number of branches" );
	    }

	    // sanity check
	    if ( root_frequencies == NULL || root_frequencies->getRevObject() == RevNullObject::getInstance() )
	    {
		throw RbException( "If you provide branch-heterogeneous substitution matrices, then you also need to provide root frequencies." );
	    }
	}
	dist->setRateMatrix( rm );
    }
    else if (q->getRevObject().isType( SiteMixtureModel::getClassTypeSpec() ) )
    {
	RevBayesCore::TypedDagNode< RevBayesCore::SiteMixtureModel >* mm = static_cast<const SiteMixtureModel &>( q->getRevObject() ).getDagNode();
	dist->setMixtureModel( mm );

	if ( root_frequencies and root_frequencies->getRevObject() != RevNullObject::getInstance() )
	    throw RbException()<<"dnPhyloCTMC: can't specify 'rootFrequencies' if Q if a SiteMixtureModel";

	if ( site_rates and site_rates->getRevObject() != RevNullObject::getInstance() )
	    throw RbException()<<"dnPhyloCTMC: can't specify 'siteRates' if Q if a SiteMixtureModel";

	if ( p_inv and p_inv->getRevObject() != RevNullObject::getInstance() )
	    throw RbException()<<"dnPhyloCTMC: can't specify 'pInv' if Q if a SiteMixtureModel";

	if ( site_matrices and site_matrices->getRevObject() != RevNullObject::getInstance() )
	    throw RbException()<<"dnPhyloCTMC: can't specify 'siteMatrices' if Q is a SiteMixtureModel";
    }
    else
    {
	RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rm = static_cast<const RateGenerator &>( q->getRevObject() ).getDagNode();
	dist->setRateMatrix( rm );
    }

    if ( site_ratesNode != NULL && site_ratesNode->getValue().size() > 0 )
    {
	dist->setSiteRates( site_ratesNode );
    }
        
    if ( site_rates_probsNode != NULL && site_rates_probsNode->getValue().size() > 0 )
    {
	dist->setSiteRatesProbs( site_rates_probsNode );
    }

    return dist;
}

RevBayesCore::TypedDistribution< RevBayesCore::AbstractHomologousDiscreteCharacterData >* Dist_phyloCTMC::createDistribution( void ) const
{
    const std::string& dt = static_cast<const RlString &>( type->getRevObject() ).getValue();

    // get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const Tree &>( tree->getRevObject() ).getDagNode();
    size_t n = size_t( static_cast<const Natural &>( nSites->getRevObject() ).getValue() );
    bool ambig = static_cast<const RlBoolean &>( treatAmbiguousAsGap->getRevObject() ).getValue();
    const std::string& code = static_cast<const RlString &>( coding->getRevObject() ).getValue();
    bool internal = static_cast<const RlBoolean &>( storeInternalNodes->getRevObject() ).getValue();
    bool gapmatch = static_cast<const RlBoolean &>( gapMatchClamped->getRevObject() ).getValue();
    bool invert = static_cast<const RlBoolean &>( invertProbability->getRevObject() ).getValue();

    if ( !(dt == "Binary" || dt == "Restriction" || dt == "Standard") && code != "all")
    {
        throw RbException( "Ascertainment bias correction only supported with Standard and Binary/Restriction datatypes" );
    }

    if ( dt == "DNA" )
    {
        RevBayesCore::PhyloCTMCSiteHomogeneousNucleotide<RevBayesCore::DnaState> *dist = new RevBayesCore::PhyloCTMCSiteHomogeneousNucleotide<RevBayesCore::DnaState>(tau, true, n, ambig, internal, gapmatch, invert);

	return setDistParameters(dist);
    }
    else if ( dt == "RNA" )
    {
        RevBayesCore::PhyloCTMCSiteHomogeneousNucleotide<RevBayesCore::RnaState> *dist = new RevBayesCore::PhyloCTMCSiteHomogeneousNucleotide<RevBayesCore::RnaState>(tau, true, n, ambig, internal, gapmatch, invert);

        return setDistParameters(dist);
    }
    else if ( dt == "AA" || dt == "Protein" )
    {
        RevBayesCore::PhyloCTMCSiteHomogeneous<RevBayesCore::AminoAcidState> *dist = new RevBayesCore::PhyloCTMCSiteHomogeneous<RevBayesCore::AminoAcidState>(tau, 20, true, n, ambig, internal, gapmatch, invert);

        return setDistParameters(dist);
    }
    else if ( dt == "Codon" )
    {
        RevBayesCore::PhyloCTMCSiteHomogeneous<RevBayesCore::CodonState> *dist = new RevBayesCore::PhyloCTMCSiteHomogeneous<RevBayesCore::CodonState>(tau, 61, true, n, ambig, internal, gapmatch, invert);
        
        return setDistParameters(dist);
    }
    else if ( dt == "Doublet" )
    {
        auto dist = new RevBayesCore::PhyloCTMCSiteHomogeneous<RevBayesCore::DoubletState>(tau, 16, true, n, ambig, internal, gapmatch, invert);

        return setDistParameters(dist);
    }
    else if ( dt == "PoMo" )
    {
        // we get the number of states from the rate matrix (we don't know, because PoMo is flexible about its rates)
	    int nChars = computeNumberOfStates();

        RevBayesCore::PhyloCTMCSiteHomogeneous<RevBayesCore::PoMoState> *dist = new RevBayesCore::PhyloCTMCSiteHomogeneous<RevBayesCore::PoMoState>(tau, nChars, !true, n, ambig, internal, gapmatch, invert);

        return setDistParameters(dist);
    }
    else if ( dt == "Standard" )
    {
        // we get the number of states from the rates matrix
	    int nChars = computeNumberOfStates();

        int cd = RevBayesCore::AscertainmentBias::ALL;
        // split the coding option on "|"
        if (code == "informative")
        {
            cd = RevBayesCore::AscertainmentBias::INFORMATIVE;
        }
        else if (code == "nstates")
        {
            cd = RevBayesCore::AscertainmentBias::NSTATES;
        }
        else if (code == "variable")
        {
            cd = RevBayesCore::AscertainmentBias::VARIABLE;
        }
        else if (code != "all")
        {
            std::stringstream ss;
            ss << "Invalid coding option \"" << code << "\"\n";
            ss << "\tAvailable Standard state codings: all, informative, nstates, variable\n";
            ss << "\tDefault: all.\n";
            throw RbException(ss.str());
        }

        RevBayesCore::PhyloCTMCSiteHomogeneous<RevBayesCore::StandardState> *dist;
        if (cd == RevBayesCore::AscertainmentBias::ALL)
        {
            dist = new RevBayesCore::PhyloCTMCSiteHomogeneous<RevBayesCore::StandardState>(tau, nChars, true, n, ambig, internal, gapmatch, invert);
        }
        else
        {
            dist = new RevBayesCore::PhyloCTMCSiteHomogeneousConditional<RevBayesCore::StandardState>(tau, nChars, true, n, ambig, RevBayesCore::AscertainmentBias::Coding(cd), internal, gapmatch, invert);
        }

        return setDistParameters(dist);
    }
    else if ( dt == "NaturalNumbers" )
    {
        // we get the number of states from the rates matrix
	    int nChars = computeNumberOfStates();

        RevBayesCore::PhyloCTMCSiteHomogeneous<RevBayesCore::NaturalNumbersState> *dist = new RevBayesCore::PhyloCTMCSiteHomogeneous<RevBayesCore::NaturalNumbersState>(tau, nChars, true, n, ambig, internal, gapmatch, invert);

	return setDistParameters(dist);
    }
    else if ( dt == "Binary" || dt == "Restriction" )
    {
        // we get the number of states from the rates matrix
	    const int nChars = computeNumberOfStates();

        // sanity check
        if ( nChars != 2 )
        {
            throw RbException( "Only binary characters allowed for type=Binary/Restriction" );
        }

        // split the coding option on "|"
        std::string s = code;
        std::vector<std::string> tokens;

        size_t pos = 0;
        while ((pos = s.find("|")) != std::string::npos)
        {
            tokens.push_back(s.substr(0, pos));
            s.erase(0, pos + 1);
        }
        tokens.push_back(s);

        // set the flags for each token
        int cd = RevBayesCore::AscertainmentBias::ALL;
        for (size_t i = 0; i < tokens.size(); i++)
        {
            if (tokens[i] == "noabsencesites")
            {
                cd |= RevBayesCore::BinaryAscertainmentBias::NOABSENCESITES;
            }
            else if (tokens[i] == "nopresencesites")
            {
                cd |= RevBayesCore::BinaryAscertainmentBias::NOPRESENCESITES;
            }
            else if (tokens[i] == "nstates")
            {
                cd |= RevBayesCore::AscertainmentBias::NSTATES;
            }
            else if (tokens[i] == "informative")
            {
                cd |= RevBayesCore::AscertainmentBias::INFORMATIVE;
            }
            else if (tokens[i] == "variable")
            {
                cd |= RevBayesCore::AscertainmentBias::VARIABLE;
            }
            else if (tokens[i] == "nosingletonpresence")
            {
                cd |= RevBayesCore::BinaryAscertainmentBias::NOSINGLETONPRESENCE;
            }
            else if (tokens[i] == "nosingletonabsence")
            {
                cd |= RevBayesCore::BinaryAscertainmentBias::NOSINGLETONABSENCE;
            }
            else if (tokens[i] == "nosingletons")
            {
                cd |= RevBayesCore::BinaryAscertainmentBias::NOSINGLETONS;
            }
            else if (tokens[i] != "all")
            {
                std::stringstream ss;
                ss << "Invalid coding option \"" << tokens[i] << "\"\n";
                ss << "\tAvailable Binary state codings: all, noabsencesites, nopresencesites, informative, variable, nosingletonpresence, nosingletonabsence, nosingletons\n";
                ss << "\tDefault: all. Codings are combined using the vertical bar \'|\'\n";
                throw RbException(ss.str());
            }
        }

        RevBayesCore::PhyloCTMCSiteHomogeneousBinary *dist = new RevBayesCore::PhyloCTMCSiteHomogeneousBinary(tau, true, (size_t)n, ambig, RevBayesCore::BinaryAscertainmentBias::Coding(cd), internal, gapmatch, invert);

        return setDistParameters(dist);
    }

    return nullptr;
}



/* Get Rev type of object */
const std::string& Dist_phyloCTMC::getClassType(void)
{

    static std::string rev_type = "Dist_phyloCTMC";

    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_phyloCTMC::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<AbstractHomologousDiscreteCharacterData>::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_phyloCTMC::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhyloCTMC";
    
    return d_name;
}


MethodTable Dist_phyloCTMC::getDistributionMethods( void ) const
{
    
    MethodTable methods = TypedDistribution<AbstractHomologousDiscreteCharacterData>::getDistributionMethods();
    
    // member functions
    ArgumentRules* siteLikelihoodsArgRules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_phyloCTMC, ModelVector<Real> >( "siteLikelihoods", variable, siteLikelihoodsArgRules, true ) );
    
    // f(Xh,rk|theta); returning a matrixreal where each row corresponds a site and each column corresponds an array (length = (num_site_rates + pInv != 0)) of rate categories [site_rate_index + pInv != 0] with the first element under inv if pInv > 0
    ArgumentRules* siteRateLikelihoodsArgRules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_phyloCTMC, MatrixReal >( "siteRateLikelihoods", variable, siteRateLikelihoodsArgRules, true ) );
    
    // f(Xh,mk|theta); returning a matrixreal where each row corresponds a site and each column corresponds an array (length = ((num_site_rates + pInv != 0) * num_site_matrices)) of mixture categories [(site_rate_index + pInv != 0) * num_site_matrices + site_matrix_index] with the first num_site_matrices elements under inv if pInv > 0
    ArgumentRules* siteMixtureLikelihoodsArgRules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_phyloCTMC, MatrixReal >( "siteMixtureLikelihoods", variable, siteMixtureLikelihoodsArgRules, true ) );
    
    ArgumentRules* siteRatesArgRules = new ArgumentRules();
    
    std::vector<std::string> optionsMethod;
    optionsMethod.push_back( "sampling" );
    optionsMethod.push_back( "weightedAverage" );
    siteRatesArgRules->push_back( new OptionRule( "estimateMethod", new RlString("sampling"), optionsMethod, "The method used to estimate the site specific rate." ) );
    
    methods.addFunction( new DistributionMemberFunction<Dist_phyloCTMC, ModelVector<RealPos> >( "siteRates", variable, siteRatesArgRules, true ) );
    
    return methods;
}


/** Return member rules (no members) */
const MemberRules& Dist_phyloCTMC::getParameterRules(void) const
{

    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        dist_member_rules.push_back( new ArgumentRule( "tree", Tree::getClassTypeSpec(), "The tree along which the process evolves.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        std::vector<TypeSpec> rateMatrixTypes;
        rateMatrixTypes.push_back( RateGenerator::getClassTypeSpec() );
        rateMatrixTypes.push_back( ModelVector<RateGenerator>::getClassTypeSpec() );
        rateMatrixTypes.push_back( SiteMixtureModel::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "Q", rateMatrixTypes, "The global, branch-specific or site-mixture rate matrices.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        // optional argument for the root frequencies
        dist_member_rules.push_back( new ArgumentRule( "rootFrequencies", Simplex::getClassTypeSpec(), "The root specific frequencies of the characters, if applicable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

        std::vector<TypeSpec> branchRateTypes;
        branchRateTypes.push_back( RealPos::getClassTypeSpec() );
        branchRateTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "branchRates", branchRateTypes, "The global or branch-specific rate multipliers.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );

        //dist_member_rules.push_back( new ArgumentRule( "siteMatrices", RlBoolean::getClassTypeSpec(), "Treat Q as vector of site mixture categories instead of branch-specific matrices?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );
        std::vector<TypeSpec> matrix_probs_types;
        matrix_probs_types.push_back(Simplex::getClassTypeSpec());
        matrix_probs_types.push_back(RlBoolean::getClassTypeSpec());

        dist_member_rules.push_back( new ArgumentRule( "siteMatrices", matrix_probs_types, "Simplex of site matrix mixture probabilities. Treats Q as vector of site mixture categories instead of branch-specific matrices.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "siteRates", ModelVector<RealPos>::getClassTypeSpec(), "The rate categories for the sites.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "siteRatesProbs", Simplex::getClassTypeSpec(), "The probability weights of rate categories for the sites.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "pInv", Probability::getClassTypeSpec(), "The probability of a site being invariant.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

        dist_member_rules.push_back( new ArgumentRule( "nSites", Natural::getClassTypeSpec(), "The number of sites, used for simulation.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural() ) );

        std::vector<std::string> options;
        options.push_back( "DNA" );
        options.push_back( "RNA" );
        options.push_back( "AA" );
        options.push_back( "Codon" );
        options.push_back( "Doublet" );
        options.push_back( "PoMo" );
        options.push_back( "Protein" );
        options.push_back( "Standard" );
        options.push_back( "NaturalNumbers" );
        options.push_back( "Binary" );
        options.push_back( "Restriction" );
        dist_member_rules.push_back( new OptionRule( "type", new RlString("DNA"), options, "The data type, used for simulation and initialization." ) );

        dist_member_rules.push_back( new ArgumentRule( "treatAmbiguousAsGap", RlBoolean::getClassTypeSpec(), "Should we treat ambiguous characters as gaps/missing?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );

        dist_member_rules.push_back( new ArgumentRule("coding", RlString::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("all") ) );

        dist_member_rules.push_back( new ArgumentRule( "storeInternalNodes", RlBoolean::getClassTypeSpec(), "Should we store internal node states in the character matrix?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );
        
        dist_member_rules.push_back( new ArgumentRule( "gapMatchClamped", RlBoolean::getClassTypeSpec(), "Should we set the simulated character to be gap or missing if the corresponding character in the clamped matrix is gap or missing?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( true ) ) );

        dist_member_rules.push_back( new ArgumentRule("invertProbability", RlBoolean::getClassTypeSpec(), "Should we return the reciprocal of the probability?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );

        rules_set = true;
    }

    return dist_member_rules;
}


const TypeSpec& Dist_phyloCTMC::getTypeSpec( void ) const
{

    static TypeSpec ts = getClassTypeSpec();

    return ts;
}


/** Print value for user */
void Dist_phyloCTMC::printValue(std::ostream& o) const
{

    o << "Character-State-Evolution-Along-Tree Process(tree=";
    if ( tree != NULL )
    {
        o << tree->getName();
    }
    else
    {
        o << "?";
    }
    o << ", Q=";
    if ( q != NULL )
    {
        o << q->getName();
    }
    else
    {
        o << "?";
    }
    o << ", branchRates=";
    if ( rate != NULL )
    {
        o << rate->getName();
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
    o << ", p_inv=";
    if ( p_inv != NULL )
    {
        o << p_inv->getName();
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
void Dist_phyloCTMC::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "Q" )
    {
        q = var;
    }
    else if ( name == "rootFrequencies" )
    {
        root_frequencies = var;
    }
    else if ( name == "branchRates" )
    {
        rate = var;
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
    else if ( name == "pInv" )
    {
        p_inv = var;
    }
    else if ( name == "nSites" )
    {
        nSites = var;
    }
    else if ( name == "type" )
    {
        type = var;
    }
    else if ( name == "treatAmbiguousAsGap" )
    {
        treatAmbiguousAsGap = var;
    }
    else if ( name == "storeInternalNodes" )
    {
        storeInternalNodes = var;
    }
    else if ( name == "gapMatchClamped" )
    {
        gapMatchClamped = var;
    }
    else if ( name == "coding" )
    {
        coding = var;
    }
    else if ( name == "invertProbability" )
    {
        invertProbability = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }

}


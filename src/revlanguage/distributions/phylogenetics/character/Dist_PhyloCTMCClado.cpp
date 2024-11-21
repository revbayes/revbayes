#include <math.h>
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "Dist_phyloCTMCClado.h"
#include "CladogeneticProbabilityMatrix.h"
#include "OptionRule.h"
#include "PhyloCTMCClado.h"
#include "Probability.h"
#include "RevNullObject.h"
#include "RlCladogeneticProbabilityMatrix.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "RlString.h"
#include "RlTree.h"
#include "StandardState.h"
#include "NaturalNumbersState.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "AbstractPhyloCTMCSiteHomogeneous.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagMemberFunction.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DiscreteTaxonData.h"
#include "DynamicNode.h"
#include "HomologousDiscreteCharacterData.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RateGenerator.h"
#include "RbBoolean.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlBoolean.h"
#include "RlConstantNode.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistribution.h"
#include "RlRateGenerator.h"
#include "RlSiteMixtureModel.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "Simplex.h"
#include "StochasticNode.h"
#include "StringUtilities.h"
#include "TaxaState.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

Dist_phyloCTMCClado::Dist_phyloCTMCClado() : TypedDistribution< AbstractHomologousDiscreteCharacterData >()
{
    
}

Dist_phyloCTMCClado::~Dist_phyloCTMCClado()
{
    
}


Dist_phyloCTMCClado* Dist_phyloCTMCClado::clone( void ) const
{
    
    return new Dist_phyloCTMCClado(*this);
}

template <typename State>
RevBayesCore::TypedDistribution< RevBayesCore::AbstractHomologousDiscreteCharacterData >* Dist_phyloCTMCClado::createDistributionForState( void ) const
{
    // get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const Tree &>( tree->getRevObject() ).getDagNode();
    size_t n = size_t( static_cast<const Natural &>( nSites->getRevObject() ).getValue() );
    bool ambig = static_cast<const RlBoolean &>( treatAmbiguousAsGap->getRevObject() ).getValue();

    bool internal = static_cast<const RlBoolean &>( storeInternalNodes->getRevObject() ).getValue();
    bool gapmatch = static_cast<const RlBoolean &>( gapMatchClamped->getRevObject() ).getValue();

    size_t nNodes = tau->getValue().getNumberOfNodes();

    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* site_ratesNode = NULL;
    if ( site_rates != NULL && site_rates->getRevObject() != RevNullObject::getInstance() )
    {
        site_ratesNode = static_cast<const ModelVector<RealPos> &>( site_rates->getRevObject() ).getDagNode();
    }
    RevBayesCore::TypedDagNode< double >* p_invNode = NULL;
    if ( p_inv != NULL && p_inv->getRevObject() != RevNullObject::getInstance() )
    {
        p_invNode = static_cast<const Probability &>( p_inv->getRevObject() ).getDagNode();
    }

    RevBayesCore::TypedDistribution< RevBayesCore::AbstractHomologousDiscreteCharacterData > *d = NULL;
    const RevBayesCore::TypedDagNode< RevBayesCore::Simplex > *rf = NULL;
    if ( root_frequencies != NULL && root_frequencies->getRevObject() != RevNullObject::getInstance() )
    {
        rf = static_cast<const Simplex &>( root_frequencies->getRevObject() ).getDagNode();
    }

    // get basic parameters that do not depend on data type
    // we get the number of states from the rates matrix
    // set the rate matrix
    size_t nChars = 1;
    size_t nCharsClado = 1;
    if ( q->getRevObject().isType( ModelVector<RateGenerator>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( q->getRevObject() ).getDagNode();
        nChars = rm->getValue()[0].getNumberOfStates();
    }
    else if ( q->getRevObject().isType( SiteMixtureModel::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::SiteMixtureModel >* mm = static_cast<const SiteMixtureModel &>( q->getRevObject() ).getDagNode();
        nChars = mm->getValue().getNumberOfStates();
    }
    else
    {
        RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rm = static_cast<const RateMatrix &>( q->getRevObject() ).getDagNode();
        nChars = rm->getValue().getNumberOfStates();
    }
    if ( cladoProbs->getRevObject().isType( ModelVector<CladogeneticProbabilityMatrix>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::CladogeneticProbabilityMatrix> >* cp = static_cast<const ModelVector<CladogeneticProbabilityMatrix> &>( q->getRevObject() ).getDagNode();
        nCharsClado = cp->getValue()[0].getNumberOfStates();
    }
    else
    {
        RevBayesCore::TypedDagNode<RevBayesCore::CladogeneticProbabilityMatrix>* cp = static_cast<const CladogeneticProbabilityMatrix &>( cladoProbs->getRevObject() ).getDagNode();
        nCharsClado = cp->getValue().getNumberOfStates();
    }
    // state space size checks
    if (rf != NULL) {
        size_t rf_size = rf->getValue().size();
        if (nChars != rf_size) {
            throw RbException("The root frequencies vector and rate matrix do not have the same number of states.\n");
        }
        if (nCharsClado != rf_size) {
            throw RbException("The root frequencies vector and cladogenetic probabilities do not have the same number of states.\n");
        }
        
    }

    if (nChars != nCharsClado) {
        throw RbException("The cladogenetic probabilities and rate matrix do not have the same number of states.\n");
    }

    RevBayesCore::PhyloCTMCClado<State> *dist = new RevBayesCore::PhyloCTMCClado<State>(tau, nChars, true, n, ambig, internal, gapmatch);

    // set the root frequencies (by default these are NULL so this is OK)
    dist->setRootFrequencies( rf );

    // set the probability for invariant site (by default this p_inv=0.0)
    dist->setPInv( p_invNode );

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

    // set the rate matrix
    if ( q->getRevObject().isType( ModelVector<RateGenerator>::getClassTypeSpec() ) )
    {
	RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( q->getRevObject() ).getDagNode();

	// sanity check
	if ( (nNodes-1) != rm->getValue().size() )
	{
	    throw RbException( "The number of substitution matrices does not match the number of branches" );
	}

	dist->setRateMatrix( rm );
    }
    else if ( q->getRevObject().isType( SiteMixtureModel::getClassTypeSpec() ) )
    {
	RevBayesCore::TypedDagNode< RevBayesCore::SiteMixtureModel >* mm = static_cast<const SiteMixtureModel &>( q->getRevObject() ).getDagNode();
	dist->setMixtureModel( mm );

	if ( root_frequencies and root_frequencies->getRevObject() != RevNullObject::getInstance() )
	    throw RbException()<<"dnPhyloCTMCClado: can't specify 'rootFrequencies' if Q if a SiteMixtureModel";

	if ( site_rates and site_rates->getRevObject() != RevNullObject::getInstance() )
	    throw RbException()<<"dnPhyloCTMCClado: can't specify 'siteRates' if Q if a SiteMixtureModel";

	if ( p_inv and p_inv->getRevObject() != RevNullObject::getInstance() )
	    throw RbException()<<"dnPhyloCTMCClado: can't specify 'pInv' if Q if a SiteMixtureModel";
    }
    else
    {
	RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rm = static_cast<const RateGenerator &>( q->getRevObject() ).getDagNode();
	dist->setRateMatrix( rm );
    }

    // set the clado probs
    if ( cladoProbs->getRevObject().isType( ModelVector<CladogeneticProbabilityMatrix>::getClassTypeSpec() ) )
    {
	RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::CladogeneticProbabilityMatrix> >* cp = static_cast<const ModelVector<CladogeneticProbabilityMatrix> &>( cladoProbs->getRevObject() ).getDagNode();

	// sanity check
	if ( (nNodes-1) != cp->getValue().size() )
	{
	    throw RbException( "The number of cladogenesis probability matrices does not match the number of branches" );
	}
	dist->setCladogenesisMatrix( cp );
    }
    else
    {
	RevBayesCore::TypedDagNode<RevBayesCore::CladogeneticProbabilityMatrix>* cp = static_cast<const CladogeneticProbabilityMatrix &>( cladoProbs->getRevObject() ).getDagNode();
	dist->setCladogenesisMatrix( cp );
    }

    if ( site_ratesNode != NULL && site_ratesNode->getValue().size() > 0 )
    {
	dist->setSiteRates( site_ratesNode );
    }

    d = dist;

    return d;
}

RevBayesCore::TypedDistribution< RevBayesCore::AbstractHomologousDiscreteCharacterData >* Dist_phyloCTMCClado::createDistribution( void ) const
{
    const std::string& dt = static_cast<const RlString &>( type->getRevObject() ).getValue();
    
    if ( dt == "NaturalNumbers" )
    {
	return createDistributionForState<RevBayesCore::NaturalNumbersState>();
    }
    else if ( dt == "Taxa" )
    {
	return createDistributionForState<RevBayesCore::TaxaState>();
    }
    else if (dt == "Standard")
    {
	return createDistributionForState<RevBayesCore::StandardState>();
    }
    else
    {
        throw RbException() << "Datatype \"" << dt << "\" not currently supported!";
    }
}



/* Get Rev type of object */
const std::string& Dist_phyloCTMCClado::getClassType(void) {
    
    static std::string rev_type = "Dist_phyloCTMCClado";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_phyloCTMCClado::getClassTypeSpec(void) {
    
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
std::string Dist_phyloCTMCClado::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhyloCTMCClado";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_phyloCTMCClado::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        // epoch model requires time tree
        dist_member_rules.push_back( new ArgumentRule( "tree"           , Tree::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        // epoch model requires vector of Q
        std::vector<TypeSpec> rateMatrixTypes;
        rateMatrixTypes.push_back( RateGenerator::getClassTypeSpec() );
        rateMatrixTypes.push_back( ModelVector<RateGenerator>::getClassTypeSpec() );
        rateMatrixTypes.push_back( SiteMixtureModel::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "Q"              , rateMatrixTypes, "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        // clado model accepts a single or vector of cladogenesis probs
        std::vector<TypeSpec> cladoMtxTypes;
        cladoMtxTypes.push_back( CladogeneticProbabilityMatrix::getClassTypeSpec() );
        cladoMtxTypes.push_back( ModelVector<CladogeneticProbabilityMatrix>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "cladoProbs"              , cladoMtxTypes, "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::DETERMINISTIC ) );
        
        // optional argument for the root frequencies
        dist_member_rules.push_back( new ArgumentRule( "rootFrequencies", Simplex::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        
        // branch rates
        std::vector<TypeSpec> branchRateTypes;
        branchRateTypes.push_back( RealPos::getClassTypeSpec() );
        branchRateTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "branchRates"    , branchRateTypes, "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        
        
        ModelVector<RealPos> *defaultSiteRates = new ModelVector<RealPos>();
        dist_member_rules.push_back( new ArgumentRule( "siteRates"      , ModelVector<RealPos>::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, defaultSiteRates ) );
        dist_member_rules.push_back( new ArgumentRule( "pInv"           , Probability::getClassTypeSpec()         , "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability(0.0) ) );
        
        dist_member_rules.push_back( new ArgumentRule( "nSites"         , Natural::getClassTypeSpec()             , "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(10) ) );
        
        
        std::vector<std::string> options;
        options.push_back( "DNA" );
        options.push_back( "RNA" );
        options.push_back( "AA" );
        options.push_back( "PoMo" );
        options.push_back( "Protein" );
        options.push_back( "Standard" );
        options.push_back( "NaturalNumbers" );
        options.push_back( "Taxa" );
        dist_member_rules.push_back( new OptionRule( "type", new RlString("NaturalNumbers"), options, "" ) );
        
        dist_member_rules.push_back( new ArgumentRule( "treatAmbiguousAsGap", RlBoolean::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );
        
        dist_member_rules.push_back( new ArgumentRule( "storeInternalNodes", RlBoolean::getClassTypeSpec(), "Should we store internal node states in the character matrix?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );
        
        dist_member_rules.push_back( new ArgumentRule( "gapMatchClamped", RlBoolean::getClassTypeSpec(), "Should we set the simulated character to gap or missing if the corresponding character in the clamped matrix is gap or missing?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( true ) ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}

const TypeSpec& Dist_phyloCTMCClado::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_phyloCTMCClado::printValue(std::ostream& o) const
{
    
    o << "Character-State-Evolution-Along-Tree-With-Cladogenesis Process(tree=";
    if ( tree != NULL ) {
        o << tree->getName();
    } else {
        o << "?";
    }
    o << ", Q=";
    o << ", cladoProbs=";
    if ( q != NULL ) {
        o << q->getName();
    } else {
        o << "?";
    }
    o << ", branchRates=";
    if ( rate != NULL ) {
        o << rate->getName();
    } else {
        o << "?";
    }
    o << ", site_rates=";
    if ( site_rates != NULL ) {
        o << site_rates->getName();
    } else {
        o << "?";
    }
    o << ", p_inv=";
    if ( p_inv != NULL ) {
        o << p_inv->getName();
    } else {
        o << "?";
    }
    o << ", nSites=";
    if ( nSites != NULL ) {
        o << nSites->getName();
    } else {
        o << "?";
    }
    o << ")";
    
}


/** Set a member variable */
void Dist_phyloCTMCClado::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "Q" )
    {
        q = var;
    }
    else if ( name == "cladoProbs" )
    {
        cladoProbs = var;
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
    else if ( name == "pInv" )
    {
        p_inv = var;
    }
    else if ( name == "nSites" )
    {
        nSites = var;
    }
    else if ( name == "nStates" )
    {
        nStates = var;
    }
    else if ( name == "storeInternalNodes" )
    {
        storeInternalNodes = var;
    }
    else if ( name == "gapMatchClamped" )
    {
        gapMatchClamped = var;
    }
    else if ( name == "type" )
    {
        type = var;
    }
    else if ( name == "treatAmbiguousAsGap" )
    {
        treatAmbiguousAsGap = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
    
}



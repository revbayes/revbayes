#include "Dist_phyloCTMCDASequence.h"

#include <cstddef>
#include <ostream>

#include "GeneralTreeHistoryCtmc.h"
#include "OptionRule.h"
#include "Probability.h"
#include "RevNullObject.h"
#include "RlBoolean.h"
#include "RlCladogeneticProbabilityMatrix.h"
#include "RlDistributionMemberFunction.h"
#include "RlRateGeneratorSequence.h"
#include "RlString.h"
#include "RlTree.h"
#include "StandardState.h"
#include "RlSimplex.h"
#include "PoMoState.h"
#include "NaturalNumbersState.h"
#include "RateGeneratorSequence.h"
#include "AminoAcidState.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DiscreteTaxonData.h"
#include "DistributionMemberFunction.h"
#include "DnaState.h"
#include "HomologousDiscreteCharacterData.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "Natural.h"
#include "PathRejectionSampleProposal.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RealPos.h"
#include "RlConstantNode.h"
#include "RlDistribution.h"
#include "RnaState.h"
#include "Simplex.h"
#include "StringUtilities.h"
#include "Tree.h"
#include "TreeHistoryCtmc.h"
#include "TypeSpec.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class CladogeneticProbabilityMatrix; }

using namespace RevLanguage;

Dist_phyloCTMCDASequence::Dist_phyloCTMCDASequence() : TypedDistribution< AbstractHomologousDiscreteCharacterData >()
{
    
}


Dist_phyloCTMCDASequence::~Dist_phyloCTMCDASequence()
{
    
}



Dist_phyloCTMCDASequence* Dist_phyloCTMCDASequence::clone( void ) const
{
    
    return new Dist_phyloCTMCDASequence(*this);
}


RevBayesCore::TypedDistribution< RevBayesCore::AbstractHomologousDiscreteCharacterData >* Dist_phyloCTMCDASequence::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const Tree &>( tree->getRevObject() ).getDagNode();
    size_t n = size_t( static_cast<const Natural &>( nSites->getRevObject() ).getValue() );
    const std::string& dt = static_cast<const RlString &>( type->getRevObject() ).getValue();
    bool ambig = static_cast<const RlBoolean &>( treatAmbiguousAsGap->getRevObject() ).getDagNode();
    size_t nNodes = tau->getValue().getNumberOfNodes();
    
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* siteRatesNode = NULL;
    if ( siteRates != NULL && siteRates->getRevObject() != RevNullObject::getInstance() )
    {
        siteRatesNode = static_cast<const ModelVector<RealPos> &>( siteRates->getRevObject() ).getDagNode();
    }
    RevBayesCore::TypedDagNode< double >* pInvNode = NULL;
    if ( pInv != NULL && pInv->getRevObject() != RevNullObject::getInstance() )
    {
        pInvNode = static_cast<const Probability &>( pInv->getRevObject() ).getDagNode();
    }
    
    RevBayesCore::TypedDistribution< RevBayesCore::AbstractHomologousDiscreteCharacterData > *d = NULL;
    const RevBayesCore::TypedDagNode< RevBayesCore::Simplex > *rf = NULL;
    if ( rootFrequencies->getRevObject() != RevNullObject::getInstance() )
    {
        rf = static_cast<const Simplex &>( rootFrequencies->getRevObject() ).getDagNode();
    }
    
    bool internal = static_cast<const RlBoolean &>( storeInternalNodes->getRevObject() ).getDagNode();
    
    if ( dt == "DNA" )
    {
        RevBayesCore::GeneralTreeHistoryCtmc<RevBayesCore::DnaState> *dist =
        new RevBayesCore::GeneralTreeHistoryCtmc<RevBayesCore::DnaState>(tau, 4, n, ambig, internal);
        
        // set the root frequencies (by default these are NULL so this is OK)
        dist->setRootFrequencies( rf );
        
        // set the clock rates
        if ( rate->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* clockRates = static_cast<const ModelVector<RealPos> &>( rate->getRevObject() ).getDagNode();
            
            // sanity check
            size_t nRates = clockRates->getValue().size();
            if ( (nNodes-1) != nRates )
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
        if ( q->getRevObject().isType( ModelVector<RateGeneratorSequence>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGeneratorSequence> >* rm = static_cast<const ModelVector<RateGeneratorSequence> &>( q->getRevObject() ).getDagNode();
            
            // sanity check
            if ( (nNodes-1) != rm->getValue().size() )
            {
                throw RbException( "The number of substitution matrices does not match the number of branches" );
            }
            
            // sanity check
            if ( rootFrequencies == NULL || rootFrequencies->getRevObject() == RevNullObject::getInstance() )
            {
                throw RbException( "If you provide branch-heterogeneous substitution matrices, then you also need to provide root frequencies." );
            }
            
            dist->setRateGenerator( rm );
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGeneratorSequence>* rm = static_cast<const RateGeneratorSequence &>( q->getRevObject() ).getDagNode();
            dist->setRateGenerator( rm );
        }
        
        if ( siteRatesNode != NULL && siteRatesNode->getValue().size() > 0 )
        {
            dist->setSiteRates( siteRatesNode );
        }
        
        d = dist;
    }
    else if ( dt == "RNA" )
    {
        RevBayesCore::GeneralTreeHistoryCtmc<RevBayesCore::RnaState> *dist = new RevBayesCore::GeneralTreeHistoryCtmc<RevBayesCore::RnaState>(tau, 4, n, ambig, internal);
        
        // set the root frequencies (by default these are NULL so this is OK)
        dist->setRootFrequencies( rf );
        
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
        if ( q->getRevObject().isType( ModelVector<RateGeneratorSequence>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGeneratorSequence> >* rm = static_cast<const ModelVector<RateGeneratorSequence> &>( q->getRevObject() ).getDagNode();
            
            // sanity check
            if ( (nNodes-1) != rm->getValue().size() )
            {
                throw RbException( "The number of substitution matrices does not match the number of branches" );
            }
            // sanity check
            if ( rootFrequencies == NULL || rootFrequencies->getRevObject() == RevNullObject::getInstance() )
            {
                throw RbException( "If you provide branch-heterogeneous substitution matrices, then you also need to provide root frequencies." );
            }
            
            
            dist->setRateGenerator( rm );
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGeneratorSequence>* rm = static_cast<const RateGeneratorSequence &>( q->getRevObject() ).getDagNode();
            dist->setRateGenerator( rm );
        }
        
        if ( siteRatesNode != NULL && siteRatesNode->getValue().size() > 0 )
        {
            dist->setSiteRates( siteRatesNode );
        }
        
        d = dist;
    }
    else if ( dt == "AA" || dt == "Protein" )
    {
        RevBayesCore::GeneralTreeHistoryCtmc<RevBayesCore::AminoAcidState> *dist = new RevBayesCore::GeneralTreeHistoryCtmc<RevBayesCore::AminoAcidState>(tau, 20, n, ambig, internal);
        // set the root frequencies (by default these are NULL so this is OK)
        dist->setRootFrequencies( rf );
        
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
        if ( q->getRevObject().isType( ModelVector<RateGeneratorSequence>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGeneratorSequence> >* rm = static_cast<const ModelVector<RateGeneratorSequence> &>( q->getRevObject() ).getDagNode();
            
            // sanity check
            if ( (nNodes-1) != rm->getValue().size() )
            {
                throw RbException( "The number of substitution matrices does not match the number of branches" );
            }
            // sanity check
            if ( rootFrequencies == NULL || rootFrequencies->getRevObject() == RevNullObject::getInstance() )
            {
                throw RbException( "If you provide branch-heterogeneous substitution matrices, then you also need to provide root frequencies." );
            }
            
            
            dist->setRateGenerator( rm );
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGeneratorSequence>* rm = static_cast<const RateGeneratorSequence &>( q->getRevObject() ).getDagNode();
            dist->setRateGenerator( rm );
        }
        
        if ( siteRatesNode != NULL && siteRatesNode->getValue().size() > 0 )
        {
            dist->setSiteRates( siteRatesNode );
        }
        
        d = dist;
    }
    else if ( dt == "PoMo" )
    {
        
        // we get the number of states from the rate matrix (we don't know, because PoMo is flexible about its rates)
        // set the rate matrix
        size_t nChars = 1;
        if ( q->getRevObject().isType( ModelVector<RateGeneratorSequence>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGeneratorSequence> >* rm = static_cast<const ModelVector<RateGeneratorSequence> &>( q->getRevObject() ).getDagNode();
            nChars = rm->getValue()[0].getNumberOfStates();
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGeneratorSequence>* rm = static_cast<const RateGeneratorSequence &>( q->getRevObject() ).getDagNode();
            nChars = rm->getValue().getNumberOfStates();
        }
        
        RevBayesCore::GeneralTreeHistoryCtmc<RevBayesCore::PoMoState> *dist = new RevBayesCore::GeneralTreeHistoryCtmc<RevBayesCore::PoMoState>(tau, nChars, n, ambig, internal);
        
        // set the root frequencies (by default these are NULL so this is OK)
        dist->setRootFrequencies( rf );
        
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
        if ( q->getRevObject().isType( ModelVector<RateGeneratorSequence>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGeneratorSequence> >* rm = static_cast<const ModelVector<RateGeneratorSequence> &>( q->getRevObject() ).getDagNode();
            
            // sanity check
            if ( (nNodes-1) != rm->getValue().size() )
            {
                throw RbException( "The number of substitution matrices does not match the number of branches" );
            }
            // sanity check
            if ( rootFrequencies == NULL || rootFrequencies->getRevObject() == RevNullObject::getInstance() )
            {
                throw RbException( "If you provide branch-heterogeneous substitution matrices, then you also need to provide root frequencies." );
            }
            
            
            dist->setRateGenerator( rm );
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGeneratorSequence>* rm = static_cast<const RateGeneratorSequence &>( q->getRevObject() ).getDagNode();
            dist->setRateGenerator( rm );
        }
        
        if ( siteRatesNode != NULL && siteRatesNode->getValue().size() > 0 )
        {
            dist->setSiteRates( siteRatesNode );
        }
        
        d = dist;
    }
    else if ( dt == "Standard" )
    {
        // we get the number of states from the rates matrix
        // set the rate matrix
        size_t nChars = 1;
        if ( q->getRevObject().isType( ModelVector<RateGeneratorSequence>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGeneratorSequence> >* rm = static_cast<const ModelVector<RateGeneratorSequence> &>( q->getRevObject() ).getDagNode();
            nChars = rm->getValue()[0].getNumberOfStates();
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGeneratorSequence>* rm = static_cast<const RateGeneratorSequence &>( q->getRevObject() ).getDagNode();
            nChars = rm->getValue().getNumberOfStates();
        }
        
        RevBayesCore::GeneralTreeHistoryCtmc<RevBayesCore::StandardState> *dist = new RevBayesCore::GeneralTreeHistoryCtmc<RevBayesCore::StandardState>(tau, nChars, n, ambig, internal);
        
        // set the root frequencies (by default these are NULL so this is OK)
        dist->setRootFrequencies( rf );
        
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
        if ( q->getRevObject().isType( ModelVector<RateGeneratorSequence>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGeneratorSequence> >* rm = static_cast<const ModelVector<RateGeneratorSequence> &>( q->getRevObject() ).getDagNode();
            
            // sanity check
            if ( (nNodes-1) != rm->getValue().size() )
            {
                throw RbException( "The number of substitution matrices does not match the number of branches" );
            }
            // sanity check
            if ( rootFrequencies == NULL || rootFrequencies->getRevObject() == RevNullObject::getInstance() )
            {
                throw RbException( "If you provide branch-heterogeneous substitution matrices, then you also need to provide root frequencies." );
            }
            
            
            dist->setRateGenerator( rm );
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGeneratorSequence>* rm = static_cast<const RateGeneratorSequence &>( q->getRevObject() ).getDagNode();
            dist->setRateGenerator( rm );
        }
        
        // set cladogenetic transiiton probs
        RevBayesCore::TypedDagNode<RevBayesCore::CladogeneticProbabilityMatrix>* cp_tdn = NULL;
        RevBayesCore::DeterministicNode<RevBayesCore::CladogeneticProbabilityMatrix>* cp_dn = NULL;
        if ( clado_probs != NULL && clado_probs->getRevObject() != RevNullObject::getInstance() )
        {

            cp_tdn = static_cast<const CladogeneticProbabilityMatrix &>( clado_probs->getRevObject() ).getDagNode();
            cp_dn  = static_cast<RevBayesCore::DeterministicNode<RevBayesCore::CladogeneticProbabilityMatrix>*>( cp_tdn );
            dist->setCladogeneticProbabilityMatrix( cp_dn );
        }
        
        if ( siteRatesNode != NULL && siteRatesNode->getValue().size() > 0 )
        {
            dist->setSiteRates( siteRatesNode );
        }
        
        d = dist;
    }
    else if ( dt == "NaturalNumbers" )
    {
        // we get the number of states from the rates matrix
        size_t nChars = 1;
        if ( q->getRevObject().isType( ModelVector<RateGeneratorSequence>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGeneratorSequence> >* rm = static_cast<const ModelVector<RateGeneratorSequence> &>( q->getRevObject() ).getDagNode();
            nChars = rm->getValue()[0].getNumberOfStates();
            // RevBayesCore::g_MAX_NAT_NUM_STATES = nChars;
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGeneratorSequence>* rm = static_cast<const RateGeneratorSequence &>( q->getRevObject() ).getDagNode();
            nChars = rm->getValue().getNumberOfStates();
            // RevBayesCore::g_MAX_NAT_NUM_STATES = nChars;
        }
        
        RevBayesCore::GeneralTreeHistoryCtmc<RevBayesCore::NaturalNumbersState> *dist = new RevBayesCore::GeneralTreeHistoryCtmc<RevBayesCore::NaturalNumbersState>(tau, nChars, n, ambig, internal);
        
        // set the root frequencies (by default these are NULL so this is OK)
        dist->setRootFrequencies( rf );
        
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
        if ( q->getRevObject().isType( ModelVector<RateGeneratorSequence>::getClassTypeSpec() ) )
        {
            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGeneratorSequence> >* rm = static_cast<const ModelVector<RateGeneratorSequence> &>( q->getRevObject() ).getDagNode();
            
            // sanity check
            if ( (nNodes-1) != rm->getValue().size() )
            {
                throw RbException( "The number of substitution matrices does not match the number of branches" );
            }
            // sanity check
            if ( rootFrequencies == NULL || rootFrequencies->getRevObject() == RevNullObject::getInstance() )
            {
                throw RbException( "If you provide branch-heterogeneous substitution matrices, then you also need to provide root frequencies." );
            }
            
            
            dist->setRateGenerator( rm );
        }
        else
        {
            RevBayesCore::TypedDagNode<RevBayesCore::RateGeneratorSequence>* rm = static_cast<const RateGeneratorSequence &>( q->getRevObject() ).getDagNode();
            dist->setRateGenerator( rm );
        }
        
        if ( siteRatesNode != NULL && siteRatesNode->getValue().size() > 0 )
        {
            dist->setSiteRates( siteRatesNode );
        }
        
        d = dist;
    }
    
    
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_phyloCTMCDASequence::getClassType(void)
{
    
    static std::string revType = "Dist_phyloCTMCDASequence";
    
    return revType;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_phyloCTMCDASequence::getClassTypeSpec(void)
{
    
    static TypeSpec revTypeSpec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
    
    return revTypeSpec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_phyloCTMCDASequence::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhyloCTMCDASequence";
    
    return d_name;
}


MethodTable Dist_phyloCTMCDASequence::getDistributionMethods( void ) const
{
    
    MethodTable methods = TypedDistribution<AbstractHomologousDiscreteCharacterData>::getDistributionMethods();
    
    // member functions
    ArgumentRules* node_state_frequencies_arg_rules = new ArgumentRules();
    node_state_frequencies_arg_rules->push_back( new ArgumentRule( "node", Natural::getClassTypeSpec(), "The index of the node.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new DistributionMemberFunction<Dist_phyloCTMCDASequence, Simplex >( "stateFrequencies", variable, node_state_frequencies_arg_rules   ) );
    
    ArgumentRules* node_character_changes_arg_rules = new ArgumentRules();
    node_character_changes_arg_rules->push_back( new ArgumentRule( "node", Natural::getClassTypeSpec(), "The index of the node.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new DistributionMemberFunction<Dist_phyloCTMCDASequence, ModelVector<Natural> >( "numCharacterChanges", variable, node_character_changes_arg_rules   ) );
    
    ArgumentRules* relative_time_arg_rules = new ArgumentRules();
    relative_time_arg_rules->push_back( new ArgumentRule( "node", Natural::getClassTypeSpec(), "The index of the node.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    relative_time_arg_rules->push_back( new ArgumentRule( "site", Natural::getClassTypeSpec(), "The index of the site.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    methods.addFunction( new DistributionMemberFunction<Dist_phyloCTMCDASequence, Simplex >( "relativeTimeInStates", variable, relative_time_arg_rules   ) );
    
    //    ArgumentRules* avgExtinctionArgRules = new ArgumentRules();
    //    methods.addFunction( new DistributionMemberFunction<Dist_heterogeneousRateBirthDeath, ModelVector<RealPos> >( "averageExtinctionRate", variable, avgExtinctionArgRules   ) );
    
    return methods;
}


/** Return member rules (no members) */
const MemberRules& Dist_phyloCTMCDASequence::getParameterRules(void) const
{
    
    static MemberRules distMemberRules;
    static bool rulesSet = false;
    
    if ( !rulesSet )
    {
        distMemberRules.push_back( new ArgumentRule( "tree", Tree::getClassTypeSpec(), "The tree along which the process evolves.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        //
        std::vector<TypeSpec> rateMatrixTypes;
        rateMatrixTypes.push_back( RateGeneratorSequence::getClassTypeSpec() );
        rateMatrixTypes.push_back( ModelVector<RateGeneratorSequence>::getClassTypeSpec() );
        distMemberRules.push_back( new ArgumentRule( "Q", rateMatrixTypes, "The global or branch-specific rate matrices.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        // cladogenetic transition probabilities
        std::vector<TypeSpec> cladoMtxTypes;
        cladoMtxTypes.push_back( CladogeneticProbabilityMatrix::getClassTypeSpec() );
//        cladoMtxTypes.push_back( ModelVector<CladogeneticProbabilityMatrix>::getClassTypeSpec() );
        distMemberRules.push_back( new ArgumentRule( "cladoProbs", cladoMtxTypes, "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL) );

        // optional argument for the root frequencies
        distMemberRules.push_back( new ArgumentRule( "rootFrequencies", Simplex::getClassTypeSpec(), "The root specific frequencies of the characters, if applicable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        
        std::vector<TypeSpec> branchRateTypes;
        branchRateTypes.push_back( RealPos::getClassTypeSpec() );
        branchRateTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        distMemberRules.push_back( new ArgumentRule( "branchRates", branchRateTypes, "The global or branch-specific rate multipliers.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        
        ModelVector<RealPos> *defaultSiteRates = new ModelVector<RealPos>();
        distMemberRules.push_back( new ArgumentRule( "siteRates", ModelVector<RealPos>::getClassTypeSpec(), "The rate categories for the sites.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, defaultSiteRates ) );
        distMemberRules.push_back( new ArgumentRule( "pInv", Probability::getClassTypeSpec(), "The probability of a site being invariant.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Probability(0.0) ) );
        
        distMemberRules.push_back( new ArgumentRule( "nSites", Natural::getClassTypeSpec(), "The number of sites, used for simulation.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(10) ) );
        
        std::vector<std::string> options;
        options.push_back( "DNA" );
        options.push_back( "RNA" );
        options.push_back( "AA" );
        options.push_back( "PoMo" );
        options.push_back( "Protein" );
        options.push_back( "Standard" );
        options.push_back( "NaturalNumbers" );
        options.push_back( "Restriction" );
        distMemberRules.push_back( new OptionRule( "type", new RlString("DNA"), options, "The data type, used for simulation and initialization." ) );
    
        distMemberRules.push_back( new ArgumentRule( "storeInternalNodes", RlBoolean::getClassTypeSpec(), "Should we store internal node states in the character matrix?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );
        
        distMemberRules.push_back( new ArgumentRule( "treatAmbiguousAsGap", RlBoolean::getClassTypeSpec(), "Should we treat ambiguous characters as gaps/missing?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );
        
        distMemberRules.push_back( new ArgumentRule("coding", RlString::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("all") ) );
        
        rulesSet = true;
    }
    
    return distMemberRules;
}


const TypeSpec& Dist_phyloCTMCDASequence::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_phyloCTMCDASequence::printValue(std::ostream& o) const
{
    
    o << "Character-State-Evolution-Along-Tree Process(tree=";
    if ( tree != NULL ) {
        o << tree->getName();
    } else {
        o << "?";
    }
    o << ", Q=";
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
    o << ", siteRates=";
    if ( siteRates != NULL ) {
        o << siteRates->getName();
    } else {
        o << "?";
    }
    o << ", pInv=";
    if ( pInv != NULL ) {
        o << pInv->getName();
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
void Dist_phyloCTMCDASequence::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
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
        clado_probs = var;
    }
    else if ( name == "rootFrequencies" )
    {
        rootFrequencies = var;
    }
    else if ( name == "branchRates" )
    {
        rate = var;
    }
    else if ( name == "siteRates" )
    {
        siteRates = var;
    }
    else if ( name == "pInv" )
    {
        pInv = var;
    }
    else if ( name == "nSites" )
    {
        nSites = var;
    }
    else if ( name == "type" )
    {
        type = var;
    }
    else if ( name == "storeInternalNodes" )
    {
        storeInternalNodes = var;
    }
    else if ( name == "treatAmbiguousAsGap" )
    {
        treatAmbiguousAsGap = var;
    }
    else if ( name == "coding" )
    {
        coding = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
    
}

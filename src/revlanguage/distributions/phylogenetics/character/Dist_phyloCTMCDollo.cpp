#include "Dist_phyloCTMCDollo.h"

#include <cstddef>
#include <ostream>

#include "PhyloCTMCSiteHomogeneousConditional.h"
#include "PhyloCTMCSiteHomogeneousDollo.h"
#include "PhyloCTMCSiteHomogeneousDolloBinary.h"
#include "OptionRule.h"
#include "RevNullObject.h"
#include "RlBoolean.h"
#include "RlRateGenerator.h"
#include "RlString.h"
#include "RlTree.h"
#include "RlSimplex.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "Natural.h"
#include "PhyloCTMCSiteHomogeneousBinary.h"
#include "RateGenerator.h"
#include "RbBoolean.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RealPos.h"
#include "RlConstantNode.h"
#include "RlDistribution.h"
#include "StringUtilities.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class Simplex; }

using namespace RevLanguage;

Dist_phyloCTMCDollo::Dist_phyloCTMCDollo() : TypedDistribution< AbstractHomologousDiscreteCharacterData >()
{

}


Dist_phyloCTMCDollo::~Dist_phyloCTMCDollo()
{

}



Dist_phyloCTMCDollo* Dist_phyloCTMCDollo::clone( void ) const
{

    return new Dist_phyloCTMCDollo(*this);
}


RevBayesCore::TypedDistribution< RevBayesCore::AbstractHomologousDiscreteCharacterData >* Dist_phyloCTMCDollo::createDistribution( void ) const
{

    // get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const Tree &>( tree->getRevObject() ).getDagNode();
    size_t n = size_t( static_cast<const Natural &>( nSites->getRevObject() ).getValue() );
    const std::string& dt = static_cast<const RlString &>( type->getRevObject() ).getValue();
    bool norm = static_cast<const RlBoolean &>( normalize->getRevObject() ).getValue();
    size_t nNodes = tau->getValue().getNumberOfNodes();
    const std::string& code = static_cast<const RlString &>( coding->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<double>* deathRateNode = static_cast<const RealPos &>( deathRate->getRevObject() ).getDagNode();

    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* siteRatesNode = NULL;
    if ( siteRates != NULL && siteRates->getRevObject() != RevNullObject::getInstance() )
    {
        siteRatesNode = static_cast<const ModelVector<RealPos> &>( siteRates->getRevObject() ).getDagNode();
    }

    RevBayesCore::TypedDistribution< RevBayesCore::AbstractHomologousDiscreteCharacterData > *d = NULL;
    const RevBayesCore::TypedDagNode< RevBayesCore::Simplex > *rf = NULL;
    if ( rootFrequencies != NULL && rootFrequencies->getRevObject() != RevNullObject::getInstance() )
    {
        rf = static_cast<const Simplex &>( rootFrequencies->getRevObject() ).getDagNode();
    }

    // we get the number of states from the rates matrix
    // set the rate matrix
    size_t nChars = 1;
    if ( q->getRevObject().isType( ModelVector<RateGenerator>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( q->getRevObject() ).getDagNode();
        nChars = rm->getValue()[0].getNumberOfStates();
    }
    else if ( q->getRevObject() != RevNullObject::getInstance() )
    {
        RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rm = static_cast<const RateGenerator &>( q->getRevObject() ).getDagNode();
        nChars = rm->getValue().getNumberOfStates();
    }

    // split the coding option on "|"
    std::string s = code;
    std::vector<std::string> tokens;

    size_t pos = 0;
    while ((pos = s.find("|")) != std::string::npos) {
        tokens.push_back(s.substr(0, pos));
        s.erase(0, pos + 1);
    }
    tokens.push_back(s);

    RevBayesCore::PhyloCTMCSiteHomogeneousDollo *dist;

    if ( dt == "Standard" )
    {
        // set the flags for each token
        int cd = RevBayesCore::DolloAscertainmentBias::NOABSENCESITES;
        for (size_t i = 0; i < tokens.size(); i++)
        {
            if (tokens[i] == "informative")
            {
                cd = RevBayesCore::BinaryAscertainmentBias::INFORMATIVE;
            }
            else if (tokens[i] == "variable")
            {
                cd = RevBayesCore::BinaryAscertainmentBias::VARIABLE;
            }
            else if (tokens[i] != "noabsencesites")
            {
                std::stringstream ss;
                ss << "Unrecognized coding option \"" << tokens[i] << "\" for datatype=Standard\n";
                ss << "\tAvailable codings: noabsencesites, informative, variable\n";
                throw RbException(ss.str());
            }
        }

        dist = new RevBayesCore::PhyloCTMCSiteHomogeneousDollo(tau, nChars, true, n, true, RevBayesCore::DolloAscertainmentBias::Coding(cd), norm);

    }
    else
    {
        // set the flags for each token
        int cd = RevBayesCore::DolloBinaryAscertainmentBias::NOABSENCESITES;
        for (size_t i = 0; i < tokens.size(); i++)
        {
            if (tokens[i] == "nopresencesites")
            {
                cd |= RevBayesCore::BinaryAscertainmentBias::NOPRESENCESITES;
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
            else if (tokens[i] != "noabsencesites")
            {
                std::stringstream ss;
                ss << "Invalid coding option \"" << tokens[i] << "\" datatype=Binary/Restriction\n";
                ss << "\tAvailable Binary state codings: noabsencesites, nopresencesites, informative, variable, nosingletonpresence, nosingletonabsence, nosingletons\n";
                ss << "\tDefault: noabsencesites. Codings are combined using the vertical bar \'|\'\n";
                throw RbException(ss.str());
            }
        }

        dist = new RevBayesCore::PhyloCTMCSiteHomogeneousDolloBinary(tau, true, n, true, RevBayesCore::DolloAscertainmentBias::Coding(cd), norm);

    }

    // set the root frequencies (by default these are NULL so this is OK)
    dist->setRootFrequencies( rf );

    if ( rate->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* clockRates = static_cast<const ModelVector<RealPos> &>( rate->getRevObject() ).getDagNode();

        // sanity check
        if ( nNodes != clockRates->getValue().size() )
        {
            throw RbException( "The number of clock rates does not match the number of nodes" );
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
        // sanity check
        if ( rootFrequencies == NULL || rootFrequencies->getRevObject() == RevNullObject::getInstance() )
        {
            throw RbException( "If you provide branch-heterogeneous substitution matrices, then you also need to provide root frequencies." );
        }


        dist->setRateMatrix( rm );
    }
    else if ( q->getRevObject() != RevNullObject::getInstance() )
    {
        RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rm = static_cast<const RateGenerator &>( q->getRevObject() ).getDagNode();
        dist->setRateMatrix( rm );
    }
    else
    {

    }

    if ( siteRatesNode != NULL && siteRatesNode->getValue().size() > 0 )
    {
        dist->setSiteRates( siteRatesNode );
    }

    dist->setDeathRate(deathRateNode);

    d = dist;


    return d;
}



/* Get Rev type of object */
const std::string& Dist_phyloCTMCDollo::getClassType(void)
{

    static std::string rev_type = "Dist_phyloCTMCDollo";

    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_phyloCTMCDollo::getClassTypeSpec(void)
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
std::string Dist_phyloCTMCDollo::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhyloCTMCDollo";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_phyloCTMCDollo::getParameterRules(void) const
{

    static MemberRules         dist_member_rules;
    static bool rulesSet = false;

    if ( !rulesSet )
    {
        dist_member_rules.push_back( new ArgumentRule( "tree", Tree::getClassTypeSpec(), "The tree along which the process evolves.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        std::vector<TypeSpec> rateMatrixTypes;
        rateMatrixTypes.push_back( RateGenerator::getClassTypeSpec() );
        rateMatrixTypes.push_back( ModelVector<RateGenerator>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "Q", rateMatrixTypes, "The global or branch-specific rate matrices.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY , NULL ) );

        dist_member_rules.push_back( new ArgumentRule( "deathRate", RealPos::getClassTypeSpec(), "The Dollo character death rate", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );

        // optional argument for the root frequencies
        dist_member_rules.push_back( new ArgumentRule( "rootFrequencies", Simplex::getClassTypeSpec(), "The root specific frequencies of the characters, if applicable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

        std::vector<TypeSpec> branchRateTypes;
        branchRateTypes.push_back( RealPos::getClassTypeSpec() );
        branchRateTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "branchRates", branchRateTypes, "The global or branch-specific rate multipliers.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );

        ModelVector<RealPos> *defaultSiteRates = new ModelVector<RealPos>();
        dist_member_rules.push_back( new ArgumentRule( "siteRates", ModelVector<RealPos>::getClassTypeSpec(), "The rate categories for the sites.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, defaultSiteRates ) );

        dist_member_rules.push_back( new ArgumentRule( "nSites", Natural::getClassTypeSpec(), "The number of sites, used for simulation.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(0L) ) );
        std::vector<std::string> options;
        options.push_back( "Standard" );
        options.push_back( "Binary" );
        options.push_back( "Restriction" );
        dist_member_rules.push_back( new OptionRule( "type", new RlString("Standard"), options, "The data type, used for simulation and initialization." ) );

        dist_member_rules.push_back( new ArgumentRule( "normalize", RlBoolean::getClassTypeSpec(), "Should we normalize the Dollo substitution matrices?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( true ) ) );

        dist_member_rules.push_back( new ArgumentRule("coding", RlString::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("noabsencesites") ) );

        rulesSet = true;
    }

    return dist_member_rules;
}


const TypeSpec& Dist_phyloCTMCDollo::getTypeSpec( void ) const
{

    static TypeSpec ts = getClassTypeSpec();

    return ts;
}


/** Print value for user */
void Dist_phyloCTMCDollo::printValue(std::ostream& o) const
{

    o << "Dollo-Character-State-Evolution-Along-Tree Process(tree=";
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
    o << ", deathRate=";
    if ( deathRate != NULL ) {
        o << deathRate->getName();
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
    o << ", nSites=";
    if ( nSites != NULL ) {
        o << nSites->getName();
    } else {
        o << "?";
    }
    o << ", type=";
    if ( type != NULL ) {
        o << type->getName();
    } else {
        o << "?";
    }
    o << ", coding=";
    if ( coding != NULL ) {
        o << coding->getName();
    } else {
        o << "?";
    }
    o << ", normalize=";
    if ( normalize != NULL ) {
        o << normalize->getName();
    } else {
        o << "?";
    }
    o << ")";

}


/** Set a member variable */
void Dist_phyloCTMCDollo::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "Q" )
    {
        q = var;
    }
    else if ( name == "deathRate" )
    {
        deathRate = var;
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
    else if ( name == "nSites" )
    {
        nSites = var;
    }
    else if ( name == "type" )
    {
        type = var;
    }
    else if ( name == "normalize" )
    {
        normalize = var;
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


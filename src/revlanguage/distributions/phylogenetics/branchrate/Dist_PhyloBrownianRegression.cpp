#include "Dist_PhyloBrownianRegression.h"

#include <stddef.h>
#include <ostream>

#include "PhyloBrownianRegressionProcess.h"
#include "RlTree.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RbException.h"
#include "RbVector.h"
#include "RealPos.h"
#include "RlDistribution.h"
#include "StringUtilities.h"
#include "Tree.h"
#include "TypeSpec.h"

using namespace RevLanguage;

Dist_PhyloBrownianRegression::Dist_PhyloBrownianRegression() : TypedDistribution< ContinuousCharacterData >()
{
    
}


Dist_PhyloBrownianRegression::~Dist_PhyloBrownianRegression()
{
    
}



Dist_PhyloBrownianRegression* Dist_PhyloBrownianRegression::clone( void ) const
{
    
    return new Dist_PhyloBrownianRegression(*this);
}


RevBayesCore::TypedDistribution< RevBayesCore::ContinuousCharacterData >* Dist_PhyloBrownianRegression::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const Tree &>( tree->getRevObject() ).getDagNode();
    size_t n_nodes = tau->getValue().getNumberOfNodes();
    
    const RevBayesCore::ContinuousCharacterData& X = static_cast<const ContinuousCharacterData &>( predictors->getRevObject() ).getValue();
    
    RevBayesCore::PhyloBrownianRegressionProcess *dist = new RevBayesCore::PhyloBrownianRegressionProcess(tau, X);

    // set the clock rates
    if ( branch_rates->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* br = static_cast<const ModelVector<RealPos> &>( branch_rates->getRevObject() ).getDagNode();
        
        // sanity check
        size_t n_rates = br->getValue().size();
        if ( (n_nodes-1) != n_rates )
        {
            throw RbException( "The number of clock rates does not match the number of branches" );
        }
        
        dist->setBranchRate( br );
        
    }
    else
    {
        RevBayesCore::TypedDagNode< double >* br = static_cast<const RealPos &>( branch_rates->getRevObject() ).getDagNode();
        
        dist->setBranchRate( br );
    }
    
    // set the slopes
    if ( slopes->getRevObject().isType( ModelVector<Real>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* sl = static_cast<const ModelVector<Real> &>( slopes->getRevObject() ).getDagNode();
        dist->setSlope( sl );
    }
    else
    {
        RevBayesCore::TypedDagNode< double >* sl = static_cast<const RealPos &>( slopes->getRevObject() ).getDagNode();
        dist->setSlope( sl );
    }
    
    // set the predictor means
    if ( means->getRevObject().isType( ModelVector<Real>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* m = static_cast<const ModelVector<Real> &>( means->getRevObject() ).getDagNode();
        dist->setMeanPredictor( m );
    }
    else
    {
        RevBayesCore::TypedDagNode< double >* m = static_cast<const RealPos &>( means->getRevObject() ).getDagNode();
        dist->setMeanPredictor( m );
    }
    
    return dist;
}



/* Get Rev type of object */
const std::string& Dist_PhyloBrownianRegression::getClassType(void)
{
    
    static std::string rev_type = "Dist_PhyloBrownianRegression";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_PhyloBrownianRegression::getClassTypeSpec(void)
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
std::string Dist_PhyloBrownianRegression::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhyloBrownianRegression";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_PhyloBrownianRegression::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        dist_member_rules.push_back( new ArgumentRule( "tree" , Tree::getClassTypeSpec(), "The tree along which the process evolves.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        std::vector<TypeSpec> branchRateTypes;
        branchRateTypes.push_back( RealPos::getClassTypeSpec() );
        branchRateTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "branchRates" , branchRateTypes, "The per branch rate-multiplier(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        
        std::vector<TypeSpec> slopes_types;
        slopes_types.push_back( Real::getClassTypeSpec() );
        slopes_types.push_back( ModelVector<Real>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "slopes" , slopes_types, "The regression slope(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        std::vector<TypeSpec> means_types;
        means_types.push_back( Real::getClassTypeSpec() );
        means_types.push_back( ModelVector<Real>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "means" , means_types, "The regression predictor mean(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        dist_member_rules.push_back( new ArgumentRule( "predictors" , ContinuousCharacterData::getClassTypeSpec(), "The continuous predictor variables.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_PhyloBrownianRegression::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_PhyloBrownianRegression::printValue(std::ostream& o) const
{
        
    o << "PhyloBrownianProcess(tree=";
    if ( tree != NULL )
    {
        o << tree->getName();
    }
    else
    {
        o << "?";
    }
    o << ", branchRates=";
    if ( branch_rates != NULL )
    {
        o << branch_rates->getName();
    }
    else
    {
        o << "?";
    }
    o << ", slopes=";
    if ( slopes != NULL )
    {
        o << slopes->getName();
    }
    else
    {
        o << "?";
    }
    o << ", means=";
    if ( means != NULL )
    {
        o << means->getName();
    }
    else
    {
        o << "?";
    }
    o << ", predictor=";
    if ( predictors != NULL )
    {
        o << predictors->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
    
}


/** Set a member variable */
void Dist_PhyloBrownianRegression::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "branchRates" )
    {
        branch_rates = var;
    }
    else if ( name == "slopes" )
    {
        slopes = var;
    }
    else if ( name == "means" )
    {
        means = var;
    }
    else if ( name == "predictors" )
    {
        predictors = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
    
}


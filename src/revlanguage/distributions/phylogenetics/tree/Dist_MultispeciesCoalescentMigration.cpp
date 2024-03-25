#include <math.h>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_MultispeciesCoalescentMigration.h"
#include "ModelVector.h"
#include "MultispeciesCoalescentMigration.h"
#include "RealPos.h"
#include "RlTaxon.h"
#include "RlTimeTree.h"
#include "StochasticNode.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "ModelObject.h"
#include "RbException.h"
#include "RbHelpReference.h"
#include "RbVector.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlRateGenerator.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

namespace RevBayesCore { class Taxon; }

using namespace RevLanguage;


/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_MultispeciesCoalescentMigration::Dist_MultispeciesCoalescentMigration() : TypedDistribution<TimeTree>()
{

}

/**
 * Clone the object
 *
 * \return a clone of the object.
 */

Dist_MultispeciesCoalescentMigration* Dist_MultispeciesCoalescentMigration::clone(void) const
{

    return new Dist_MultispeciesCoalescentMigration(*this);

}




/**
 * Create a new internal distribution object.
 *
 * This function simply dynamically allocates a new internal distribution object that can be
 * associated with the variable. The internal distribution object is created by calling its
 * constructor and passing the distribution-parameters (other DAG nodes) as arguments of the
 * constructor. The distribution constructor takes care of the proper hook-ups.
 *
 * \return A new internal distribution object.
 */
RevBayesCore::MultispeciesCoalescentMigration* Dist_MultispeciesCoalescentMigration::createDistribution( void ) const
{

    // Get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* st = static_cast<const TimeTree &>( species_tree->getRevObject() ).getDagNode();
    const std::vector<RevBayesCore::Taxon>      &t  = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();

    // get the number of nodes for the tree
    size_t n_nodes = st->getValue().getNumberOfNodes();

    RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>*    q = static_cast<const RateGenerator &   >( Q->getRevObject()        ).getDagNode();
    RevBayesCore::TypedDagNode<double>*                         d = static_cast<const RealPos &         >( delta->getRevObject()    ).getDagNode();

    RevBayesCore::MultispeciesCoalescentMigration*   dist = new RevBayesCore::MultispeciesCoalescentMigration( st, t, q, d );

    if ( Ne->getRequiredTypeSpec().isDerivedOf( ModelVector<RealPos>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* ne_node = static_cast<const ModelVector<RealPos> &>( Ne->getRevObject() ).getDagNode();

        // sanity check
        if ( n_nodes != ne_node->getValue().size() )
        {
            throw RbException( "The number of effective population sizes does not match the number of branches." );
        }

        dist->setNes( ne_node );
    }
    else
    {
        RevBayesCore::TypedDagNode<double>* ne_node = static_cast<const RealPos &>( Ne->getRevObject() ).getDagNode();
        dist->setNe( ne_node );
    }
    dist->redrawValue();

    return dist;
}



/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_MultispeciesCoalescentMigration::getClassType(void)
{

    static std::string rev_type = "Dist_MultispeciesCoalescentMigration";

    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_MultispeciesCoalescentMigration::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<TimeTree>::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_MultispeciesCoalescentMigration::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "MultispeciesCoalescentMigration";

    return d_name;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_MultispeciesCoalescentMigration::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
      a_names.push_back( "MSCM" );

    return a_names;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the Multispecies Coalescent process are:
 * (1) Species tree.
 * (2) Population size.
 * (3) Gene name to species name map.
 * (4) the number of taxa.
 *
 * \return The member rules.
 */
const MemberRules& Dist_MultispeciesCoalescentMigration::getParameterRules(void) const
{

    static MemberRules member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
        member_rules.push_back( new ArgumentRule( "speciesTree", TimeTree::getClassTypeSpec(), "The species tree in which the gene trees evolve.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        std::vector<TypeSpec> branch_Ne_types;
        branch_Ne_types.push_back( RealPos::getClassTypeSpec() );
        branch_Ne_types.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        member_rules.push_back( new ArgumentRule( "Ne",     branch_Ne_types, "The effective population size(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule( "Q",      RateGenerator::getClassTypeSpec(), "The migration rate matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule( "rate",   RealPos::getClassTypeSpec(), "The rate scalar of the migration rate matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule( "taxa",   ModelVector<Taxon>::getClassTypeSpec(), "The vector of taxa which have species and individual names.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return member_rules;
}


/**
 * Set a member variable.
 *
 * Sets a member variable with the given name and store the pointer to the variable.
 * The value of the variable might still change but this function needs to be called again if the pointer to
 * the variable changes. The current values will be used to create the distribution object.
 *
 * \param[in]    name     Name of the member variable.
 * \param[in]    var      Pointer to the variable.
 */
void Dist_MultispeciesCoalescentMigration::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "speciesTree" )
    {
        species_tree = var;
    }
    else if ( name == "Ne" )
    {
        Ne = var;
    }
    else if ( name == "Q" )
    {
        Q = var;
    }
    else if ( name == "rate" )
    {
        delta = var;
    }
    else if ( name == "taxa" )
    {
        taxa = var;
    }
    else
    {
        TypedDistribution<TimeTree>::setConstParameter(name, var);
    }

}


/** Get type spec */
const TypeSpec& Dist_MultispeciesCoalescentMigration::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}

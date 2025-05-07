#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Move_NarrowExchangeRateMatrix.h"
#include "NarrowExchangeRateMatrixProposal.h"
#include "RbException.h"
#include "RealPos.h"
#include "RlRateGenerator.h"
#include "RlTimeTree.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "DagNode.h"
#include "ModelObject.h"
#include "Move.h"
#include "RbVector.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"

namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { class RateGenerator; }
namespace RevBayesCore { class Tree; }


using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Move_NarrowExchangeRateMatrix::Move_NarrowExchangeRateMatrix() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the move.
 */
Move_NarrowExchangeRateMatrix* Move_NarrowExchangeRateMatrix::clone(void) const
{
    
    return new Move_NarrowExchangeRateMatrix(*this);
}


/**
 * Create a new internal move object.
 *
 * This function simply dynamically allocates a new internal move object that is
 * associated with the variable (DAG-node). The internal move object is created by calling its
 * constructor and passing the move-parameters (the variable and other parameters) as arguments of the
 * constructor. The move constructor takes care of the proper hook-ups.
 *
 * \return A new internal distribution object.
 */
void Move_NarrowExchangeRateMatrix::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tmp = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );
//    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* rm = static_cast<const ModelVector<RateGenerator> &>( rate_matrices->getRevObject() ).getDagNode();

    std::vector< RevBayesCore::StochasticNode<RevBayesCore::RateGenerator> *> rm;
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::RateGenerator> >* deterministic_vector = static_cast<const ModelVector<RateGenerator> &>( rate_matrices->getRevObject() ).getDagNode();
    std::vector<const RevBayesCore::DagNode*> parents = deterministic_vector->getParents();
    for (std::vector<const RevBayesCore::DagNode*>::const_iterator it = parents.begin(); it != parents.end(); ++it)
    {
        const RevBayesCore::DagNode *tmp_node = *it;
        const RevBayesCore::StochasticNode<RevBayesCore::RateGenerator> *the_node = dynamic_cast< const RevBayesCore::StochasticNode<RevBayesCore::RateGenerator>* >( tmp_node );
        if ( the_node != NULL )
        {
            rm.push_back( const_cast<RevBayesCore::StochasticNode<RevBayesCore::RateGenerator> *>( the_node ) );
        }
        else
        {
            throw RbException("Could not add the node because it isn't a vector of stochastic nodes of rate matrices.");
        }
    }
    
    RevBayesCore::Proposal *p = new RevBayesCore::NarrowExchangeRateMatrixProposal(n,rm);
    value = new RevBayesCore::MetropolisHastingsMove(p,w);
    
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Move_NarrowExchangeRateMatrix::getClassType(void)
{
    
    static std::string rev_type = "Move_NarrowExchangeRateMatrix";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Move_NarrowExchangeRateMatrix::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_NarrowExchangeRateMatrix::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "NarrowExchangeRateMatrix";
    
    return c_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the scale move are:
 * (1) the variable which must be a time-tree.
 *
 * \return The member rules.
 */
const MemberRules& Move_NarrowExchangeRateMatrix::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        memberRules.push_back( new ArgumentRule( "tree", TimeTree::getClassTypeSpec(), "The tree variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        memberRules.push_back( new ArgumentRule( "rateMatrices", ModelVector<RateGenerator>::getClassTypeSpec(), "The branch-specific rate matrices.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::DETERMINISTIC ) );
        

        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        memberRules.insert( memberRules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return memberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Move_NarrowExchangeRateMatrix::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/**
 * Print the value for the user.
 */
void Move_NarrowExchangeRateMatrix::printValue(std::ostream &o) const
{
    
    o << "Narrow(";
    if (tree != NULL)
    {
        o << tree->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
    
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
void Move_NarrowExchangeRateMatrix::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "rateMatrices" )
    {
        rate_matrices = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
    
}





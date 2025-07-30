#include "Move_IndependentTopology.h"

#include <cstddef>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "IndependentTopologyProposal.h"
#include "ModelVector.h"
#include "MetropolisHastingsMove.h"
#include "RbException.h"
#include "RealPos.h"
#include "RlClade.h"
#include "RlTimeTree.h"
#include "RlTypedDistribution.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "DagNode.h"
#include "ModelObject.h"
#include "Move.h"
#include "RevNullObject.h"
#include "RlDistribution.h"
#include "RlTree.h"
#include "StochasticNode.h"
#include "StringUtilities.h"
#include "TypedDistribution.h"

namespace RevBayesCore { class Clade; }
namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class RbVector; }


using namespace RevLanguage;


Move_IndependentTopology::Move_IndependentTopology() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_IndependentTopology* Move_IndependentTopology::clone(void) const
{
    
    return new Move_IndependentTopology(*this);
}


void Move_IndependentTopology::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new empirical tree move
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    size_t del = static_cast<const Natural &>( delay->getRevObject() ).getValue();
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tmp = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );

    const Distribution& rlDistribution                        = static_cast<const Distribution &>( proposal_distribution->getRevObject() );
    RevBayesCore::TypedDistribution<RevBayesCore::Tree>* prop = static_cast<RevBayesCore::TypedDistribution<RevBayesCore::Tree>* >( rlDistribution.createDistribution() );

    const RevBayesCore::Clade& c      = static_cast<const Clade&>( outgroup->getRevObject() ).getValue();

    RevBayesCore::Proposal *p;

    if( substitution_rates->getRevObject() != RevNullObject::getInstance() )
    {
        RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* tmp = static_cast<const ModelVector<RealPos> &>( substitution_rates->getRevObject() ).getDagNode();
        if( tmp->isStochastic() == true )
        {
            RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> >* stoch = dynamic_cast<RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> >* >( tmp );

            p = new RevBayesCore::IndependentTopologyProposal(t, prop, stoch, c);
        }
        else
        {

            std::vector<const RevBayesCore::DagNode*> parents = tmp->getParents();
            std::vector< RevBayesCore::StochasticNode<double> *> n;
            for (std::vector<const RevBayesCore::DagNode*>::const_iterator it = parents.begin(); it != parents.end(); ++it)
            {
                const RevBayesCore::StochasticNode<double> *the_node = dynamic_cast< const RevBayesCore::StochasticNode<double>* >( *it );
                if ( the_node != NULL )
                {
                    n.push_back( const_cast< RevBayesCore::StochasticNode<double>* >( the_node ) );
                }
                else
                {
                    throw RbException("Substitution rates vector isn't stochastic or a vector of stochastic nodes.");
                }
            }

            p = new RevBayesCore::IndependentTopologyProposal(t, prop, n, c);
        }

    }
    else
    {
        p = new RevBayesCore::IndependentTopologyProposal(t, prop, NULL, c);
    }

    value = new RevBayesCore::MetropolisHastingsMove(p, w, del, false);
    
}


/** Get Rev type of object */
const std::string& Move_IndependentTopology::getClassType(void)
{
    static std::string rev_type = "Move_IndependentTopology";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_IndependentTopology::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_IndependentTopology::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "IndependentTopology";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_IndependentTopology::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        move_member_rules.push_back( new ArgumentRule( "tree", TimeTree::getClassTypeSpec(), "The stochastic time tree variable on which this moves operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "proposal", TypedDistribution<Tree>::getClassTypeSpec(), "The proposal topology distribution.",  ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        move_member_rules.push_back( new ArgumentRule( "rates", ModelVector<RealPos>::getClassTypeSpec(), "The vector of branch-specific substitution rates.",  ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, NULL ) );
        move_member_rules.push_back( new ArgumentRule( "outgroup", Clade::getClassTypeSpec(), "Outgroup used to root the proposal topologies.",  ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Clade() ) );

        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_IndependentTopology::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_IndependentTopology::printValue(std::ostream &o) const
{
    
    o << "Move_IndependentTopology(";
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


/** Set a member variable */
void Move_IndependentTopology::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if( name == "proposal" )
    {
        proposal_distribution = var;
    }
    else if( name == "rates" )
    {
        substitution_rates = var;
    }
    else if( name == "outgroup" )
    {
        outgroup = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}


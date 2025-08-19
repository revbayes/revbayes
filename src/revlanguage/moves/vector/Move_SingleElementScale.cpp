#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Move_SingleElementScale.h"
#include "RealPos.h"
#include "SingleElementScaleProposal.h"
#include "DagNode.h"
#include "ModelObject.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RbException.h"
#include "Real.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlMove.h"
#include "StochasticNode.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevLanguage;

Move_SingleElementScale::Move_SingleElementScale() : Move () {}

Move_SingleElementScale* Move_SingleElementScale::clone(void) const {

    return new Move_SingleElementScale(*this);
}


void Move_SingleElementScale::constructInternalObject( void )
{
    // we free the memory first
    delete value;

    // now allocate a new element scaling move
    double l = static_cast<const RealPos &>( lambda->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* tmp = static_cast<const ModelVector<RealPos> &>( x->getRevObject() ).getDagNode();
    std::vector<const RevBayesCore::DagNode*> p = tmp->getParents();
    std::vector< RevBayesCore::StochasticNode<double> *> n;
    for (std::vector<const RevBayesCore::DagNode*>::const_iterator it = p.begin(); it != p.end(); ++it)
    {
        const RevBayesCore::StochasticNode<double> *the_node = dynamic_cast< const RevBayesCore::StochasticNode<double>* >( *it );
        if ( the_node != NULL )
        {
            n.push_back( const_cast< RevBayesCore::StochasticNode<double>* >( the_node ) );
        }
        else
        {
            throw RbException("Could not create a mvVectorElementScale because the node isn't a vector of stochastic nodes.");
        }
    }

    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();

    RevBayesCore::Proposal *prop = new RevBayesCore::SingleElementScaleProposal(n,l);
    value = new RevBayesCore::MetropolisHastingsMove(prop,w,t);
}


const std::string& Move_SingleElementScale::getClassType(void) {

    static std::string rev_type = "Move_SingleElementScale";

    return rev_type;
}


const TypeSpec& Move_SingleElementScale::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );

    return rev_type_spec;
}


std::string Move_SingleElementScale::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "VectorSingleElementScale";

    return c_name;
}

/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the SingleElementScale move are:
 * (1) the variable which must be a vector of positive real nodes.
 * (2) the tuning parameter lambda that defines the size of the proposal (positive real)
 * (3) a flag whether auto-tuning should be used.
 *
 * \return The member rules.
 */
const MemberRules& Move_SingleElementScale::getParameterRules(void) const
{

    static MemberRules move_member_rules;
    static bool rules_set = false;

    if ( !rules_set )
    {

        move_member_rules.push_back( new ArgumentRule( "x"     , ModelVector<RealPos>::getClassTypeSpec(), "The variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::DETERMINISTIC ) );
        move_member_rules.push_back( new ArgumentRule( "lambda", RealPos::getClassTypeSpec()             , "The scaling factor (strength) of this move.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new Real(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune"  , RlBoolean::getClassTypeSpec()           , "Should we tune the scaling factor during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );

        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );

        rules_set = true;
    }

    return move_member_rules;
}


const TypeSpec& Move_SingleElementScale::getTypeSpec( void ) const {

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


void Move_SingleElementScale::printValue(std::ostream &o) const {

    o << "SingleElementScale(";
    if (x != NULL)
    {
        o << x->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


void Move_SingleElementScale::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {

    if ( name == "x" )
    {
        x = var;
    }
    else if ( name == "lambda" )
    {
        lambda = var;
    }
    else if ( name == "tune" )
    {
        tune = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}

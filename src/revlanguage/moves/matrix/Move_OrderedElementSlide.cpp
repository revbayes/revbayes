#include "Move_OrderedElementSlide.h"

#include <cstddef>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "OrderedElementSlideProposal.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlMatrixReal.h"
#include "TypeSpec.h"
#include "Move.h"
#include "StochasticNode.h"

namespace RevBayesCore { class MatrixReal; }
namespace RevBayesCore { class Proposal; }

using namespace RevLanguage;

Move_OrderedElementSlide::Move_OrderedElementSlide() : Move()
{
}


Move_OrderedElementSlide* Move_OrderedElementSlide::clone(void) const
{
    return new Move_OrderedElementSlide(*this);
}


void Move_OrderedElementSlide::constructInternalObject( void )
{
    delete value;

    RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal>* tmp = static_cast<const MatrixReal &>( x->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::MatrixReal> *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::MatrixReal> *>( tmp );

    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();

    RevBayesCore::Proposal *p = new RevBayesCore::OrderedElementSlideProposal( n );
    value = new RevBayesCore::MetropolisHastingsMove(p, w, false);
}


const std::string& Move_OrderedElementSlide::getClassType(void)
{
    static std::string rev_type = "Move_OrderedElementSlide";

    return rev_type;
}


const TypeSpec& Move_OrderedElementSlide::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );

    return rev_type_spec;
}


std::string Move_OrderedElementSlide::getMoveName( void ) const
{
    static std::string c_name = "OrderedElementSlide";

    return c_name;
}


const MemberRules& Move_OrderedElementSlide::getParameterRules(void) const
{
    static MemberRules move_member_rules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        move_member_rules.push_back( new ArgumentRule( "x", MatrixReal::getClassTypeSpec(), "The variable whose parts are ordered vectors.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );

        const MemberRules& inherited_rules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inherited_rules.begin(), inherited_rules.end() );

        rules_set = true;
    }

    return move_member_rules;
}


const TypeSpec& Move_OrderedElementSlide::getTypeSpec( void ) const
{
    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


void Move_OrderedElementSlide::printValue(std::ostream &o) const
{
    o << "OrderedElementSlide(";
    if (x != NULL) { o << x->getName(); } else { o << "?"; }
    o << ")";
}


void Move_OrderedElementSlide::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "x" )
    {
        x = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}

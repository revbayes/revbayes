#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Move_MatrixSingleElementSwap.h"
#include "MatrixRealSingleElementSwapProposal.h"
#include "MatrixReal.h"
#include "Natural.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RbException.h"
#include "RevNullObject.h"
#include "RlMatrixReal.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"


using namespace RevLanguage;

Move_MatrixSingleElementSwap::Move_MatrixSingleElementSwap() : Move()
{

}


Move_MatrixSingleElementSwap* Move_MatrixSingleElementSwap::clone(void) const
{

    return new Move_MatrixSingleElementSwap(*this);
}


void Move_MatrixSingleElementSwap::constructInternalObject( void )
{
    // we free the memory first
    delete value;

    RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal>* tmp = static_cast<const MatrixReal &>( x->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::MatrixReal> *m = static_cast<RevBayesCore::StochasticNode<RevBayesCore::MatrixReal> *>( tmp );

    double we = static_cast<const RealPos &>( weight->getRevObject() ).getValue();

    bool has_margin = ( margin->getRevObject() != RevNullObject::getInstance() );
    bool has_row    = ( row->getRevObject() != RevNullObject::getInstance() );
    bool has_col    = ( col->getRevObject() != RevNullObject::getInstance() );

    if ( has_margin + has_row + has_col > 1 )
    {
        throw RbException("mvMatrixElementSwap takes a margin, a row or a column, not more than one.");
    }

    // row= and col= pin a line and imply the margin; Rev indices are 1-based, and -1 draws a line.
    // With none of them the whole matrix is the pool.
    std::int64_t mg = 0;
    std::int64_t i  = -1;

    if ( has_margin == true )
    {
        mg = static_cast<const Natural &>( margin->getRevObject() ).getValue();

        if ( mg != 1 && mg != 2 )
        {
            throw RbException("mvMatrixElementSwap takes margin 1 (within a row) or 2 (within a column).");
        }
    }
    else if ( has_row == true )
    {
        mg = 1;
        i = static_cast<const Natural &>( row->getRevObject() ).getValue() - 1;
    }
    else if ( has_col == true )
    {
        mg = 2;
        i = static_cast<const Natural &>( col->getRevObject() ).getValue() - 1;
    }

    RevBayesCore::MatrixRealSingleElementSwapProposal *p = new RevBayesCore::MatrixRealSingleElementSwapProposal( m, mg, i );

    value = new RevBayesCore::MetropolisHastingsMove(p, we, false);
}


/** Get Rev type of object */
const std::string& Move_MatrixSingleElementSwap::getClassType(void)
{

    static std::string rev_type = "Move_MatrixSingleElementSwap";

    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Move_MatrixSingleElementSwap::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 */
std::string Move_MatrixSingleElementSwap::getMoveName( void ) const
{
    std::string c_name = "MatrixElementSwap";

    return c_name;
}


/** Return member rules */
const MemberRules& Move_MatrixSingleElementSwap::getParameterRules(void) const
{

    static MemberRules move_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
        move_member_rules.push_back( new ArgumentRule( "x", MatrixReal::getClassTypeSpec(), "The matrix on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "margin", Natural::getClassTypeSpec(), "Swap within a row (1) or within a column (2), as in R's MARGIN. Omit to swap any two elements.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
        move_member_rules.push_back( new ArgumentRule( "row", Natural::getClassTypeSpec(), "Swap within this row rather than a random line.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
        move_member_rules.push_back( new ArgumentRule( "col", Natural::getClassTypeSpec(), "Swap within this column rather than a random line.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );

        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );

        rules_set = true;
    }

    return move_member_rules;
}


/** Get type spec */
const TypeSpec& Move_MatrixSingleElementSwap::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


/** Get type spec */
void Move_MatrixSingleElementSwap::printValue(std::ostream &o) const
{

    o << "Move_MatrixSingleElementSwap(";
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


/** Set a member variable */
void Move_MatrixSingleElementSwap::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "x" )
    {
        x = var;
    }
    else if ( name == "margin" )
    {
        margin = var;
    }
    else if ( name == "row" )
    {
        row = var;
    }
    else if ( name == "col" )
    {
        col = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}

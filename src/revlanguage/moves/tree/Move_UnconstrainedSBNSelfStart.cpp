#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "SubtreePruneRegraftProposal.h"
#include "NearestNeighborInterchange_nonClockProposal.h"
#include "UnrootedUnconstrainedSBNProposal.h"
#include "RevObject.h"
#include "RbException.h"
#include "RealPos.h"
#include "ModelVector.h"
#include "UnconstrainedSBNMetaMove.h"
#include "Move_UnconstrainedSBNSelfStart.h"
#include "RlBranchLengthTree.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"


using namespace RevLanguage;

Move_UnconstrainedSBNSelfStart::Move_UnconstrainedSBNSelfStart() : Move()
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_UnconstrainedSBNSelfStart* Move_UnconstrainedSBNSelfStart::clone(void) const
{

	return new Move_UnconstrainedSBNSelfStart(*this);
}


void Move_UnconstrainedSBNSelfStart::constructInternalObject( void )
{
    // we free the memory first
    delete value;

    // Get variable
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const BranchLengthTree &>( tree->getRevObject() ).getDagNode();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );

    // Get options
    // RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* rmw = static_cast<const ModelVector<RealPos> &>( relativeMoveWeights->getRevObject() ).getDagNode();
    const RevBayesCore::RbVector<double> &rmw = static_cast<const ModelVector<RealPos> &>( relativeMoveWeights->getRevObject() ).getValue();

    // const RevBayesCore::RbVector<long> &e = static_cast<const ModelVector<Natural> &>( inidices->getRevObject() ).getValue();

    int n0    = static_cast<const Natural &>( burnin->getRevObject() ).getValue();
    int c0    = static_cast<const Natural &>( waitBeforeFitting->getRevObject() ).getValue();

    // relativeMoveWeights( m ),
    // waitBeforeRecording( n0 ),
    // waitBeforeFitting( c0 ),


    // The proposals
    // RevBayesCore::SBNParameters sbn = RevBayesCore::SBNParameters();
    RevBayesCore::UnrootedUnconstrainedSBNProposal *u = new RevBayesCore::UnrootedUnconstrainedSBNProposal(t);
std::cout << "made UnrootedUnconstrainedSBNProposal object" << std::endl;

    RevBayesCore::Proposal *s = new RevBayesCore::SubtreePruneRegraftProposal(t);

    RevBayesCore::Proposal *n = new RevBayesCore::NearestNeighborInterchange_nonClockProposal(t);

    // TODO: need to get rid of tuning, this move must tune
    bool tune = true;

    // Bundle proposals
    value = new RevBayesCore::UnconstrainedSBNMetaMove(u,n,s,rmw,n0,c0,w,tune);

}


/** Get class name of object */
const std::string& Move_UnconstrainedSBNSelfStart::getClassName(void)
{

    static std::string rbClassName = "Move_UnconstrainedSBNSelfStart";

	return rbClassName;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_UnconstrainedSBNSelfStart::getClassTypeSpec(void)
{

    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Move::getClassTypeSpec() ) );

	return rbClass;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_UnconstrainedSBNSelfStart::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "UnconstrainedSBNSelfStart";

    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_UnconstrainedSBNSelfStart::getParameterRules(void) const
{

    static MemberRules USBNSSMemberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {

        USBNSSMemberRules.push_back( new ArgumentRule( "tree",                BranchLengthTree::getClassTypeSpec(), "The tree variable this move operates on.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        USBNSSMemberRules.push_back( new ArgumentRule( "relativeMoveWeights", ModelVector<RealPos>::getClassTypeSpec(), "The relative weights on the SBN, NNI, and SPR moves respectively.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        USBNSSMemberRules.push_back( new ArgumentRule( "burnin",              Natural::getClassTypeSpec()  , "The number of move attempts to wait before tracking the tree.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(5000) ) );
        USBNSSMemberRules.push_back( new ArgumentRule( "waitBeforeFitting",   Natural::getClassTypeSpec()  , "The number of move attempts post-burnin to wait before fitting the SBN.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(10000) ) );

        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        USBNSSMemberRules.insert( USBNSSMemberRules.end(), inheritedRules.begin(), inheritedRules.end() );

        rules_set = true;
    }

    return USBNSSMemberRules;
}

/** Get type spec */
const TypeSpec& Move_UnconstrainedSBNSelfStart::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}



/** Get type spec */
void Move_UnconstrainedSBNSelfStart::printValue(std::ostream &o) const
{

    o << "UnconstrainedSBNSelfStart(";
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


/** Set a NearestNeighborInterchange variable */
void Move_UnconstrainedSBNSelfStart::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "relativeMoveWeights" )
    {
        relativeMoveWeights = var;
    }
    else if ( name == "burnin" )
    {
        burnin = var;
    }
    else if ( name == "waitBeforeFitting" )
    {
        waitBeforeFitting = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}

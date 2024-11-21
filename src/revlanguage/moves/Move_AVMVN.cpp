#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Move_AVMVN.h"
#include "AVMVNProposal.h"
#include "Probability.h"
#include "RbException.h"
#include "RbMathLogic.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "ContinuousStochasticNode.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "Move.h"
#include "Natural.h"
#include "RbBoolean.h"
#include "Real.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "RlUtils.h"
#include "StochasticNode.h"

namespace RevBayesCore { class Simplex; }
namespace RevBayesCore { template <class valueType> class RbVector; }


using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Move_AVMVN::Move_AVMVN() : Move()
{

    // add member methods

    // first, the argument rules
    ArgumentRules* addRealArgRules                  = new ArgumentRules();
    ArgumentRules* addRealPosArgRules               = new ArgumentRules();
    ArgumentRules* addSimplexArgRules               = new ArgumentRules();
    ArgumentRules* addModelVectorArgRules           = new ArgumentRules();
    ArgumentRules* removeRealArgRules               = new ArgumentRules();
    ArgumentRules* removeRealPosArgRules            = new ArgumentRules();
    ArgumentRules* removeSimplexArgRules            = new ArgumentRules();
    ArgumentRules* removeModelVectorArgRules        = new ArgumentRules();


    // next, set the specific arguments
    addRealArgRules->push_back(                     new ArgumentRule( "var"        , Real::getClassTypeSpec(),                 "The variable to move"             , ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
    addRealPosArgRules->push_back(                  new ArgumentRule( "var"        , RealPos::getClassTypeSpec(),              "The variable to move"             , ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
    addSimplexArgRules->push_back(                  new ArgumentRule( "var"        , Simplex::getClassTypeSpec(),              "The variable to move"             , ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
    addModelVectorArgRules->push_back(              new ArgumentRule( "var"        , ModelVector<Real>::getClassTypeSpec(),    "The variable to move"             , ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
    removeRealArgRules->push_back(                  new ArgumentRule( "var"        , Real::getClassTypeSpec(),                 "The variable to move"             , ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
    removeRealPosArgRules->push_back(               new ArgumentRule( "var"        , RealPos::getClassTypeSpec(),              "The variable to move"             , ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
    removeSimplexArgRules->push_back(               new ArgumentRule( "var"        , Simplex::getClassTypeSpec(),              "The variable to move"             , ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
    removeModelVectorArgRules->push_back(           new ArgumentRule( "var"        , ModelVector<Real>::getClassTypeSpec(),    "The variable to move"             , ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );


    // finally, create the methods
    methods.addFunction( new MemberProcedure( "addVariable",    RlUtils::Void, addRealArgRules) );
    methods.addFunction( new MemberProcedure( "addVariable",    RlUtils::Void, addRealPosArgRules) );
    methods.addFunction( new MemberProcedure( "addVariable",    RlUtils::Void, addSimplexArgRules) );
    methods.addFunction( new MemberProcedure( "addVariable",    RlUtils::Void, addModelVectorArgRules) );
    methods.addFunction( new MemberProcedure( "removeVariable", RlUtils::Void, removeRealArgRules) );
    methods.addFunction( new MemberProcedure( "removeVariable", RlUtils::Void, removeRealPosArgRules) );
    methods.addFunction( new MemberProcedure( "removeVariable", RlUtils::Void, removeSimplexArgRules) );
    methods.addFunction( new MemberProcedure( "removeVariable", RlUtils::Void, removeModelVectorArgRules) );

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
Move_AVMVN* Move_AVMVN::clone(void) const
{

    return new Move_AVMVN(*this);
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
void Move_AVMVN::constructInternalObject( void )
{

    // we free the memory first
    delete value;

    // now allocate a new up-down-scale move
    double s = static_cast<const RealPos &>( sigmaSquared->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    double e = static_cast<const Probability &>( epsilon->getRevObject() ).getValue();
    int n0   = static_cast<const Natural &>( waitBeforeLearning->getRevObject() ).getValue();
    int c0   = static_cast<const Natural &>( waitBeforeUsing->getRevObject() ).getValue();
    int m    = static_cast<const Natural &>( maxUpdates->getRevObject() ).getValue();

    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();

    // finally create the internal move object
    RevBayesCore::AVMVNProposal *prop = new RevBayesCore::AVMVNProposal(s, e, n0, c0, m);

    value = new RevBayesCore::MetropolisHastingsMove(prop,w,t);

}

RevPtr<RevVariable> Move_AVMVN::executeMethod(const std::string& name, const std::vector<Argument>& args, bool &found)
{

    if ( name == "addVariable" )
    {
        found = true;

        Real* uReal = dynamic_cast<Real *>( &args[0].getVariable()->getRevObject() );
        Real* upReal = dynamic_cast<RealPos *>( &args[0].getVariable()->getRevObject() );
        Simplex* sim = dynamic_cast<Simplex *>( &args[0].getVariable()->getRevObject() );
        ModelVector<Real>* uVector = dynamic_cast<ModelVector<Real> *>( &args[0].getVariable()->getRevObject() );

        // Handle scalar variables with possible transforms
        if ( upReal != NULL )
        {
            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::AVMVNProposal &prop = static_cast<RevBayesCore::AVMVNProposal&>( m->getProposal() );

            // Logit-transform or log-transform?
            RevBayesCore::ContinuousStochasticNode *n = dynamic_cast<RevBayesCore::ContinuousStochasticNode *>( upReal->getDagNode() );
            if ( n != NULL && RevBayesCore::RbMath::isFinite(n->getMin()) && RevBayesCore::RbMath::isFinite(n->getMax()) )
            {
                prop.addLogitScalar(n);
            }
            else
            {
                RevBayesCore::StochasticNode<double> *n2 = dynamic_cast<RevBayesCore::StochasticNode<double> *>( upReal->getDagNode() );
                prop.addLogScalar(n2);

            }

        }
        else if ( uReal != NULL )
        {
            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::AVMVNProposal &prop = static_cast<RevBayesCore::AVMVNProposal&>( m->getProposal() );

            // Logit-transform?
            RevBayesCore::ContinuousStochasticNode *n = dynamic_cast<RevBayesCore::ContinuousStochasticNode *>( uReal->getDagNode() );
            if ( n != NULL && RevBayesCore::RbMath::isFinite(n->getMin()) && RevBayesCore::RbMath::isFinite(n->getMax()) )
            {
                prop.addLogitScalar(n);
            }
            else
            {
                RevBayesCore::StochasticNode<double> *n2 = dynamic_cast<RevBayesCore::StochasticNode<double> *>( uReal->getDagNode() );
                prop.addUntransformedScalar(n2);

            }

        }
        else if ( sim != NULL )
        {
            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::AVMVNProposal &prop = static_cast<RevBayesCore::AVMVNProposal&>( m->getProposal() );

            RevBayesCore::StochasticNode<RevBayesCore::Simplex> *n = dynamic_cast<RevBayesCore::StochasticNode<RevBayesCore::Simplex > *>( sim->getDagNode() );
            prop.addLogConstrainedSumVector(n);

        }
        else if ( uVector != NULL && uVector->getDagNode()->isStochastic() == true )
        {
            RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> > *the_node = dynamic_cast< RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> > * >( uVector->getDagNode() );
            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::AVMVNProposal &prop = static_cast<RevBayesCore::AVMVNProposal&>( m->getProposal() );

            if ( the_node != NULL )
            {
                prop.addUntransformedVector( the_node );
            }
            else
            {
                throw RbException("Could not add the node because it isn't a stochastic nodes.");
            }
        }
        else
        {
            throw RbException() << "A problem occured when trying to add " <<  args[0].getVariable()->getName() << " to the move.";
        }

        return NULL;
    }
    else if ( name == "removeVariable" )
    {
        found = true;

        Real* uReal = dynamic_cast<Real *>( &args[0].getVariable()->getRevObject() );
        Real* upReal = dynamic_cast<RealPos *>( &args[0].getVariable()->getRevObject() );
        Simplex* sim = dynamic_cast<Simplex *>( &args[0].getVariable()->getRevObject() );
        ModelVector<Real>* uVector = dynamic_cast<ModelVector<Real> *>( &args[0].getVariable()->getRevObject() );

        // Handle scalar variables with possible transforms
        if ( upReal != NULL )
        {

            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::AVMVNProposal &prop = static_cast<RevBayesCore::AVMVNProposal&>( m->getProposal() );

            // Logit-transform or log-transform?
            RevBayesCore::ContinuousStochasticNode *n = dynamic_cast<RevBayesCore::ContinuousStochasticNode *>( upReal );
            if ( n != NULL && RevBayesCore::RbMath::isFinite(n->getMin()) && RevBayesCore::RbMath::isFinite(n->getMax()) )
            {
                prop.removeLogitScalar(n);
            }
            else
            {
                RevBayesCore::StochasticNode<double> *n2 = dynamic_cast<RevBayesCore::StochasticNode<double> *>( upReal );
                prop.removeLogScalar(n2);

            }

        }
        else if ( uReal != NULL )
        {

            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::AVMVNProposal &prop = static_cast<RevBayesCore::AVMVNProposal&>( m->getProposal() );

            // Logit-transform?
            RevBayesCore::ContinuousStochasticNode *n = dynamic_cast<RevBayesCore::ContinuousStochasticNode *>( uReal );
            if ( n != NULL && RevBayesCore::RbMath::isFinite(n->getMin()) && RevBayesCore::RbMath::isFinite(n->getMax()) )
            {
                prop.removeLogitScalar(n);
            }
            else
            {
                RevBayesCore::StochasticNode<double> *n2 = dynamic_cast<RevBayesCore::StochasticNode<double> *>( uReal );
                prop.removeUntransformedScalar(n2);

            }

        }
        else if ( sim != NULL )
        {

            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::AVMVNProposal &prop = static_cast<RevBayesCore::AVMVNProposal&>( m->getProposal() );

//            RevBayesCore::TypedDagNode< RevBayesCore::Simplex >* tmp = static_cast<const Simplex &>( args[0].getVariable()->getRevObject() ).getDagNode();
//            RevBayesCore::StochasticNode< RevBayesCore::Simplex > *n2 = static_cast<RevBayesCore::StochasticNode< RevBayesCore::Simplex > *>( tmp );

            RevBayesCore::StochasticNode<RevBayesCore::Simplex> *n = dynamic_cast<RevBayesCore::StochasticNode<RevBayesCore::Simplex > *>( sim->getDagNode() );
            prop.removeLogConstrainedSumVector(n);

        }
        else if ( uVector != NULL && uVector->getDagNode()->isStochastic() == true )
        {
            RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> > *the_node = dynamic_cast< RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> > * >( uVector->getDagNode() );
            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::AVMVNProposal &prop = static_cast<RevBayesCore::AVMVNProposal&>( m->getProposal() );

            if ( the_node != NULL )
            {
                prop.addUntransformedVector( the_node );
            }
            else
            {
                throw RbException("Could not remove the node because it isn't a stochastic nodes.");
            }
        }
        else
        {
            throw RbException() << "A problem occured when trying to remove " <<  args[0].getVariable()->getName() << " from the move.";
        }

        return NULL;
    }

    return Move::executeMethod( name, args, found );
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Move_AVMVN::getClassType(void)
{

    static std::string rev_type = "Move_AVMVN";

    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Move_AVMVN::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );

    return rev_type_spec;
}

/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_AVMVN::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "AVMVN";

    return c_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the scale move are:
 * (1) the variable which must be an integer.
 *
 * \return The member rules.
 */
const MemberRules& Move_AVMVN::getParameterRules(void) const
{

    static MemberRules memberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        memberRules.push_back( new ArgumentRule( "sigmaSquared"        , RealPos::getClassTypeSpec()  , "The scaling factor (strength) of the proposal.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RealPos(1.0) ) );
        memberRules.push_back( new ArgumentRule( "epsilon"             , RealPos::getClassTypeSpec()  , "The mixture weight of the post-learning move on a simple identity matrix.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new Probability(0.05) ) );
        memberRules.push_back( new ArgumentRule( "waitBeforeLearning"  , Natural::getClassTypeSpec()  , "The number of move attempts to wait before tracking the covariance of the variables.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(2500) ) );
        memberRules.push_back( new ArgumentRule( "waitBeforeUsing"     , Natural::getClassTypeSpec()  , "The number of move attempts to wait before using the learned covariance matrix.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(5000) ) );
        memberRules.push_back( new ArgumentRule( "maxUpdates"          , Natural::getClassTypeSpec()  , "The maximum number of updates to the empirical covariance matrix (matrix is only updated when MCMC tunes).", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(10000) ) );
        memberRules.push_back( new ArgumentRule( "tune"                , RlBoolean::getClassTypeSpec(), "Should we tune the scaling factor during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );

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
const TypeSpec& Move_AVMVN::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}



void Move_AVMVN::printValue(std::ostream &o) const
{

    o << "Move_AVMVN(?)";

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
void Move_AVMVN::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "sigmaSquared" )
    {
        sigmaSquared = var;
    }
    else if ( name == "epsilon" )
    {
        epsilon = var;
    }
    else if ( name == "waitBeforeLearning" )
    {
        waitBeforeLearning = var;
    }
    else if ( name == "waitBeforeUsing" )
    {
        waitBeforeUsing = var;
    }
    else if ( name == "maxUpdates" )
    {
        maxUpdates = var;
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

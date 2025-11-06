#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Move_UpDownTreeScale.h"
#include "UpDownScaleProposal.h"
#include "RbException.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "RlTimeTree.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "DagNode.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "Move.h"
#include "RbBoolean.h"
#include "Real.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "RlUtils.h"
#include "StochasticNode.h"

namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class RbVector; }


using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Move_UpDownTreeScale::Move_UpDownTreeScale() : Move()
{
    
    // add member methods
    
    // first, the argument rules
    ArgumentRules* addTreeArgRules                  = new ArgumentRules();
    ArgumentRules* addScalarArgRules                = new ArgumentRules();
    ArgumentRules* addModelVectorArgRules           = new ArgumentRules();
    ArgumentRules* addModelVectorPosArgRules        = new ArgumentRules();
    ArgumentRules* addCompositeVectorArgRules       = new ArgumentRules();
    ArgumentRules* addCompositeVectorPosArgRules    = new ArgumentRules();
    ArgumentRules* removeTreeArgRules               = new ArgumentRules();
    ArgumentRules* removeScalarArgRules             = new ArgumentRules();
    ArgumentRules* removeModelVectorArgRules        = new ArgumentRules();
    ArgumentRules* removeModelVectorPosArgRules     = new ArgumentRules();
    ArgumentRules* removeCompositeVectorArgRules    = new ArgumentRules();
    ArgumentRules* removeCompositeVectorPosArgRules = new ArgumentRules();
    
    
    // next, set the specific arguments
    addTreeArgRules->push_back(                     new ArgumentRule( "tree", TimeTree::getClassTypeSpec(),             "The tree variable to scale.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC) );
    addTreeArgRules->push_back(                     new ArgumentRule( "up"  , RlBoolean::getClassTypeSpec(),            "Scaling up or down?", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    addScalarArgRules->push_back(                   new ArgumentRule( "var" , Real::getClassTypeSpec(),                 "The variable to scale", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
    addScalarArgRules->push_back(                   new ArgumentRule( "up"  , RlBoolean::getClassTypeSpec(),            "Scaling up or down?", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    addModelVectorArgRules->push_back(              new ArgumentRule( "var" , ModelVector<Real>::getClassTypeSpec(),    "The variable to scale", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
    addModelVectorArgRules->push_back(              new ArgumentRule( "up"  , RlBoolean::getClassTypeSpec(),            "Scaling up or down?", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    addModelVectorPosArgRules->push_back(           new ArgumentRule( "var" , ModelVector<RealPos>::getClassTypeSpec(), "The variable to scale", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
    addModelVectorPosArgRules->push_back(           new ArgumentRule( "up"  , RlBoolean::getClassTypeSpec(),            "Scaling up or down?", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    addCompositeVectorArgRules->push_back(          new ArgumentRule( "var" , ModelVector<Real>::getClassTypeSpec(),    "The variable to scale", ArgumentRule::BY_REFERENCE, ArgumentRule::DETERMINISTIC ) );
    addCompositeVectorArgRules->push_back(          new ArgumentRule( "up"  , RlBoolean::getClassTypeSpec(),            "Scaling up or down?", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    addCompositeVectorPosArgRules->push_back(       new ArgumentRule( "var" , ModelVector<RealPos>::getClassTypeSpec(), "The variable to scale", ArgumentRule::BY_REFERENCE, ArgumentRule::DETERMINISTIC ) );
    addCompositeVectorPosArgRules->push_back(       new ArgumentRule( "up"  , RlBoolean::getClassTypeSpec(),            "Scaling up or down?", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    removeTreeArgRules->push_back(                  new ArgumentRule( "tree", TimeTree::getClassTypeSpec(),             "The tree variable to scale.",ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
    removeTreeArgRules->push_back(                  new ArgumentRule( "up"  , RlBoolean::getClassTypeSpec(),            "The variable to scale", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    removeScalarArgRules->push_back(                new ArgumentRule( "var" , Real::getClassTypeSpec(),                 "Scaling up or down?", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
    removeScalarArgRules->push_back(                new ArgumentRule( "up"  , RlBoolean::getClassTypeSpec(),            "The variable to scale", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    removeModelVectorArgRules->push_back(           new ArgumentRule( "var" , ModelVector<Real>::getClassTypeSpec(),    "Scaling up or down?", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
    removeModelVectorArgRules->push_back(           new ArgumentRule( "up"  , RlBoolean::getClassTypeSpec(),            "The variable to scale", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    removeModelVectorPosArgRules->push_back(        new ArgumentRule( "var" , ModelVector<RealPos>::getClassTypeSpec(), "Scaling up or down?", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
    removeModelVectorPosArgRules->push_back(        new ArgumentRule( "up"  , RlBoolean::getClassTypeSpec(),            "The variable to scale", ArgumentRule::BY_VALUE, ArgumentRule::ANY) );
    removeCompositeVectorArgRules->push_back(       new ArgumentRule( "var" , ModelVector<Real>::getClassTypeSpec(),    "Scaling up or down?", ArgumentRule::BY_REFERENCE, ArgumentRule::DETERMINISTIC ) );
    removeCompositeVectorArgRules->push_back(       new ArgumentRule( "up"  , RlBoolean::getClassTypeSpec(),            "The variable to scale", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    removeCompositeVectorPosArgRules->push_back(    new ArgumentRule( "var" , ModelVector<RealPos>::getClassTypeSpec(), "Scaling up or down?", ArgumentRule::BY_REFERENCE, ArgumentRule::DETERMINISTIC ) );
    removeCompositeVectorPosArgRules->push_back(    new ArgumentRule( "up"  , RlBoolean::getClassTypeSpec(),            "The variable to scale", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    
    
    // finally, create the methods
    methods.addFunction( new MemberProcedure( "addVariable", RlUtils::Void, addTreeArgRules) );
    methods.addFunction( new MemberProcedure( "addVariable", RlUtils::Void, addScalarArgRules) );
    methods.addFunction( new MemberProcedure( "addVariable", RlUtils::Void, addModelVectorArgRules) );
    methods.addFunction( new MemberProcedure( "addVariable", RlUtils::Void, addModelVectorPosArgRules) );
    methods.addFunction( new MemberProcedure( "addVariable", RlUtils::Void, addCompositeVectorArgRules) );
    methods.addFunction( new MemberProcedure( "addVariable", RlUtils::Void, addCompositeVectorPosArgRules) );
    methods.addFunction( new MemberProcedure( "removeVariable", RlUtils::Void, removeTreeArgRules) );
    methods.addFunction( new MemberProcedure( "removeVariable", RlUtils::Void, removeScalarArgRules) );
    methods.addFunction( new MemberProcedure( "removeVariable", RlUtils::Void, removeModelVectorArgRules) );
    methods.addFunction( new MemberProcedure( "removeVariable", RlUtils::Void, removeModelVectorPosArgRules) );
    methods.addFunction( new MemberProcedure( "removeVariable", RlUtils::Void, removeCompositeVectorArgRules) );
    methods.addFunction( new MemberProcedure( "removeVariable", RlUtils::Void, removeCompositeVectorPosArgRules) );
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
Move_UpDownTreeScale* Move_UpDownTreeScale::clone(void) const
{
    
    return new Move_UpDownTreeScale(*this);
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
void Move_UpDownTreeScale::constructInternalObject( void )
{
    
    // we free the memory first
    delete value;
    
    // now allocate a new up-down-scale move
    double l = static_cast<const RealPos &>( lambda->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    
    // finally create the internal move object
    RevBayesCore::UpDownScaleProposal *prop = new RevBayesCore::UpDownScaleProposal(l);
    
    value = new RevBayesCore::MetropolisHastingsMove(prop,w,t);
    
}


RevPtr<RevVariable> Move_UpDownTreeScale::executeMethod(const std::string& name, const std::vector<Argument>& args, bool &found)
{
    
    if ( name == "addVariable" )
    {
        found = true;
        
        TimeTree* uTree = dynamic_cast<TimeTree *>( &args[0].getVariable()->getRevObject() );
        Real* uReal = dynamic_cast<Real *>( &args[0].getVariable()->getRevObject() );
        ModelVector<Real>* uVector = dynamic_cast<ModelVector<Real> *>( &args[0].getVariable()->getRevObject() );
        ModelVector<RealPos>* upVector = dynamic_cast<ModelVector<RealPos> *>( &args[0].getVariable()->getRevObject() );
        bool up = static_cast<RlBoolean &>( args[1].getVariable()->getRevObject() ).getValue();
        
        if ( uTree != NULL )
        {
            RevBayesCore::StochasticNode<RevBayesCore::Tree> *the_node = dynamic_cast< RevBayesCore::StochasticNode<RevBayesCore::Tree> * >( uTree->getDagNode() );
            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::UpDownScaleProposal &prop = static_cast<RevBayesCore::UpDownScaleProposal&>( m->getProposal() );
            
            if ( the_node != NULL )
            {
                prop.addVariable( the_node, up );
            }
            else
            {
                throw RbException("Could not add the node because it isn't a stochastic nodes.");
            }
        }
        else if ( uReal != NULL )
        {
            RevBayesCore::StochasticNode<double> *the_node = dynamic_cast< RevBayesCore::StochasticNode<double> * >( uReal->getDagNode() );
            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::UpDownScaleProposal &prop = static_cast<RevBayesCore::UpDownScaleProposal&>( m->getProposal() );
            
            if ( the_node != NULL )
            {
                prop.addVariable( the_node, up );
            }
            else
            {
                throw RbException("Could not add the node because it isn't a stochastic nodes.");
            }
        }
        else if ( uVector != NULL && uVector->getDagNode()->isStochastic() == true )
        {
            
            RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> > *the_node = dynamic_cast< RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> > * >( uVector->getDagNode() );
            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::UpDownScaleProposal &prop = static_cast<RevBayesCore::UpDownScaleProposal&>( m->getProposal() );
            
            if ( the_node != NULL )
            {
                prop.addVariable( the_node, up );
            }
            else
            {
                throw RbException("Could not add the node because it isn't a stochastic nodes.");
            }
        }
        else if ( uVector != NULL && uVector->getDagNode()->isStochastic() == false )
        {
            
            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::UpDownScaleProposal &prop = static_cast<RevBayesCore::UpDownScaleProposal&>( m->getProposal() );
            
            RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* deterministic_vector = uVector->getDagNode();
            std::vector<const RevBayesCore::DagNode*> parents = deterministic_vector->getParents();
            for (std::vector<const RevBayesCore::DagNode*>::const_iterator it = parents.begin(); it != parents.end(); ++it)
            {
                const RevBayesCore::DagNode *tmp_node = *it;
                const RevBayesCore::StochasticNode<double> *the_node = dynamic_cast< const RevBayesCore::StochasticNode<double>* >( tmp_node );
                if ( the_node != NULL )
                {
                    prop.addVariable( const_cast<RevBayesCore::StochasticNode<double> *>( the_node ), up );
                }
                else
                {
                    throw RbException("Could not add the node because it isn't a vector of stochastic nodes.");
                }
            }
            
        }
        else if ( upVector != NULL && upVector->getDagNode()->isStochastic() == true )
        {
            
            RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> > *the_node = dynamic_cast< RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> > * >( upVector->getDagNode() );
            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::UpDownScaleProposal &prop = static_cast<RevBayesCore::UpDownScaleProposal&>( m->getProposal() );
            
            if ( the_node != NULL )
            {
                prop.addVariable( the_node, up );
            }
            else
            {
                throw RbException("Could not add the node because it isn't a stochastic nodes.");
            }
            
        }
        else if ( upVector != NULL && upVector->getDagNode()->isStochastic() == false )
        {
            
            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::UpDownScaleProposal &prop = static_cast<RevBayesCore::UpDownScaleProposal&>( m->getProposal() );
            
            RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* deterministic_vector = upVector->getDagNode();
            std::vector<const RevBayesCore::DagNode*> parents = deterministic_vector->getParents();
            for (std::vector<const RevBayesCore::DagNode*>::const_iterator it = parents.begin(); it != parents.end(); ++it)
            {
                const RevBayesCore::DagNode *tmp_node = *it;
                const RevBayesCore::StochasticNode<double> *the_node = dynamic_cast< const RevBayesCore::StochasticNode<double>* >( tmp_node );
                if ( the_node != NULL )
                {
                    prop.addVariable( const_cast<RevBayesCore::StochasticNode<double> *>( the_node ), up );
                }
                else
                {
                    throw RbException("Could not add the node because it isn't a vector of stochastic nodes.");
                }
            }
            
        }
        else
        {
            throw RbException() << "A problem occured when trying to add " << args[0].getVariable()->getName() << " to the move.";
        }
        
        return NULL;
    }
    else if ( name == "removeVariable" )
    {
        found = true;
        
        TimeTree* uTree = dynamic_cast<TimeTree *>( &args[0].getVariable()->getRevObject() );
        Real* uReal = dynamic_cast<Real *>( &args[0].getVariable()->getRevObject() );
        ModelVector<Real>* uVector = dynamic_cast<ModelVector<Real> *>( &args[0].getVariable()->getRevObject() );
        ModelVector<RealPos>* upVector = dynamic_cast<ModelVector<RealPos> *>( &args[0].getVariable()->getRevObject() );
        bool up = static_cast<RlBoolean &>( args[1].getVariable()->getRevObject() ).getValue();
        
        if ( uTree != NULL )
        {
            RevBayesCore::StochasticNode<RevBayesCore::Tree> *the_node = dynamic_cast< RevBayesCore::StochasticNode<RevBayesCore::Tree> * >( uTree->getDagNode() );
            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::UpDownScaleProposal &prop = static_cast<RevBayesCore::UpDownScaleProposal&>( m->getProposal() );
            
            if ( the_node != NULL )
            {
                prop.removeVariable( the_node, up );
            }
            else
            {
                throw RbException("Could not add the node because it isn't a stochastic nodes.");
            }
        }
        else if ( uReal != NULL )
        {
            RevBayesCore::StochasticNode<double> *the_node = dynamic_cast< RevBayesCore::StochasticNode<double> * >( uReal->getDagNode() );
            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::UpDownScaleProposal &prop = static_cast<RevBayesCore::UpDownScaleProposal&>( m->getProposal() );
            
            if ( the_node != NULL )
            {
                prop.removeVariable( the_node, up );
            }
            else
            {
                throw RbException("Could not add the node because it isn't a stochastic nodes.");
            }
        }
        else if ( uVector != NULL && uVector->getDagNode()->isStochastic() == true )
        {
            RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> > *the_node = dynamic_cast< RevBayesCore::StochasticNode<RevBayesCore::RbVector<double> > * >( uVector->getDagNode() );
            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::UpDownScaleProposal &prop = static_cast<RevBayesCore::UpDownScaleProposal&>( m->getProposal() );
            
            if ( the_node != NULL )
            {
                prop.removeVariable( the_node, up );
            }
            else
            {
                throw RbException("Could not add the node because it isn't a stochastic nodes.");
            }
        }
        else if ( uVector != NULL && uVector->getDagNode()->isStochastic() == false )
        {
            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::UpDownScaleProposal &prop = static_cast<RevBayesCore::UpDownScaleProposal&>( m->getProposal() );
            
            RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* deterministic_vector = uVector->getDagNode();
            std::vector<const RevBayesCore::DagNode*> parents = deterministic_vector->getParents();
            for (std::vector<const RevBayesCore::DagNode*>::const_iterator it = parents.begin(); it != parents.end(); ++it)
            {
                const RevBayesCore::StochasticNode<double> *the_node = dynamic_cast< const RevBayesCore::StochasticNode<double>* >( *it );
                if ( the_node != NULL )
                {
                    prop.removeVariable( const_cast<RevBayesCore::StochasticNode<double> *>( the_node ), up );
                }
                else
                {
                    throw RbException("Could not remove the node because it isn't a vector of stochastic nodes.");
                }
            }

        }
        else if ( upVector != NULL && upVector->getDagNode()->isStochastic() == false )
        {
            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::UpDownScaleProposal &prop = static_cast<RevBayesCore::UpDownScaleProposal&>( m->getProposal() );
            
            RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* deterministic_vector = upVector->getDagNode();
            std::vector<const RevBayesCore::DagNode*> parents = deterministic_vector->getParents();
            for (std::vector<const RevBayesCore::DagNode*>::const_iterator it = parents.begin(); it != parents.end(); ++it)
            {
                const RevBayesCore::StochasticNode<double> *the_node = dynamic_cast< const RevBayesCore::StochasticNode<double>* >( *it );
                if ( the_node != NULL )
                {
                    prop.removeVariable( const_cast<RevBayesCore::StochasticNode<double> *>( the_node ), up );
                }
                else
                {
                    throw RbException("Could not remove the node because it isn't a vector of stochastic nodes.");
                }
            }
            
        }
        else
        {
            throw RbException() << "A problem occured when trying to add " << args[0].getVariable()->getName() << " to the move.";
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
const std::string& Move_UpDownTreeScale::getClassType(void)
{
    
    static std::string rev_type = "Move_UpDownTreeScale";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Move_UpDownTreeScale::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_UpDownTreeScale::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "UpDownScale";
    
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
const MemberRules& Move_UpDownTreeScale::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {

        move_member_rules.push_back( new ArgumentRule( "lambda"      , RealPos::getClassTypeSpec()  , "The scaling factor (strength) of the proposal.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RealPos(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune"        , RlBoolean::getClassTypeSpec(), "Should we tune the scaling factor during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Move_UpDownTreeScale::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



void Move_UpDownTreeScale::printValue(std::ostream &o) const
{
    
    o << "Move_UpDownTreeScale(?)";
    
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
void Move_UpDownTreeScale::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "lambda" )
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

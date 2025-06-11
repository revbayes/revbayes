#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Move_RateAgeSubtreeProposal.h"
#include "RateAgeSubtreeProposal.h"
#include "Probability.h"
#include "RbException.h"
#include "RbMathLogic.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "RlSimplex.h"
#include "RlString.h"
#include "RlTree.h"
#include "RlTimeTree.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "WorkspaceVector.h"
#include "VectorFunction.h"


using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Move_RateAgeSubtreeProposal::Move_RateAgeSubtreeProposal() : Move()
{
    
    // first, the argument rules
    ArgumentRules* addModelVectorArgRules           = new ArgumentRules();
//    ArgumentRules* removeModelVectorArgRules        = new ArgumentRules();
    
    addModelVectorArgRules->push_back(              new ArgumentRule( "var"        , ModelVector<RealPos>::getClassTypeSpec(),    "The variable to move"             , ArgumentRule::BY_REFERENCE, ArgumentRule::DETERMINISTIC ) );
//    removeModelVectorArgRules->push_back(           new ArgumentRule( "var"        , ModelVector<Real>::getClassTypeSpec(),    "The variable to move"             , ArgumentRule::BY_REFERENCE, ArgumentRule::DETERMINISTIC ) );
    
    // finally, create the methods
    methods.addFunction( new MemberProcedure( "addRates",    RlUtils::Void, addModelVectorArgRules) );
//    methods.addFunction( new MemberProcedure( "removeRates", RlUtils::Void, removeModelVectorArgRules) );
    
}

RevPtr<RevVariable> Move_RateAgeSubtreeProposal::executeMethod(const std::string& name, const std::vector<Argument>& args, bool &found)
{
    
    if ( name == "addRates" )
    {
        found = true;
        
        RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
        RevBayesCore::RateAgeSubtreeProposal &prop = static_cast<RevBayesCore::RateAgeSubtreeProposal&>( m->getProposal() );

        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* tmpRates = static_cast<const ModelVector<RealPos> &>( args[0].getVariable()->getRevObject() ).getDagNode();
        std::vector< RevBayesCore::StochasticNode<double> *> rates;
        RevBayesCore::DeterministicNode< RevBayesCore::RbVector<double> >*dnode = static_cast< RevBayesCore::DeterministicNode< RevBayesCore::RbVector<double> > *>( tmpRates );
        
        RevBayesCore::VectorFunction<double>* funcVec = dynamic_cast<RevBayesCore::VectorFunction<double>*>( &dnode->getFunction() );
        if ( funcVec == NULL )
        {
            throw RbException("Wrong argument type for the rates vector. We expect a vector of iid elements.");
        }
        const std::vector<const RevBayesCore::TypedDagNode<double>* >& pars = funcVec->getVectorParameters();
        
        for (std::vector<const RevBayesCore::TypedDagNode<double>* >::const_iterator it = pars.begin(); it != pars.end(); ++it)
        {
            rates.push_back( const_cast<RevBayesCore::StochasticNode<double>* >(static_cast<const RevBayesCore::StochasticNode<double>* >( *it ) ) );
        }

        prop.addRates( rates );
        
//        prop.addRates( rates );
        
//        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* tmpRates = static_cast<const ModelVector<RealPos> &>( rates->getRevObject() ).getDagNode();
        
//        ModelVector<RealPos>* uVector = dynamic_cast<ModelVector<RealPos> *>( &args[0].getVariable()->getRevObject() );
//
//        if ( uVector != NULL && uVector->getDagNode()->isStochastic() == true )
//        {
//
//            RevBayesCore::DeterministicNode<RevBayesCore::RbVector<double> > *the_node = dynamic_cast< RevBayesCore::DeterministicNode<RevBayesCore::RbVector<double> > * >( uVector->getDagNode() );
//            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
//            RevBayesCore::RateAgeProposal &prop = static_cast<RevBayesCore::RateAgeProposal&>( m->getProposal() );
//
//            RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* tmpRates = static_cast<const ModelVector<RealPos> &>( rates->getRevObject() ).getDagNode();
//            std::vector< RevBayesCore::StochasticNode<double> *> rates;
//            RevBayesCore::DeterministicNode< RevBayesCore::RbVector<double> >*dnode = static_cast< RevBayesCore::DeterministicNode< RevBayesCore::RbVector<double> > *>( tmpRates );
//
//            if ( the_node != NULL )
//            {
//                prop.addRates( the_node );
//            }
//            else
//            {
//                throw RbException("Could not add the node because it isn't a stochastic nodes.");
//            }
//        }
//        else
//        {
//            throw RbException() << "A problem occured when trying to add " << args[0].getVariable()->getName() << " to the move.";
//        }
        
        return NULL;
    }
    else if ( name == "removeRates" )
    {
        
    }
    
    return Move::executeMethod( name, args, found );
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
Move_RateAgeSubtreeProposal* Move_RateAgeSubtreeProposal::clone(void) const
{
    
    return new Move_RateAgeSubtreeProposal(*this);
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
void Move_RateAgeSubtreeProposal::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new element scaling move
    double a = static_cast<const RealPos &>( alpha->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    double r = static_cast<const RealPos &>( tune_target->getRevObject() ).getValue();
    size_t del = static_cast<const Natural &>( delay->getRevObject() ).getValue();
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();

    // get the tree
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tmp = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );

    // finally create the internal move object
    RevBayesCore::RateAgeSubtreeProposal *prop = new RevBayesCore::RateAgeSubtreeProposal(n, a, r);
    
    value = new RevBayesCore::MetropolisHastingsMove(prop,w,del,t);
    
}

/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Move_RateAgeSubtreeProposal::getClassType(void)
{
    
    static std::string rev_type = "Move_RateAgeSubtreeProposal";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Move_RateAgeSubtreeProposal::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}

/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_RateAgeSubtreeProposal::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "RateAgeSubtreeProposal";
    
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
const MemberRules& Move_RateAgeSubtreeProposal::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new ArgumentRule( "tree",   Tree::getClassTypeSpec(),      "The tree to which the branch-rates belong.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "alpha",  RealPos::getClassTypeSpec(),   "The scaling factor (strength) of this move.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new Real(1.0) ) );
        memberRules.push_back( new ArgumentRule( "tune"  , RlBoolean::getClassTypeSpec(), "Should we tune the scaling factor during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );
        
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
const TypeSpec& Move_RateAgeSubtreeProposal::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



void Move_RateAgeSubtreeProposal::printValue(std::ostream &o) const
{
    
    o << "Move_RateAgeSubtreeProposal(?)";

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
void Move_RateAgeSubtreeProposal::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "alpha" )
    {
        alpha = var;
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

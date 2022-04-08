#include <stddef.h>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Move_SubtreeScaleRelativeBranchLengthSimplexSingleTreeLength.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlBranchLengthTree.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"


namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Move_SubtreeScaleRelativeBranchLengthSimplexSingleTreeLength::Move_SubtreeScaleRelativeBranchLengthSimplexSingleTreeLength() : Move()
{
    
    // first, the argument rules
    ArgumentRules* addTreeLengthScalarRules = new ArgumentRules();
    
    // next, set the specific arguments
    addTreeLengthScalarRules->push_back( new ArgumentRule( "treeLength" , RealPos::getClassTypeSpec(), "The tree length scalar.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
    
    // finally, create the methods
    methods.addFunction( new MemberProcedure( "addTreeLengthScalar", RlUtils::Void, addTreeLengthScalarRules) );
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the move.
 */
Move_SubtreeScaleRelativeBranchLengthSimplexSingleTreeLength* Move_SubtreeScaleRelativeBranchLengthSimplexSingleTreeLength::clone(void) const
{
    
    return new Move_SubtreeScaleRelativeBranchLengthSimplexSingleTreeLength(*this);
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
void Move_SubtreeScaleRelativeBranchLengthSimplexSingleTreeLength::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    double l = static_cast<const RealPos &>( lambda->getRevObject() ).getValue();
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    double r = static_cast<const Probability &>( tuneTarget->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp_tree = static_cast<const BranchLengthTree &>( tree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *tr = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp_tree );
    
    RevBayesCore::TypedDagNode< RevBayesCore::Simplex>* tmp_relbls = static_cast<const Simplex &>( relative_branch_lengths->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode< RevBayesCore::Simplex> *relbls = static_cast<RevBayesCore::StochasticNode< RevBayesCore::Simplex > *>( tmp_relbls );
    
    RevBayesCore::Proposal *p = new RevBayesCore::SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal(tr, relbls, l, r);
    value = new RevBayesCore::MetropolisHastingsMove(p, w, t);
    
}


RevPtr<RevVariable> Move_SubtreeScaleRelativeBranchLengthSimplexSingleTreeLength::executeMethod(const std::string& name, const std::vector<Argument>& args, bool &found)
{
    
    if ( name == "addTreeLengthScalar" )
    {
        found = true;
        
        RealPos* tmp_tl = dynamic_cast<RealPos *>( &args[0].getVariable()->getRevObject() );
        
        if ( tmp_tl != NULL )
        {
            RevBayesCore::StochasticNode<double> *tl = dynamic_cast< RevBayesCore::StochasticNode<double> * >( tmp_tl->getDagNode() );
            
            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal &prop = static_cast<RevBayesCore::SubtreeScaleRelativeBranchLengthSimplexSingleTreeLengthProposal&>( m->getProposal() );
            
            if ( tl != NULL )
            {
                prop.addTreeLengthScalar( tl );
            }
            else
            {
                throw RbException("Could not add the node because it isn't a stochastic nodes.");
            }
        }
        
        else
        {
            throw RbException("A problem occured when trying to add " + args[0].getVariable()->getName() + " to the move.");
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
const std::string& Move_SubtreeScaleRelativeBranchLengthSimplexSingleTreeLength::getClassType(void)
{
    
    static std::string rev_type = "Move_SubtreeScaleRelativeBranchLengthSimplexSingleTreeLength";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Move_SubtreeScaleRelativeBranchLengthSimplexSingleTreeLength::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_SubtreeScaleRelativeBranchLengthSimplexSingleTreeLength::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "SubtreeScaleRelativeBranchLengthSimplexSingleTreeLength";
    
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
const MemberRules& Move_SubtreeScaleRelativeBranchLengthSimplexSingleTreeLength::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        memberRules.push_back( new ArgumentRule( "tree", BranchLengthTree::getClassTypeSpec(), "The tree variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        memberRules.push_back( new ArgumentRule( "relativeBranchLengths", Simplex::getClassTypeSpec(), "The simplex of relative branch lengths.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC)  );
        
        memberRules.push_back( new ArgumentRule( "lambda", RealPos::getClassTypeSpec(), "The strength of the proposal.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Real(1.0) ) );
        memberRules.push_back( new ArgumentRule( "tune", RlBoolean::getClassTypeSpec(), "Should we tune this move during burnin?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( true ) ) );
        
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
const TypeSpec& Move_SubtreeScaleRelativeBranchLengthSimplexSingleTreeLength::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/**
 * Print the value for the user.
 */
void Move_SubtreeScaleRelativeBranchLengthSimplexSingleTreeLength::printValue(std::ostream &o) const
{
    
    o << "SubtreeScaleRelativeBranchLengthSimplexSingleTreeLength(";
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
void Move_SubtreeScaleRelativeBranchLengthSimplexSingleTreeLength::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" ) {
        tree = var;
    }
    else if ( name == "relativeBranchLengths" ) {
        relative_branch_lengths = var;
    }
    else if ( name == "lambda" ) {
        lambda = var;
    }
    else if ( name == "tune" ) {
        tune = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
    
}





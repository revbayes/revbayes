#include "Func_chronoToPhylo.h"
#include "ChronoToPhyloFunction.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlBranchLengthTree.h"
#include "RlDeterministicNode.h"
#include "RlTimeTree.h"
#include "RlTree.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_chronoToPhylo::Func_chronoToPhylo( void ) : TypedFunction<BranchLengthTree>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_chronoToPhylo* Func_chronoToPhylo::clone( void ) const
{
    
    return new Func_chronoToPhylo( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::Tree>* Func_chronoToPhylo::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* chrono = static_cast<const Tree&>( args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* rates = static_cast<const ModelVector<RealPos> &>( args[1].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::ChronoToPhyloFunction* f = new RevBayesCore::ChronoToPhyloFunction( chrono, rates );
    
    return f;
}


///** Execute function */
//RevPtr<RevVariable> Func_chronoToPhylo::execute( void )
//{
//    RevBayesCore::TypedFunction<RevBayesCore::Tree>* d = createFunction();
//    RevBayesCore::DeterministicNode<RevBayesCore::Tree>* rv;
//    rv = new DeterministicNode<RevBayesCore::Tree>("", d, this->clone());
//    
//    RevObject& ro = args[0].getVariable()->getRevObject();
//    if ( ro.isType( TimeTree::getClassTypeSpec() ) )
//    {
//        return new RevVariable( new TimeTree(rv) );
//    }
//    else if ( ro.isType( BranchLengthTree::getClassTypeSpec() ) )
//    {
//        return new RevVariable( new BranchLengthTree(rv) );
//    }
//    
//    return new RevVariable( new Tree(rv) );
//}


/* Get argument rules */
const ArgumentRules& Func_chronoToPhylo::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "tree",  TimeTree::getClassTypeSpec(), "The tree variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rates" , ModelVector<RealPos>::getClassTypeSpec() , "The branch rates we use to multiply/rescale the branches with.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_chronoToPhylo::getClassType(void)
{
    
    static std::string rev_type = "Func_chronoToPhylo";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_chronoToPhylo::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}



/**
 * Get the primary Rev name for this function.
 */
std::string Func_chronoToPhylo::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnChronoToPhylo";
    
    return f_name;
}


const TypeSpec& Func_chronoToPhylo::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

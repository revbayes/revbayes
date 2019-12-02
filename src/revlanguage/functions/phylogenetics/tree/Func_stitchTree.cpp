#include "Func_stitchTree.h"
#include "ModelVector.h"
#include "StitchTreeFunction.h"
#include "RlBoolean.h"
#include "Real.h"
#include "RealPos.h"
#include "RlTree.h"
#include "RlTimeTree.h"
#include "RlDeterministicNode.h"
#include "RlTaxon.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_stitchTree::Func_stitchTree( void ) : TypedFunction<Tree>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_stitchTree* Func_stitchTree::clone( void ) const
{
    
    return new Func_stitchTree( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::Tree>* Func_stitchTree::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* bt = static_cast<const Tree&>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::Tree> >* pc = static_cast<const ModelVector<TimeTree> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::RbVector<RevBayesCore::RbVector<RevBayesCore::Taxon> > pt = static_cast<const ModelVector<ModelVector<Taxon> > &>( this->args[2].getVariable()->getRevObject() ).getValue();
  
    RevBayesCore::StitchTreeFunction* f = new RevBayesCore::StitchTreeFunction( bt, pc, pt );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_stitchTree::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "backboneTree",              Tree::getClassTypeSpec(), "The backbone tree variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY) );
        argumentRules.push_back( new ArgumentRule( "patchClades" , ModelVector<TimeTree>::getClassTypeSpec(), "The patch clade variables.",  ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY) );
        argumentRules.push_back( new ArgumentRule( "patchTaxa"  , ModelVector<ModelVector<Taxon> >::getClassTypeSpec(), "The taxa per patch clade to be stitched on to the backbone tree. Shares order with patchClades.",  ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY) );
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_stitchTree::getClassType(void)
{
    
    static std::string rev_type = "Func_stitchTree";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_stitchTree::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}



/**
 * Get the primary Rev name for this function.
 */
std::string Func_stitchTree::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnStitchTree";
    
    return f_name;
}


const TypeSpec& Func_stitchTree::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

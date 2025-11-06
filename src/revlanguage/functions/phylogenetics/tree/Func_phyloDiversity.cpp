#include "Func_phyloDiversity.h"

#include <cstddef>

#include "ModelVector.h"
#include "PhyloDiversityFunction.h"
#include "RealPos.h"
#include "RlBoolean.h"
#include "RlClade.h"
#include "RlTypedFunction.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlTree.h"
#include "TypeSpec.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class Clade; }
namespace RevBayesCore { class Tree; }


using namespace RevLanguage;


Func_phyloDiversity::Func_phyloDiversity( void ) : TypedFunction<RealPos>( )
{
    
}


Func_phyloDiversity* Func_phyloDiversity::clone( void ) const
{
    
    return new Func_phyloDiversity( *this );
}


RevBayesCore::TypedFunction<double>* Func_phyloDiversity::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const Tree &>(             this->args[0].getVariable()->getRevObject() ).getDagNode();
    const RevBayesCore::Clade& sample                             = static_cast<const Clade &>(                this->args[1].getVariable()->getRevObject() ).getValue();
    bool includeRoot                                              = static_cast<const RlBoolean &>(            this->args[2].getVariable()->getRevObject() ).getValue();
    
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* w = NULL;
    if ( this->args[3].getVariable() != NULL && this->args[3].getVariable()->getRevObject() != RevNullObject::getInstance() )
    {
        w = static_cast<const ModelVector<RealPos> &>( this->args[3].getVariable()->getRevObject() ).getDagNode();
    }
    
    RevBayesCore::PhyloDiversityFunction* f = new RevBayesCore::PhyloDiversityFunction( tau, sample, includeRoot, w );
    return f;
}


const ArgumentRules& Func_phyloDiversity::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "tree",        Tree::getClassTypeSpec(),                 "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "sample",      Clade::getClassTypeSpec(),                "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "includeRoot", RlBoolean::getClassTypeSpec(),            "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        argumentRules.push_back( new ArgumentRule( "weights",     ModelVector<RealPos>::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new ModelVector<RealPos>() ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_phyloDiversity::getClassType(void)
{
    
    static std::string rev_type = "Func_phyloDiversity";
    
    return rev_type;
}


const TypeSpec& Func_phyloDiversity::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedFunction<RealPos>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_phyloDiversity::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnPD";
    
    return f_name;
}


const TypeSpec& Func_phyloDiversity::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


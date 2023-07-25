#ifndef Func_collapseSA_h
#define Func_collapseSA_h

#include <string>
#include <iosfwd>
#include <vector>

#include "RlTree.h"
#include "RlTypedFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RevPtr.h"
#include "RlDeterministicNode.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "GenericFunction.h"

RevBayesCore::Tree* collapseSampledAncestorsFunc(const RevBayesCore::Tree& tree)
{
    auto tree2 = tree.clone();

    tree2->collapseSampledAncestors();

    return tree2;
}

namespace RevLanguage {

template <typename T>
class Func_collapseSA :  public TypedFunction<T> {

public:

    // Basic utility functions
    Func_collapseSA<T>*                                                 clone(void) const                                          //!< Clone the object
    {
        return new Func_collapseSA<T>( *this );
    }

    static const std::string&                                           getClassType(void)                                         //!< Get Rev type
    {
        static std::string rev_type = "Func_collapseSA";
    
        return rev_type;
    }

    static const TypeSpec&                                              getClassTypeSpec(void)                                     //!< Get class type spec
    {
        static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
        return rev_type_spec;
    }

    std::string                                                         getFunctionName(void) const                                //!< Get the primary name of the function in Rev
    {
        // create a name variable that is the same for all instance of this class
        std::string f_name = "fnCollapseSA";
    
        return f_name;

    }
    const TypeSpec&                                                     getTypeSpec(void) const                                    //!< Get the type spec of the instance
    {
        static TypeSpec type_spec = getClassTypeSpec();
    
        return type_spec;
    }

    // Function functions you have to override
    RevBayesCore::TypedFunction< RevBayesCore::Tree> *                  createFunction(void) const                                 //!< Create internal function object
    {
        RevBayesCore::TypedDagNode< RevBayesCore::Tree >* tree = static_cast<const Tree &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

        return RevBayesCore::generic_function_ptr< RevBayesCore::Tree >( collapseSampledAncestorsFunc, tree);

    }

    const ArgumentRules&                                                getArgumentRules(void) const                               //!< Get argument rules
    {
        static ArgumentRules argumentRules = ArgumentRules();
        static bool          rules_set = false;
    
        if ( !rules_set )
        {

            argumentRules.push_back( new ArgumentRule( "tree", T::getClassTypeSpec(), "The tree variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
            rules_set = true;
        }

        return argumentRules;
    }
};

};


#endif /* Func_collapseSA_h */

#include "ModelVector.h"
#include "Natural.h"
#include "RbUtil.h"
#include "RlBranchLengthTree.h"
#include "RlMemberFunction.h"
#include "RlString.h"
#include "RealPos.h"
#include "TypeSpec.h"

#include <sstream>

using namespace RevLanguage;

/** Default constructor */
BranchLengthTree::BranchLengthTree(void) : Tree()
{
    
    ArgumentRules* rerootArgRules = new ArgumentRules();
    rerootArgRules->push_back( new ArgumentRule("leaf", RlString::getClassTypeSpec(), "The outgroup leaf.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    
    methods.addFunction( new MemberProcedure( "reroot", RlUtils::Void,  rerootArgRules       ) );
    
}

/** Construct from bool */
BranchLengthTree::BranchLengthTree(RevBayesCore::Tree *t) : Tree( t )
{
    
    ArgumentRules* rerootArgRules = new ArgumentRules();
    rerootArgRules->push_back( new ArgumentRule("leaf", RlString::getClassTypeSpec(), "The outgroup leaf.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    
    methods.addFunction( new MemberProcedure( "reroot", RlUtils::Void,  rerootArgRules       ) );
    
}

/** Construct from bool */
BranchLengthTree::BranchLengthTree(const RevBayesCore::Tree &t) : Tree( new RevBayesCore::Tree( t ) )
{
    
    ArgumentRules* rerootArgRules = new ArgumentRules();
    rerootArgRules->push_back( new ArgumentRule("leaf", RlString::getClassTypeSpec(), "The outgroup leaf.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    
    methods.addFunction( new MemberProcedure( "reroot", RlUtils::Void,  rerootArgRules       ) );
    
}

/** Construct from bool */
BranchLengthTree::BranchLengthTree(RevBayesCore::TypedDagNode<RevBayesCore::Tree> *n) : Tree( n )
{
    
    ArgumentRules* rerootArgRules = new ArgumentRules();
    rerootArgRules->push_back( new ArgumentRule("leaf", RlString::getClassTypeSpec(), "The outgroup leaf.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    
    methods.addFunction( new MemberProcedure( "reroot", RlUtils::Void,  rerootArgRules       ) );
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
BranchLengthTree* BranchLengthTree::clone(void) const
{
    
	return new BranchLengthTree(*this);
}


/* Map calls to member methods */
RevLanguage::RevPtr<RevLanguage::RevVariable> BranchLengthTree::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "reroot")
    {
        
        const RevObject& st = args[0].getVariable()->getRevObject();
        if ( st.isType( RlString::getClassTypeSpec() ) )
        {
            std::string n = std::string( static_cast<const RlString&>( st ).getValue() );
            this->dag_node->getValue().reroot(n, true);
        }
        return NULL;
        
    }
    
    return Tree::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& BranchLengthTree::getClassType(void)
{
    
    static std::string rev_type = "BranchLengthTree";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& BranchLengthTree::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Tree::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/** Get type spec */
const TypeSpec& BranchLengthTree::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


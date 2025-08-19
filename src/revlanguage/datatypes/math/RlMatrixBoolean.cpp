#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

#include "ConstantNode.h"
#include "MatrixBoolean.h"
#include "Natural.h"
#include "ModelVector.h"
#include "RlBoolean.h"
#include "RlMatrixBoolean.h"
#include "RlMemberFunction.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MemberFunction.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlDeterministicNode.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

/* Default constructor */
MatrixBoolean::MatrixBoolean(void) : ModelObject<RevBayesCore::MatrixBoolean>( new RevBayesCore::MatrixBoolean( 1, 1, false ) )
{
    
    // initialize the member methods
    initializeMethods();
}

MatrixBoolean::MatrixBoolean(const RevBayesCore::MatrixBoolean &v) : ModelObject<RevBayesCore::MatrixBoolean>( v.clone() )
{
    
    // initialize the member methods
    initializeMethods();
}

MatrixBoolean::MatrixBoolean(RevBayesCore::MatrixBoolean *m) : ModelObject<RevBayesCore::MatrixBoolean>( m )
{
    
    // initialize the member methods
    initializeMethods();
}

MatrixBoolean::MatrixBoolean( RevBayesCore::TypedDagNode<RevBayesCore::MatrixBoolean> * mat ) : ModelObject<RevBayesCore::MatrixBoolean>( mat )
{
    
    // initialize the member methods
    initializeMethods();
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
MatrixBoolean* MatrixBoolean::clone(void) const
{
    
    return new MatrixBoolean(*this);
}


/** Get Rev type of object */
const std::string& MatrixBoolean::getClassType(void)
{
    
    static std::string rev_type = "MatrixBoolean";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& MatrixBoolean::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/** Get type spec */
const TypeSpec& MatrixBoolean::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/* Map calls to member methods */
RevPtr<RevVariable> MatrixBoolean::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "column")
    {
        found = true;
        
        const Natural& index = static_cast<const Natural&>( args[0].getVariable()->getRevObject() );
        size_t i = index.getValue() - 1;
        
        RevBayesCore::RbVector<RevBayesCore::Boolean> m = this->dag_node->getValue().getColumn( i );
        ModelVector< RlBoolean > *bl_vect = new ModelVector< RlBoolean >;
        for(size_t j = 0; j < m.size(); ++j)
        {
            bl_vect->push_back( RlBoolean(m[j]) );
        }
        
        return new RevVariable( bl_vect );
    }
    else if (name == "ncolumns")
    {
        found = true;
        
        size_t s = this->dag_node->getValue().getNumberOfColumns();
        return new RevVariable( new Natural( s ) );
    }
    else if (name == "nrows")
    {
        found = true;
        
        size_t s = this->dag_node->getValue().getNumberOfRows();
        return new RevVariable( new Natural( s ) );
    }
    else if (name == "size")
    {
        found = true;
        
        size_t s = this->dag_node->getValue().size();
        return new RevVariable( new Natural( s ) );
    }
    
    return ModelObject<RevBayesCore::MatrixBoolean>::executeMethod( name, args, found );
}


void MatrixBoolean::initializeMethods( void )
{
    // Add method for call "x[]" as a function
    ArgumentRules* squareBracketArgRules = new ArgumentRules();
    squareBracketArgRules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the row.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberFunction<MatrixBoolean, ModelVector< RlBoolean > >("[]", this, squareBracketArgRules ) );
    
    // add method for call "column" as a function
    ArgumentRules* columnArgRules = new ArgumentRules();
    columnArgRules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the column.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "column", ModelVector< RlBoolean >::getClassTypeSpec(), columnArgRules ) );
    
    // add method for call "size" as a function
    ArgumentRules* sizeArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "size", Natural::getClassTypeSpec(), sizeArgRules) );
    
    // add method for call "ncolumns" as a function
    ArgumentRules* nColumnsArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "ncolumns", Natural::getClassTypeSpec(), nColumnsArgRules) );

    // add method for call "nrows" as a function
    ArgumentRules* nRowsArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nrows", Natural::getClassTypeSpec(), nRowsArgRules) );

}

#include <sstream>
#include <string>
#include <vector>

#include "ConstantNode.h"
#include "MatrixReal.h"
#include "Natural.h"
#include "ModelVector.h"
#include "RlMatrixRealPos.h"
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
#include "RbVector.h"
#include "RealPos.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlDeterministicNode.h"
#include "RlMatrixReal.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

/* Default constructor */
MatrixRealPos::MatrixRealPos(void) : MatrixReal()
{
    
    // initialize the member methods
    initializeMethods();
}

MatrixRealPos::MatrixRealPos(const RevBayesCore::MatrixReal& from) : MatrixReal( from )
{
    
    // initialize the member methods
    initializeMethods();
}

MatrixRealPos::MatrixRealPos(RevBayesCore::MatrixReal* m) : MatrixReal( m )
{
    
    // initialize the member methods
    initializeMethods();
}

MatrixRealPos::MatrixRealPos( RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal> * mat ) : MatrixReal( mat )
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
MatrixRealPos* MatrixRealPos::clone(void) const
{
    
    return new MatrixRealPos(*this);
}


/* Map calls to member methods */
RevPtr<RevVariable> MatrixRealPos::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "min")
    {
        
        found = true;
        
        double m = this->dag_node->getValue().getMin();
        return new RevVariable( new RealPos( m ) );
    }
    else if (name == "max")
    {
        found = true;
        
        double m = this->dag_node->getValue().getMax();
        return new RevVariable( new RealPos( m ) );
    }
    else if (name == "column")
    {
        found = true;
        
        const Natural& index = static_cast<const Natural&>( args[0].getVariable()->getRevObject() );
        int i = (int)index.getValue() - 1;
        
        RevBayesCore::RbVector<double> m = this->dag_node->getValue().getColumn( i );

        return new RevVariable( new ModelVector<RealPos>( m ) );
    }
    
    return MatrixReal::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& MatrixRealPos::getClassType(void)
{
    
    static std::string rev_type = "MatrixRealPos";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& MatrixRealPos::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( MatrixReal::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


///** Return member rules (no members) */
//const MemberRules& MatrixRealPos::getParameterRules(void) const
//{
//
//    static MemberRules memberRules;
//    static bool rules_set = false;
//
//    if ( !rules_set )
//    {
//
//        memberRules.push_back( new ArgumentRule("x", RevObject::getClassTypeSpec(), ArgumentRule::BY_CONSTANT_REFERENCE ) );
//        memberRules.push_back( new Ellipsis( RevObject::getClassTypeSpec() ) );
//
//        rules_set = true;
//    }
//
//    return memberRules;
//}


/** Get type spec */
const TypeSpec& MatrixRealPos::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


void MatrixRealPos::initializeMethods( void )
{
    
    methods.eraseFunction( "[]" );
    methods.eraseFunction( "min" );
    methods.eraseFunction( "max" );
    methods.eraseFunction( "column" );
    
    // Add method for call "x[]" as a function
    ArgumentRules* squareBracketArgRules = new ArgumentRules();
    squareBracketArgRules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the row.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberFunction<MatrixRealPos,ModelVector<RealPos> >("[]", this, squareBracketArgRules ) );
    
    // add method for call "min" as a function
    ArgumentRules* minArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "min", RealPos::getClassTypeSpec(), minArgRules) );
    
    // add method for call "max" as a function
    ArgumentRules* maxArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "max", RealPos::getClassTypeSpec(), maxArgRules) );
    
    // add method for call "column" as a function
    ArgumentRules* columnArgRules = new ArgumentRules();
    columnArgRules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the column.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "column", ModelVector<RealPos>::getClassTypeSpec(), columnArgRules ) );
    
}



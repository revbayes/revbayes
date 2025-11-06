#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

#include "ConstantNode.h"
#include "MatrixReal.h"
#include "Natural.h"
#include "ModelVector.h"
#include "Real.h"
#include "RlMatrixReal.h"
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
MatrixReal::MatrixReal(void) : ModelObject<RevBayesCore::MatrixReal>( new RevBayesCore::MatrixReal(1,1,0) )
{
    
    // initialize the member methods
    initializeMethods();
}

MatrixReal::MatrixReal(const RevBayesCore::MatrixReal& from) : ModelObject<RevBayesCore::MatrixReal>( from.clone() )
{
    
    // initialize the member methods
    initializeMethods();
}

MatrixReal::MatrixReal(RevBayesCore::MatrixReal* m) : ModelObject<RevBayesCore::MatrixReal>( m )
{
    
    // initialize the member methods
    initializeMethods();
}

MatrixReal::MatrixReal( RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal> * mat ) : ModelObject<RevBayesCore::MatrixReal>( mat )
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
MatrixReal* MatrixReal::clone(void) const
{
    
	return new MatrixReal(*this);
}


/* Map calls to member methods */
RevPtr<RevVariable> MatrixReal::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "min")
    {
        
        found = true;
        
        double m = this->dag_node->getValue().getMin();
        return new RevVariable( new Real( m ) );
    }
    else if (name == "max")
    {
        found = true;
        
        double m = this->dag_node->getValue().getMax();
        return new RevVariable( new Real( m ) );
    }
    else if (name == "column")
    {
        found = true;
        
        const Natural& index = static_cast<const Natural&>( args[0].getVariable()->getRevObject() );
        int i = index.getValue() - 1;
        
        RevBayesCore::RbVector<double> m = this->dag_node->getValue().getColumn( i );
        return new RevVariable( new ModelVector<Real>( m ) );
    }
    else if (name == "nrows")
    {
        found = true;
        
        size_t s = this->dag_node->getValue().getNumberOfRows();
        return new RevVariable( new Natural( s ) );
    }
    else if (name == "ncolumns")
    {
        found = true;
        
        size_t s = this->dag_node->getValue().getNumberOfColumns();
        return new RevVariable( new Natural( s ) );
    }
    else if (name == "size")
    {
        found = true;
        
        size_t s = this->dag_node->getValue().size();
        return new RevVariable( new Natural( s ) );
    }
    
    return ModelObject<RevBayesCore::MatrixReal>::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& MatrixReal::getClassType(void)
{
    
    static std::string rev_type = "MatrixReal";
    
	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& MatrixReal::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


///** Return member rules (no members) */
//const MemberRules& MatrixReal::getParameterRules(void) const
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
const TypeSpec& MatrixReal::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


void MatrixReal::initializeMethods( void )
{
    // Add method for call "x[]" as a function
    ArgumentRules* squareBracketArgRules = new ArgumentRules();
    squareBracketArgRules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the row.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberFunction<MatrixReal,ModelVector<Real> >("[]", this, squareBracketArgRules ) );
    
    // add method for call "min" as a function
    ArgumentRules* minArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "min", Real::getClassTypeSpec(), minArgRules) );
    
    // add method for call "max" as a function
    ArgumentRules* maxArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "max", Real::getClassTypeSpec(), maxArgRules) );
    
    // add method for call "column" as a function
    ArgumentRules* columnArgRules = new ArgumentRules();
    columnArgRules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the column.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "column", ModelVector<Real>::getClassTypeSpec(), columnArgRules ) );
    
    // add method for call "size" as a function
    ArgumentRules* sizeArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "size", Natural::getClassTypeSpec(), sizeArgRules) );

    // add method for call "size" as a function
    ArgumentRules* nRowsArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "nrows", Natural::getClassTypeSpec(), nRowsArgRules) );

    // add method for call "size" as a function
    ArgumentRules* nColumnsArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "ncolumns", Natural::getClassTypeSpec(), nColumnsArgRules) );

}

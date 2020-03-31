#include <sstream>
#include <string>
#include <vector>

#include "ConstantNode.h"
#include "MatrixReal.h"
#include "Natural.h"
#include "ModelVector.h"
#include "RlStochasticMatrix.h"
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
#include "RlMatrixRealPos.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

/* Default constructor */
StochasticMatrix::StochasticMatrix(void) : MatrixRealPos()
{
    
    // initialize the member methods
    initializeMethods();
}

StochasticMatrix::StochasticMatrix(const RevBayesCore::MatrixReal& from) : MatrixRealPos( from )
{
    
    // initialize the member methods
    initializeMethods();
}

StochasticMatrix::StochasticMatrix(RevBayesCore::MatrixReal* m) : MatrixRealPos( m )
{
    
    // initialize the member methods
    initializeMethods();
}

StochasticMatrix::StochasticMatrix( RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal> * mat ) : MatrixRealPos( mat )
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
StochasticMatrix* StochasticMatrix::clone(void) const
{
    
    return new StochasticMatrix(*this);
}


/* Map calls to member methods */
RevPtr<RevVariable> StochasticMatrix::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
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
const std::string& StochasticMatrix::getClassType(void)
{
    
    static std::string rev_type = "StochasticMatrix";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& StochasticMatrix::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( MatrixReal::getClassTypeSpec() ) );
    
    return rev_type_spec;
}

/** Get type spec */
const TypeSpec& StochasticMatrix::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


void StochasticMatrix::initializeMethods( void )
{
    
    methods.eraseFunction( "[]" );
    methods.eraseFunction( "min" );
    methods.eraseFunction( "max" );
    methods.eraseFunction( "column" );
    
    // Add method for call "x[]" as a function
    ArgumentRules* squareBracketArgRules = new ArgumentRules();
    squareBracketArgRules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the row.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberFunction<StochasticMatrix, ModelVector<RealPos> >("[]", this, squareBracketArgRules ) );
    
    // add method for call "column" as a function
    ArgumentRules* columnArgRules = new ArgumentRules();
    columnArgRules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the column.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "column", ModelVector<RealPos>::getClassTypeSpec(), columnArgRules ) );
    
}



#include "ConstantNode.h"
#include "Integer.h"
#include "Natural.h"
#include "ModelVector.h"
#include "RlBoolean.h"
#include "Probability.h"
#include "RlMatrixRealSymmetric.h"
#include "RbUtil.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "RlMemberFunction.h"

#include <iomanip>
#include <sstream>

using namespace RevLanguage;

/* Default constructor */
MatrixRealSymmetric::MatrixRealSymmetric(void) : MatrixReal()
{
    
    // initialize the member methods
    initializeMethods();
}


/* Construct from reference */
MatrixRealSymmetric::MatrixRealSymmetric( const RevBayesCore::MatrixReal &mat ) : MatrixReal( mat.clone() )
{
    
    // initialize the member methods
    initializeMethods();
}


/* Construct from pointer */
MatrixRealSymmetric::MatrixRealSymmetric( RevBayesCore::MatrixReal *mat ) : MatrixReal( mat )
{
    
    // initialize the member methods
    initializeMethods();
}


/* Construct from DAG node */
MatrixRealSymmetric::MatrixRealSymmetric( RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal> * mat ) : MatrixReal( mat )
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
MatrixRealSymmetric* MatrixRealSymmetric::clone(void) const
{
    
	return new MatrixRealSymmetric(*this);
}


/** Convert to type. The caller manages the returned object. */
RevObject* MatrixRealSymmetric::convertTo( const TypeSpec& type ) const
{
    
    return RevObject::convertTo( type );
}

RevPtr<RevVariable> MatrixRealSymmetric::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found) {
    
//    if (name == "precision")
//    {
//        found = true;
//        
//        //int n = (int)this->dag_node->getValue().getNumberOfCharacters();
//
//        RevBayesCore::MatrixReal& m = this->dag_node->getValue();
//        std::cout << m[0][0] << std::endl;
//
//        // get the member with give index
//        /*const RevBayesCore::MatrixReal& m = static_cast<const RevBayesCore::MatrixReal&>( args[0].getVariable()->getRevObject() );
//        
//        if (this->dag_node->getValue().getNumberOfStates() < (size_t)(index.getValue()) ) {
//            throw RbException("Index out of bounds in []");
//            }
//        
//        const std::vector<double>& element = this->dag_node->getValue()[ size_t(index.getValue()) - 1];
//        RevBayesCore::RbVector<double> elementVector;
//        for (size_t i=0; i < this->dag_node->getValue().size(); ++i) {
//            elementVector.push_back( element[i] );
//            }
//        
//        return new RevVariable( new ModelVector<Real>( elementVector ) );*/
//    }
//    
//    if (name == "diagonal")
//    {
//        found = true;
//        RevBayesCore::RbVector<double> d = this->dag_node->getValue().getDiagonal();
//        return new RevVariable( new ModelVector<Real>( d ) );
//    }
    
    return ModelObject<RevBayesCore::MatrixReal>::executeMethod( name, args, found );
    
}

/** Get Rev type of object */
const std::string& MatrixRealSymmetric::getClassType(void)
{
    
    static std::string rev_type = "MatrixRealSymmetric";
    
	return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& MatrixRealSymmetric::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( MatrixReal::getClassTypeSpec() ) );
    
	return rev_type_spec;
}

/** Get type spec */
const TypeSpec& MatrixRealSymmetric::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

void MatrixRealSymmetric::initializeMethods( void )
{

    ArgumentRules* precisionArgRules = new ArgumentRules();
    methods.addFunction( new MemberFunction<MatrixRealSymmetric, MatrixRealSymmetric>( "precision", this, precisionArgRules) );
    
    ArgumentRules* diagonalArgRules = new ArgumentRules();
    methods.addFunction( new MemberFunction<MatrixRealSymmetric, ModelVector<RealPos> >( "diagonal", this, diagonalArgRules) );
    
    ArgumentRules* upperTriangleArgRules = new ArgumentRules();
    methods.addFunction( new MemberFunction<MatrixRealSymmetric, ModelVector<RealPos> >( "upperTriangle", this, upperTriangleArgRules) );
    

}



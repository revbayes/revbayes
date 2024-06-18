#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "AverageDistanceMatrix.h"
#include "DistanceMatrix.h"
#include "MatrixReal.h"
#include "ModelVector.h"
#include "MemberProcedure.h"
#include "Natural.h"
#include "RbVector.h"
#include "RlAverageDistanceMatrix.h"
#include "RlDistanceMatrix.h"
#include "RlMatrixBoolean.h"
#include "RlString.h"
#include "RevVariable.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MemberFunction.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "RbException.h"
#include "RbVectorImpl.h"
#include "Real.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RlConstantNode.h"
#include "RlMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlUtils.h"
#include "Taxon.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

/* Default constructor */
AverageDistanceMatrix::AverageDistanceMatrix(void) : ModelObject<RevBayesCore::AverageDistanceMatrix>( )
{
    
    // initialize the member methods
    initializeMethods();
}

AverageDistanceMatrix::AverageDistanceMatrix( RevBayesCore::AverageDistanceMatrix *a) : ModelObject<RevBayesCore::AverageDistanceMatrix>( a )
{
    
    // initialize the member methods
    initializeMethods();
}

AverageDistanceMatrix::AverageDistanceMatrix(const RevBayesCore::AverageDistanceMatrix &v) : ModelObject<RevBayesCore::AverageDistanceMatrix>( v.clone() )
{
    
    // initialize the member methods
    initializeMethods();
}

AverageDistanceMatrix::AverageDistanceMatrix( RevBayesCore::TypedDagNode<RevBayesCore::AverageDistanceMatrix> * mat ) : ModelObject<RevBayesCore::AverageDistanceMatrix>( mat )
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
AverageDistanceMatrix* AverageDistanceMatrix::clone() const
{
    
    return new AverageDistanceMatrix( *this );
}


/** Get Rev type of object */
const std::string& AverageDistanceMatrix::getClassType(void)
{
    
    static std::string rev_type = "AverageDistanceMatrix";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& AverageDistanceMatrix::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/** Get the type spec of this class. We return a member variable because instances might have different element types. */
const TypeSpec& AverageDistanceMatrix::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


size_t AverageDistanceMatrix::size( void ) const
{
    return this->dag_node->getValue().size();
}


/* Map calls to member methods */
RevPtr<RevVariable> AverageDistanceMatrix::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "distanceMatrix")
    {
        found = true;

        return new RevVariable(new DistanceMatrix((RevBayesCore::DistanceMatrix) (this->dag_node->getValue().getDistanceMatrix() ) ) ) ;
    }
    else if (name == "mask")
    {
        found = true;

        return new RevVariable(new MatrixBoolean((RevBayesCore::MatrixBoolean) (this->dag_node->getValue().getMask() ) ) ) ;

    }
    else if (name == "names")
    {
        found = true;
        
        std::vector<RevBayesCore::Taxon> names = this->dag_node->getValue().getTaxa();
        ModelVector<RlString> *n = new ModelVector<RlString>();
        for (size_t i = 0; i < names.size(); ++i)
        {
            n->push_back( names[i].getName() );
        }
        return new RevVariable( n );
    }
    else if (name == "size")
    {
        found = true;
        
        size_t s = this->dag_node->getValue().getSize();
        return new RevVariable( new Natural( s ) );
    }
    else if (name == "completeness")
    {
        found = true;
        
        double cmpl = this->dag_node->getValue().getCompleteness();
        return new RevVariable( new Real( cmpl ) );
    }
    else if ( name == "ultrametricImputation" )
    {
        found = true;
            
        this->dag_node->getValue().ultrametricImputation();
        return NULL;
    }
    else if ( name == "getElement" )
    {
        
        if ( args.size() > 1 && args[0].getVariable()->getRevObject().isType( Natural::getClassTypeSpec() ) && args[1].getVariable()->getRevObject().isType( Natural::getClassTypeSpec() ) )
        {
            found = true;
            // get the member with given indices
            const Natural& i = static_cast<const Natural&>( args[0].getVariable()->getRevObject() );
            const Natural& j = static_cast<const Natural&>( args[1].getVariable()->getRevObject() );

            if ( size() < (size_t)(i.getValue()) )
            {
                throw RbException("Index i out of bounds in getElement");
            }
            if ( size() < (size_t)(j.getValue()) )
            {
                throw RbException("Index j out of bounds in getElement");
            }

            std::pair<double, bool> element = static_cast< RevBayesCore::AverageDistanceMatrix& >( this->dag_node->getValue() ).getElement(size_t(i.getValue()) - 1, size_t(j.getValue()) - 1);
            
            ModelVector<RlString> *n = new ModelVector<RlString>();
            
            std::stringstream s;
            s << element.first << ", ";
            if (element.second)
            {
                s << "TRUE";
            }
            else
            {
                s << "FALSE";
            }
            n->push_back( s.str() );

            return new RevVariable( n );

        }

    }
    
    return ModelObject<RevBayesCore::AverageDistanceMatrix>::executeMethod( name, args, found );
}


void AverageDistanceMatrix::initializeMethods( void )
{
    // add method for call "distanceMatrix" as a function
    ArgumentRules* distanceMatrixArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "distanceMatrix", DistanceMatrix::getClassTypeSpec(), distanceMatrixArgRules ) );
    
    // add method for call "mask" as a function
    ArgumentRules* maskArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "mask", MatrixBoolean::getClassTypeSpec(), maskArgRules ) );
    
    // add method for call "names" as a function
    ArgumentRules* namesArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "names", ModelVector<RlString>::getClassTypeSpec(), namesArgRules ) );
    
    // add method for call "size" as a function
    ArgumentRules* sizeArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "size", Natural::getClassTypeSpec(), sizeArgRules ) );
    
    // add method for call "completeness" as a function
    ArgumentRules* complArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "completeness", Real::getClassTypeSpec(), complArgRules) );
    
    // add method for call "ultrametricImputation" as a function
    ArgumentRules* ultraImputeArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "ultrametricImputation", RlUtils::Void, ultraImputeArgRules) );
    
    // add method for call "getElement" as a function
    ArgumentRules* getElementArgRules = new ArgumentRules();
    getElementArgRules->push_back( new ArgumentRule( "i" , Natural::getClassTypeSpec(), "The row.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    getElementArgRules->push_back( new ArgumentRule( "j" , Natural::getClassTypeSpec(), "The column.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "getElement", ModelVector<Real>::getClassTypeSpec(), getElementArgRules) );

}


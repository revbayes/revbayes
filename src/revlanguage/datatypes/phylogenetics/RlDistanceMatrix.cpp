#include <cstddef>
#include <vector>
#include <iosfwd>
#include <string>

#include "ArgumentRule.h"
#include "DistanceMatrix.h"
#include "MatrixReal.h"
#include "ModelVector.h"
#include "MemberProcedure.h"
#include "Natural.h"
#include "RbVector.h"
#include "RlDistanceMatrix.h"
#include "RlMatrixReal.h"
#include "RlMatrixRealSymmetric.h"
#include "RlString.h"
#include "RevVariable.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "RbException.h"
#include "RbVectorImpl.h"
#include "Real.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RlConstantNode.h"
#include "RlUtils.h"
#include "Taxon.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

DistanceMatrix::DistanceMatrix(void) : ModelObject<RevBayesCore::DistanceMatrix>( )
{

    ArgumentRules* namesArgRules               = new ArgumentRules();
    ArgumentRules* matrixArgRules              = new ArgumentRules();
    ArgumentRules* symmetricMatrixArgRules     = new ArgumentRules();
    ArgumentRules* elementArgRules             = new ArgumentRules();
    ArgumentRules* setElementArgRules          = new ArgumentRules();
    ArgumentRules* sizeArgRules                = new ArgumentRules();
    
    elementArgRules->push_back( new ArgumentRule( "i" , Natural::getClassTypeSpec(), "The row.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    elementArgRules->push_back( new ArgumentRule( "j" , Natural::getClassTypeSpec(), "The column.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    
    setElementArgRules->push_back( new ArgumentRule( "i" , Natural::getClassTypeSpec(), "The row.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    setElementArgRules->push_back( new ArgumentRule( "j" , Natural::getClassTypeSpec(), "The column.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    setElementArgRules->push_back( new ArgumentRule( "v" , Real::getClassTypeSpec(), "The value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    
    methods.addFunction( new MemberProcedure( "names", ModelVector<RlString>::getClassTypeSpec(), namesArgRules           ) );
    methods.addFunction( new MemberProcedure( "matrix", MatrixReal::getClassTypeSpec(),            matrixArgRules          ) );
    methods.addFunction( new MemberProcedure( "symmetricMatrix", MatrixRealSymmetric::getClassTypeSpec(),            symmetricMatrixArgRules          ) );
    methods.addFunction( new MemberProcedure( "getElement", Real::getClassTypeSpec(),                  elementArgRules         ) );
    methods.addFunction( new MemberProcedure( "setElement", RlUtils::Void,                  		      setElementArgRules      ) );
    methods.addFunction( new MemberProcedure( "size", Natural::getClassTypeSpec(), sizeArgRules           ) );

	
}


DistanceMatrix::DistanceMatrix( RevBayesCore::DistanceMatrix *v) : ModelObject<RevBayesCore::DistanceMatrix>( v )
{

    ArgumentRules* namesArgRules               = new ArgumentRules();
    ArgumentRules* matrixArgRules              = new ArgumentRules();
    ArgumentRules* symmetricMatrixArgRules     = new ArgumentRules();
    ArgumentRules* elementArgRules             = new ArgumentRules();
    ArgumentRules* setElementArgRules          = new ArgumentRules();
    ArgumentRules* sizeArgRules                = new ArgumentRules();
    
    elementArgRules->push_back( new ArgumentRule( "i" , Natural::getClassTypeSpec(), "The row.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    elementArgRules->push_back( new ArgumentRule( "j" , Natural::getClassTypeSpec(), "The column.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    
    setElementArgRules->push_back( new ArgumentRule( "i" , Natural::getClassTypeSpec(), "The row.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    setElementArgRules->push_back( new ArgumentRule( "j" , Natural::getClassTypeSpec(), "The column.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    setElementArgRules->push_back( new ArgumentRule( "v" , Real::getClassTypeSpec(), "The value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    
    methods.addFunction( new MemberProcedure( "names", ModelVector<RlString>::getClassTypeSpec(), namesArgRules           ) );
    methods.addFunction( new MemberProcedure( "matrix", MatrixReal::getClassTypeSpec(),            matrixArgRules          ) );
    methods.addFunction( new MemberProcedure( "symmetricMatrix", MatrixRealSymmetric::getClassTypeSpec(),            symmetricMatrixArgRules          ) );
    methods.addFunction( new MemberProcedure( "getElement", Real::getClassTypeSpec(),                  elementArgRules         ) );
    methods.addFunction( new MemberProcedure( "setElement", RlUtils::Void,                  		      setElementArgRules      ) );
    methods.addFunction( new MemberProcedure( "size", Natural::getClassTypeSpec(), sizeArgRules           ) );


}


DistanceMatrix::DistanceMatrix( const RevBayesCore::DistanceMatrix &v) : ModelObject<RevBayesCore::DistanceMatrix>( v.clone() )
{
    
    ArgumentRules* namesArgRules               = new ArgumentRules();
    ArgumentRules* matrixArgRules              = new ArgumentRules();
    ArgumentRules* symmetricMatrixArgRules     = new ArgumentRules();
    ArgumentRules* elementArgRules             = new ArgumentRules();
    ArgumentRules* setElementArgRules          = new ArgumentRules();
    ArgumentRules* sizeArgRules                = new ArgumentRules();
    
    elementArgRules->push_back( new ArgumentRule( "i" , Natural::getClassTypeSpec(), "The row.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    elementArgRules->push_back( new ArgumentRule( "j" , Natural::getClassTypeSpec(), "The column.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    
    setElementArgRules->push_back( new ArgumentRule( "i" , Natural::getClassTypeSpec(), "The row.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    setElementArgRules->push_back( new ArgumentRule( "j" , Natural::getClassTypeSpec(), "The column.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    setElementArgRules->push_back( new ArgumentRule( "v" , Real::getClassTypeSpec(), "The value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    
    methods.addFunction( new MemberProcedure( "names", ModelVector<RlString>::getClassTypeSpec(), namesArgRules           ) );
    methods.addFunction( new MemberProcedure( "matrix", MatrixReal::getClassTypeSpec(),            matrixArgRules          ) );
    methods.addFunction( new MemberProcedure( "symmetricMatrix", MatrixRealSymmetric::getClassTypeSpec(),            symmetricMatrixArgRules          ) );
    methods.addFunction( new MemberProcedure( "getElement", Real::getClassTypeSpec(),                  elementArgRules         ) );
    methods.addFunction( new MemberProcedure( "setElement", RlUtils::Void,                  		      setElementArgRules      ) );
    methods.addFunction( new MemberProcedure( "size", Natural::getClassTypeSpec(), sizeArgRules           ) );

}


DistanceMatrix::DistanceMatrix( RevBayesCore::TypedDagNode<RevBayesCore::DistanceMatrix> *m) : ModelObject<RevBayesCore::DistanceMatrix>( m )
{

	ArgumentRules* namesArgRules               = new ArgumentRules();
	ArgumentRules* matrixArgRules              = new ArgumentRules();
    ArgumentRules* symmetricMatrixArgRules     = new ArgumentRules();
	ArgumentRules* elementArgRules             = new ArgumentRules();
	ArgumentRules* setElementArgRules          = new ArgumentRules();
    ArgumentRules* sizeArgRules                = new ArgumentRules();

    elementArgRules->push_back( new ArgumentRule( "i" , Natural::getClassTypeSpec(), "The row.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    elementArgRules->push_back( new ArgumentRule( "j" , Natural::getClassTypeSpec(), "The column.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    
    setElementArgRules->push_back( new ArgumentRule( "i" , Natural::getClassTypeSpec(), "The row.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    setElementArgRules->push_back( new ArgumentRule( "j" , Natural::getClassTypeSpec(), "The column.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    setElementArgRules->push_back( new ArgumentRule( "v" , Real::getClassTypeSpec(), "The value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    
    methods.addFunction( new MemberProcedure( "names", ModelVector<RlString>::getClassTypeSpec(), namesArgRules           ) );
    methods.addFunction( new MemberProcedure( "matrix", MatrixReal::getClassTypeSpec(),            matrixArgRules          ) );
    methods.addFunction( new MemberProcedure( "symmetricMatrix", MatrixRealSymmetric::getClassTypeSpec(),            symmetricMatrixArgRules          ) );
    methods.addFunction( new MemberProcedure( "getElement", Real::getClassTypeSpec(),                  elementArgRules         ) );
    methods.addFunction( new MemberProcedure( "setElement", RlUtils::Void,                  		      setElementArgRules      ) );
    methods.addFunction( new MemberProcedure( "size", Natural::getClassTypeSpec(), sizeArgRules           ) );

}


DistanceMatrix* DistanceMatrix::clone() const
{
    return new DistanceMatrix( *this );
}


/* Map calls to member methods */
RevPtr<RevVariable> DistanceMatrix::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "matrix")
    {
        found = true;

        return new RevVariable(new MatrixReal((RevBayesCore::MatrixReal) (this->dag_node->getValue().getMatrix() ) ) ) ;
    }
    else if (name == "symmetricMatrix")
    {
        found = true;
        //TODO: make sure the matrix really is symmetric
        return new RevVariable(new MatrixRealSymmetric((RevBayesCore::MatrixReal) (this->dag_node->getValue().getMatrix() ) ) ) ;

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
	else if ( name == "size" )
	{
		found = true;
		
		// return a new RevVariable with the size of this container
		return RevPtr<RevVariable>( new RevVariable( new Natural( size() ), "" ) );
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

			double element = static_cast< RevBayesCore::DistanceMatrix& >( this->dag_node->getValue() ).getElement(size_t(i.getValue()) - 1, size_t(j.getValue()) - 1);
			
			return new RevVariable( new Real( element ) );

		}

	}
	else if ( name == "setElement" )
	{
		
		if ( args.size() > 2 && args[0].getVariable()->getRevObject().isType( Natural::getClassTypeSpec() ) && args[1].getVariable()->getRevObject().isType( Natural::getClassTypeSpec() ) && args[2].getVariable()->getRevObject().isType( Real::getClassTypeSpec() ) )
		{
			found = true;
			// get the member with given indices
			const Natural& i = static_cast<const Natural&>( args[0].getVariable()->getRevObject() );
			const Natural& j = static_cast<const Natural&>( args[1].getVariable()->getRevObject() );
			const Real& v = static_cast<const Real&>( args[2].getVariable()->getRevObject() );

			if ( size() < (size_t)(i.getValue()) )
			{
				throw RbException("Index i out of bounds in getElement");
			}
			if ( size() < (size_t)(j.getValue()) )
			{
				throw RbException("Index j out of bounds in getElement");
			}
			
			static_cast< RevBayesCore::DistanceMatrix& >( this->dag_node->getValue() ).getElement(size_t(i.getValue()) - 1, size_t(j.getValue()) - 1) = double(v.getValue() );
			
		}
		return NULL;
		
	}

	
    return ModelObject<RevBayesCore::DistanceMatrix>::executeMethod( name, args, found );
}





const Real* DistanceMatrix::getElement(size_t idx, size_t idy) const
{
	double element = static_cast< RevBayesCore::DistanceMatrix& >( this->dag_node->getValue() ).getElement(idx - 1, idy - 1);
	
	return new Real( element ) ;

}


void DistanceMatrix::setElement(size_t idx, size_t idy, double& value) {

	static_cast< RevBayesCore::DistanceMatrix& >( this->dag_node->getValue() ).getElement(idx - 1, idy - 1) = value;
	
	return ;
}


/*
template <typename rlType>
rlType* ModelVector<rlType>::getElement(size_t idx) const
{
	return new rlType( this->getValue()[ idx ] );
}



 RevObject*   RlDistanceMatrix::getElement(size_t idx) const {
	
	return new RevVariable( this->getValue()[ idx ] );
}*/

/**
 * Add an element to the end of the vector.
 */

/*void RlDistanceMatrix::push_back(const RevObject &x)
{
	//We do nothing, in principle we could implement a push_back...
}
*/

/**
 * Size of the matrix.
 */
size_t DistanceMatrix::size( void ) const
{
	return this->dag_node->getValue().size();
}




/* Get Rev type of object */
const std::string& DistanceMatrix::getClassType(void)
{
    
    static std::string rev_type = "DistanceMatrix";
    
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& DistanceMatrix::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/** Get the type spec of this class. We return a member variable because instances might have different element types. */
const TypeSpec& DistanceMatrix::getTypeSpec(void) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}


#include "RlSimplex.h"

#include "ConstantNode.h"
#include "ModelVector.h"
#include "RlMemberFunction.h"
#include "Natural.h"
#include "Probability.h"
#include "RbException.h"
#include "Real.h"
#include "RealPos.h"
#include "WorkspaceVector.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Ellipsis.h"
#include "MemberFunction.h"
#include "MethodTable.h"
#include "RbConstIterator.h"
#include "RbConstIteratorImpl.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevVariable.h"
#include "RlDeterministicNode.h"
#include "RlTypedFunction.h"
#include "SimplexFunction.h"
#include "SimplexFromVectorFunction.h"
#include "TypeSpec.h"
#include "WorkspaceToCoreWrapperObject.h"

namespace RevLanguage { class Argument; }

using namespace RevLanguage;


/**
 * Construct an empty simplex.
 * We simply rely on the base class.
 */
Simplex::Simplex( void ) : ModelObject<RevBayesCore::Simplex>( new RevBayesCore::Simplex() )
{
    initMethods();
}


/**
 * Construct simplex from a vector of real values (doubles). We
 * rescale the values here and make sure they are all positive.
 * Just in case.
 */
Simplex::Simplex( const RevBayesCore::Simplex& v ) : ModelObject<RevBayesCore::Simplex>()
{
    
    // Now set the constant value of the simplex
    RevBayesCore::ConstantNode< RevBayesCore::Simplex >* newNode = new RevBayesCore::ConstantNode< RevBayesCore::Simplex >( "", v.clone() );
    this->setDagNode( newNode );
    
    initMethods();
}


/**
 * Make sure that a distribution or a function associated with
 * a dynamic node are of the right type.
 *
 * @todo Make sure we actually have a simplex stored in n (or an
 *       NA value)
 */
Simplex::Simplex( RevBayesCore::TypedDagNode<RevBayesCore::Simplex>* n ) : ModelObject<RevBayesCore::Simplex>( n )
{

    initMethods();
}


/**
 * Virtual destructor to ensure proper handling
 * of deletion from base class pointers or
 * references in case of future derived classes.
 */
Simplex::~Simplex()
{
}


/** Generate a type-safe clone */
Simplex* Simplex::clone( void ) const
{
    return new Simplex( *this );
}


void Simplex::constructInternalObject( void )
{
  // we free the memory first
    if ( dag_node != NULL )
    {
        if ( dag_node->decrementReferenceCount() == 0 )
        {
            delete dag_node;
        }
    }

    RevBayesCore::TypedFunction<RevBayesCore::Simplex> * f = NULL;

    if ( args.size() == 1 )
    {
        RevObject& b = args[0]->getRevObject();

        if ( b.isType( ModelVector<RealPos>::getClassTypeSpec() ) )
        {
            const RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* vec;
            vec = static_cast< const ModelVector<RealPos>& >( args[0]->getRevObject() ).getDagNode();

            f = new RevBayesCore::SimplexFromVectorFunction( vec );
        }
        else
        {
            throw( RbException("Must provide at least 2 RealPos values to Simplex constructor") );
        }
    }
    else if ( args.size() > 1 )
    {
        std::vector<const RevBayesCore::TypedDagNode<double>* > params;

        for ( size_t i = 0; i < args.size(); i++ )
        {
            RevObject& b = args[i]->getRevObject();
            if ( b.isType( RealPos::getClassTypeSpec() ) )
            {
                const RealPos &val = static_cast<const RealPos &>( args[i]->getRevObject() );
                params.push_back( val.getDagNode() );
            }
            else
            {
                throw( RbException("Can only provide one RealPos[] vector to Simplex constructor") );
            }
        }

        f = new RevBayesCore::SimplexFunction( params );
    }
    else
    {
        throw( RbException("Must provide at least 2 RealPos values or exactly one RealPos[] vector to Simplex constructor") );
    }

    dag_node = new RevBayesCore::DeterministicNode<RevBayesCore::Simplex>("", f);
    dag_node->incrementReferenceCount();

}


/**
 * In this function we check whether this type is convertible to some other
 * Rev object type. Here we focus entirely on supporting conversion to
 * other generic vectors with compatible elements. This is not done automatically
 * because of the templating: a vector of RealPos does not inherit from a vector
 * of Real, for example.
 */
RevObject* Simplex::convertTo( const TypeSpec& type ) const
{
    
    if ( type == ModelVector<Real>::getClassTypeSpec() )
    {
        return new ModelVector<Real>( RevBayesCore::RbVector<double>(dag_node->getValue()) );
    }
    else if ( type == ModelVector<RealPos>::getClassTypeSpec() )
    {
        return new ModelVector<RealPos>( RevBayesCore::RbVector<double>(dag_node->getValue()) );
    }

    
    return ModelObject<RevBayesCore::Simplex>::convertTo( type );
}


/**
 * Map calls to member methods.
 */
RevPtr<RevVariable> Simplex::executeMethod( std::string const &name, const std::vector<Argument> &args, bool &found )
{
    
    if ( name == "size" )
    {
        found = true;
        
        // return a new RevVariable with the size of this container
        return RevPtr<RevVariable>( new RevVariable( new Natural( size() ), "" ) );
    }
    
    return ModelObject<RevBayesCore::Simplex>::executeMethod( name, args, found );
}



/** Get Rev type of object */
const std::string& Simplex::getClassType(void)
{
    
    static std::string rev_type = "Simplex";
    
	return rev_type; 
}


/** Get class type spec describing type of object */
const TypeSpec& Simplex::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &ModelObject<RevBayesCore::Simplex>::getClassTypeSpec(), &Probability::getClassTypeSpec() );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Simplex::getConstructorFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "Simplex";

    return f_name;
}


/**
 * Get the aliases of the Rev name for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Simplex::getConstructorFunctionAliases( void ) const
{
    // create a constructor function name alias variable that is the same for all instance of this class
    std::vector<std::string> aliases;

    aliases.push_back("simplex");

    return aliases;
}


/** Return member rules (no members) */
const MemberRules& Simplex::getParameterRules(void) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "x", ModelVector<RealPos>::getClassTypeSpec(), "A vector of values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        argumentRules.push_back( new ArgumentRule( "x1", RealPos::getClassTypeSpec(), "first value", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        argumentRules.push_back( new ArgumentRule( "x2", RealPos::getClassTypeSpec(), "second value", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        argumentRules.push_back( new Ellipsis ( "additional values", RealPos::getClassTypeSpec() ) );
        rules_set = true;
    }

    return argumentRules;
}



/** Get the type spec of this class */
const TypeSpec& Simplex::getTypeSpec(void) const
{
    return getClassTypeSpec();
}


Probability* Simplex::getElement(size_t idx) const
{
    return new Probability( this->getValue()[ idx ] );
}


/**
 * Initialize the methods.
 */
void Simplex::initMethods( void )
{
	ArgumentRules* sizeArgRules = new ArgumentRules();
	this->methods.addFunction( new MemberProcedure( "size", Natural::getClassTypeSpec(), sizeArgRules) );

}


/**
 * In this function we check whether this type is convertible to some other
 * Rev object type. Here we focus entirely on supporting conversion to
 * other generic vectors with compatible elements. This is not done automatically
 * because of the templating: a vector of RealPos does not inherit from a vector
 * of Real, for example.
 */
double Simplex::isConvertibleTo( const TypeSpec& type, bool once ) const
{
    
    if ( once == true && type.getParentType() == getClassTypeSpec().getParentType() )
    {
        // We want to convert to another model vector
        
        // Simply check whether our elements can convert to the desired element type
        RevBayesCore::RbConstIterator<double> i;
        double penalty = 0.0;
        for ( i = this->getValue().begin(); i != this->getValue().end(); ++i )
        {
            const double& orgInternalElement = *i;
            Probability orgElement = Probability( orgInternalElement );
            
            // Test whether this element is already of the desired element type or can be converted to it
            if ( type.getElementTypeSpec() != NULL && orgElement.getTypeSpec() != *type.getElementTypeSpec() )
            {
                
                double element_penalty = orgElement.isConvertibleTo( *type.getElementTypeSpec(), once );
                if ( element_penalty == -1 )
                {
                    // we cannot convert this element
                    return -1;
                }
                penalty += element_penalty;
            }
            
        }
        
        return penalty;
    }
    else if ( type == WorkspaceVector<RevObject>::getClassTypeSpec() )
    {
        // yes we can
        return 0.0;
    }
    
    return ModelObject<RevBayesCore::Simplex>::isConvertibleTo( type, once );
}


/**
 * Add an element to the end of the vector.
 */
void Simplex::push_back(const double &x)
{
    return this->dag_node->getValue().push_back( x );
}


/**
 * Add an element to the end of the vector.
 */
void Simplex::push_back(const Probability &x)
{
    return this->dag_node->getValue().push_back( x.getValue() );
}


/**
 * Push a Rev object element onto the back of the vector.
 */
void Simplex::push_back( const RevObject &x )
{
    
    // cast the object
    const Probability *x_converted = dynamic_cast< const Probability* >( &x );
    
    if ( x_converted == NULL )
    {
        throw RbException() << "Could not append an element of type " << x.getType() << " to a vector of type " <<  this->getType() ;
    }
    
    // Push it onto the back of the elements vector
    this->push_back( *x_converted );
}


/**
 * Print the value of the vector, respecting the formatting of the
 * model object elements. We do this by retrieving the elements, one
 * by one, and printing them using their own printValue implementation.
 * Among other things, this takes care of proper formatting.
 *
 * We make a perfectly safe const cast here, since we only utilize the
 * const printValue function of the element.
 */
void Simplex::printValue( std::ostream& o, bool user ) const
{
    
    this->getDagNode()->printValue( o, ",", -1, true, user, true );
    
}


/**
 * Size of the vector.
 */
size_t Simplex::size( void ) const
{
    return this->dag_node->getValue().size();
}


/** Throw an error if someone attempts to sort us */
void Simplex::sort( void )
{
    throw RbException( "Illegal attempt to sort a Simplex variable" );
}


/** Throw an error if someone attempts to apply unique */
void Simplex::unique( void )
{
    throw RbException( "Illegal attempt to apply unique to a Simplex variable" );
}

/** Set a member variable */
void Simplex::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "x" || name == "x1" || name == "x2" || name == "additional values" || name == "" )
    {
        if ( var->getRevObject() != RevNullObject::getInstance() )
        {
            args.push_back(var);
        }
    }
    else
    {
        RevObject::setConstParameter(name, var);
    }
}


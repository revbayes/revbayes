#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <sstream> // IWYU pragma: keep

#include "ArgumentRules.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "OptionRule.h"
#include "RbException.h"
#include "RbHelpType.h"
#include "RevObject.h"
#include "RlUtils.h"
#include "TypeSpec.h"
#include "ArgumentRule.h"
#include "RbHelpArgument.h"
#include "RbHelpEntry.h"
#include "RbHelpFunction.h"
#include "RevMemberObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"

namespace RevBayesCore { class DagNode; }
namespace RevLanguage { class Argument; }
namespace RevLanguage { class UserFunction; }

using namespace RevLanguage;


/**
 * Default constructor.
 */
RevObject::RevObject( bool includeMemberMethods ) : RevMemberObject()
{
    
    if ( includeMemberMethods == true )
    {
        ArgumentRules* getMethodsArgRules = new ArgumentRules();
    
        // Add the 'methods()' method
        methods.addFunction( new MemberProcedure( "methods", RlUtils::Void, getMethodsArgRules) );
    }
    
}


/** 
 * Destructor. 
 * Nothing to do here because we don't own anything.
 */
RevObject::~RevObject( void )
{
    
}


/**
 * Addition operation.
 * Since we don't know the types and thus don't know the special behavior we simply throw and error.
 *
 * \param[in]   rhs     The right hand side value of the division operation.
 *
 * \return              NULL
 */
RevObject* RevObject::add(const RevObject &rhs) const
{
    throw RbException() << "Cannot add a value of type '" << this->getType() << "' to a value of type '" << rhs.getType() << "'.";
    
    return NULL;
}


void RevObject::addMethods(const MethodTable &mt)
{
    methods.insertInheritedMethods( mt );
}


/** The default implementation does nothing because we don't have an internal object */
void RevObject::constructInternalObject( void )
{
    // nothing to do
}



/** Convert to type and dim. The caller manages the returned object. */
RevObject* RevObject::convertTo(const TypeSpec& typeSpec) const
{
        
    throw RbException() << "Failed conversion from type '" << getType() << "' to type '" << typeSpec.getType() << "'" ;
    
    return NULL;
}


/**
 * Decrement operation.
 * Since we don't know the types and thus don't know the special behavior we simply throw and error.
 */
void RevObject::decrement( void ) 
{
    throw RbException() << "Cannot decrement a value of type '" << this->getType() << "'.";

}


/**
 * Division operation.
 * Since we don't know the types and thus don't know the special behavior we simply throw and error.
 *
 * \param[in]   rhs     The right hand side value of the division operation.
 *
 * \return              NULL
 */
RevObject* RevObject::divide(const RevObject &rhs) const
{
    throw RbException() << "Cannot divide a value of type '" << this->getType() << "' to a value of type '" << rhs.getType() << "'.";
    
    return NULL;
}


/** 
 * Execute simple method. 
 */
RevPtr<RevVariable> RevObject::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if ( name == "methods" )
    {
        found = true;
        
        // just print the method names (including inherited methods)
        const MethodTable &m = getMethods();
        m.printValue(std::cout, true);
        
        return NULL;
    }
    
    found = false;
    
    return NULL;
}


/** Get class vector describing type of object */
const std::string& RevObject::getClassType(void)
{
    
    static std::string rev_type = "RevObject";
	return rev_type; 
}


/**
 * Get class vector describing type of object. We are
 * the base class of all other classes, so parent is NULL.
 */
const TypeSpec& RevObject::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), NULL );
	
    return rev_type_spec; 
}



/** Return member rules (no members) */
const MemberRules& RevObject::getParameterRules(void) const
{
    
    static const MemberRules rules = MemberRules();
    
    return rules;
}


/**
 * Get common member methods.
 */
const MethodTable& RevObject::getMethods( void ) const
{
    
    return methods;
}


/**
 * Get the Rev type of the object
 */
const std::string& RevObject::getType( void ) const
{
    
    return getTypeSpec().getType();
}


/**
 * Get the value as a DAG node. This default implementation throws an error.
 */
RevBayesCore::DagNode* RevObject::getDagNode( void ) const
{
    throw RbException() << "Workspace objects cannot be used inside DAG's! You tried to access the DAG node of type '" << getClassType() << "'.";
    
    return NULL;
}


/**
 * Get the value as a DAG node. This default implementation throws an error.
 */
bool RevObject::hasDagNode( void ) const
{
    return false;
}


RevBayesCore::RbHelpEntry* RevObject::constructTypeSpecificHelp( void ) const
{
    
    return new RevBayesCore::RbHelpType();
}


void RevObject::addSpecificHelpFields(RevBayesCore::RbHelpEntry *e) const
{
    
    RevBayesCore::RbHelpType *help_entry = static_cast<RevBayesCore::RbHelpType*>( e );
    
    // create the constructor
    RevBayesCore::RbHelpFunction help_constructor = RevBayesCore::RbHelpFunction();
    
    // usage
    help_constructor.setUsage( getConstructorUsage() );
    
    // arguments
    const MemberRules& rules = getParameterRules();
    std::vector<RevBayesCore::RbHelpArgument> arguments = std::vector<RevBayesCore::RbHelpArgument>();
    
    for ( size_t i=0; i<rules.size(); ++i )
    {
        const ArgumentRule &the_rule = rules[i];
        
        RevBayesCore::RbHelpArgument argument = RevBayesCore::RbHelpArgument();
        
        argument.setLabel( the_rule.getArgumentLabel() );
        argument.setDescription( the_rule.getArgumentDescription() );
        
        std::string type = "<any>";
        if ( the_rule.getArgumentDagNodeType() == ArgumentRule::CONSTANT )
        {
            type = "<constant>";
        }
        else if ( the_rule.getArgumentDagNodeType() == ArgumentRule::STOCHASTIC )
        {
            type = "<stochastic>";
        }
        else if ( the_rule.getArgumentDagNodeType() == ArgumentRule::DETERMINISTIC )
        {
            type = "<deterministic>";
        }
        argument.setArgumentDagNodeType( type );
        
        std::string passing_method = "value";
        if ( the_rule.getEvaluationType() == ArgumentRule::BY_CONSTANT_REFERENCE )
        {
            passing_method = "const reference";
        }
        else if ( the_rule.getEvaluationType() == ArgumentRule::BY_REFERENCE )
        {
            passing_method = "reference";
        }
        argument.setArgumentPassingMethod(  passing_method );
        
        argument.setValueType( the_rule.getArgumentTypeSpec()[0].getType() );
        
        if ( the_rule.hasDefault() )
        {
            std::stringstream ss;
            the_rule.getDefaultVariable().getRevObject().printValue( ss, true);
            argument.setDefaultValue( ss.str() );
        }
        else
        {
            argument.setDefaultValue( "" );
        }
        
        // loop options
        std::vector<std::string> options = std::vector<std::string>();
        const OptionRule *opt_rule = dynamic_cast<const OptionRule*>( &the_rule );
        if ( opt_rule != NULL )
        {
            options = opt_rule->getOptions();
        }
        argument.setOptions( options );
        
        // add the argument to the argument list
        arguments.push_back( argument );
    }
    
    help_constructor.setArguments( arguments );
    
    //
    std::vector<RevBayesCore::RbHelpFunction> constructors;
    constructors.push_back( help_constructor );
    help_entry->setConstructors( constructors );
    
    
    help_entry->setMethods( getHelpMethods() );

}




/**
 * Get usage information for help system
 */
std::string RevObject::getConstructorUsage( void ) const
{
    
    std::ostringstream o;
    
    if ( getConstructorFunctionName() == "" )
    {
        o << "<unnamed>(";
    }
    else
    {
        o << "" << getConstructorFunctionName() << "(";
    }
    
    const ArgumentRules& argRules = getParameterRules();
    for (size_t i=0; i<argRules.size(); ++i)
    {
        if (i != 0)
        {
            o << ", ";
        }
        argRules[i].printValue(o);
    }
    o << ")";
    
    return o.str();
}


/**
 * Get the help entry for this class
 */
RevBayesCore::RbHelpEntry* RevObject::getHelpEntry( void ) const
{
    // create the help function entry that we will fill with some values
    RevBayesCore::RbHelpEntry *help = constructTypeSpecificHelp();
    RevBayesCore::RbHelpEntry &help_entry = *help;
    
    const TypeSpec& type_spec = getTypeSpec();

    // name
    std::string name = getConstructorFunctionName();
    if ( name == "c_name" )
    {
        name = type_spec.getType();
    }
    help_entry.setName( name );
    
    // aliases
    std::vector<std::string> aliases = getConstructorFunctionAliases();
    help_entry.setAliases( aliases );

    const TypeSpec* parentTypeSpec = type_spec.getParentTypeSpec();

    std::vector<std::string> typeSpec;
    while(parentTypeSpec != NULL)
    {
        typeSpec.push_back( parentTypeSpec->getType() );
        parentTypeSpec = parentTypeSpec->getParentTypeSpec();
    }

    // type spec
    help_entry.setTypeSpec( typeSpec );

    // now add the specific help stuff
    addSpecificHelpFields( help );
    
    return help;
    
}


/** Get the help entry for this class */
std::vector<RevBayesCore::RbHelpFunction> RevObject::getHelpMethods( void ) const
{
    
    // construct the vector
    std::vector<RevBayesCore::RbHelpFunction> help_methods;
    
    const MethodTable &methods = getMethods();
    
    for (std::multimap<std::string, Function*>::const_iterator it = methods.begin(); it != methods.end(); ++it)
    {
        
        Function &the_function = *it->second;
        
        // create the method
        RevBayesCore::RbHelpFunction help_method = RevBayesCore::RbHelpFunction();
        
        // name
        help_method.setName( it->first );
        
        // usage
        help_method.setUsage( the_function.getHelpUsage() );
        
        // arguments
        const ArgumentRules& rules = the_function.getArgumentRules();
        std::vector<RevBayesCore::RbHelpArgument> arguments = std::vector<RevBayesCore::RbHelpArgument>();
        
        for ( size_t i=0; i<rules.size(); ++i )
        {
            const ArgumentRule &the_rule = rules[i];
            
            RevBayesCore::RbHelpArgument argument = RevBayesCore::RbHelpArgument();
            
            argument.setLabel( the_rule.getArgumentLabel() );
            argument.setDescription( the_rule.getArgumentDescription() );
            
            std::string type = "<any>";
            if ( the_rule.getArgumentDagNodeType() == ArgumentRule::CONSTANT )
            {
                type = "<constant>";
            }
            else if ( the_rule.getArgumentDagNodeType() == ArgumentRule::STOCHASTIC )
            {
                type = "<stochastic>";
            }
            else if ( the_rule.getArgumentDagNodeType() == ArgumentRule::DETERMINISTIC )
            {
                type = "<deterministic>";
            }
            argument.setArgumentDagNodeType( type );
            
            std::string passing_method = "value";
            if ( the_rule.getEvaluationType() == ArgumentRule::BY_CONSTANT_REFERENCE )
            {
                passing_method = "const reference";
            }
            else if ( the_rule.getEvaluationType() == ArgumentRule::BY_REFERENCE )
            {
                passing_method = "reference";
            }
            argument.setArgumentPassingMethod(  passing_method );
            
            argument.setValueType( the_rule.getArgumentTypeSpec()[0].getType() );
            
            if ( the_rule.hasDefault() )
            {
                std::stringstream ss;
                the_rule.getDefaultVariable().getRevObject().printValue( ss, true);
                argument.setDefaultValue( ss.str() );
            }
            else
            {
                argument.setDefaultValue( "" );
            }
            
            // loop options
            std::vector<std::string> options = std::vector<std::string>();
            const OptionRule *opt_rule = dynamic_cast<const OptionRule*>( &the_rule );
            if ( opt_rule != NULL )
            {
                options = opt_rule->getOptions();
            }
            argument.setOptions( options );
            
            // add the argument to the argument list
            arguments.push_back( argument );
        }
        
        help_method.setArguments( arguments );
        
        // return value
        help_method.setReturnType( the_function.getReturnType().getType() );
        
        //
        help_methods.push_back( help_method );
        
    }
    
    // virtual dummy function
    
    return help_methods;
}


/**
 * Increment operation.
 * Since we don't know the types and thus don't know the special behavior we simply throw and error.
 */
void RevObject::increment( void ) 
{
    
    throw RbException() << "Cannot increment a value of type '" << this->getType() << "'.";
    
}


/**
 * The default implementation is that the object is not an abstract type object.
 * Only abstract type objects need to override this.
 */
bool RevObject::isAbstract( void ) const
{
    return false;
}


/**
 * Is the object or any of its upstream members or elements
 * modifiable by the user through assignment? By default we
 * return true, which is the safest option if a derived class
 * forgets to override this function.
 */
bool RevObject::isAssignable( void ) const
{
    return true;
}


/**
 * The default implementation is that the variable is constant. Only variables which actually store
 * internally DAG nodes have to ask the DAG nodes if they are constant.
 */
bool RevObject::isConstant( void ) const
{

    return true;
}


/** 
 * Is convertible to type? 
 * -1 represent false and any positive number the cost of conversion.
 */
double RevObject::isConvertibleTo(const TypeSpec& type, bool once) const
{
    
    return -1.0;
}


/* Are we of specified language type? */
bool RevObject::isType(const TypeSpec& typeSpec) const
{
    
    return getTypeSpec().isDerivedOf( typeSpec );
}


/**
 * The default implementation is that the variable is not a model object. Only variables which actually store
 * internally DAG nodes are model objects.
 */
bool RevObject::isModelObject( void ) const
{
    
    return false;
}


/**
  * Converting the value of the internal variable to a constant. The default implementation does nothing because we don't have a DAG node as our internal variable.
  * Note, RevLanguage types which can be used as types in the DAG should overwrite this method.
  */
void RevObject::makeConstantValue( void )
{
    // do nothing

}


/**
 * Make a new object that is an indirect deterministic reference to the object.
 * The default implementation throws an error.
 */
RevObject* RevObject::makeIndirectReference(void)
{
    std::ostringstream msg;
    msg << "The type '" << getClassType() << "' not supported in indirect reference assignments (yet)";
    throw RbException(msg.str());
}


/**
 * Convert the object to a deterministic object with a user-defined Rev function inside it.
 */
void RevObject::makeUserFunctionValue( UserFunction* fxn )
{
    std::ostringstream msg;
    msg << "The type '" << getClassType() << "' not supported in user-defined function nodes (yet)";
    throw RbException( msg.str() );
}


/**
 * Multiplication operation.
 * Since we don't know the types and thus don't know the special behavior we simply throw and error.
 *
 * \param[in]   rhs     The right hand side value of the division operation.
 *
 * \return              NULL
 */
RevObject* RevObject::multiply(const RevObject &rhs) const
{
    throw RbException() << "Cannot multiply a value of type '" << this->getType() << "' to a value of type '" << rhs.getType() << "'.";
    
    return NULL;
}


/** Set a member variable */
void RevObject::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    throw RbException() << "No parameter with name \"" << name << "\" found to set.";
}


/* Set a member variable.
 * In this default implementation, we delegate to setConstParameter.
 * Derived classes of MemberObject who need non-const RevVariable should overwrite this function.
 * If you don't care if the variable is const, then you should only overwrite the setConstParameter.
 */
void RevObject::setParameter(const std::string& name, const RevPtr<RevVariable> &var)
{
    
    setConstParameter(name, RevPtr<const RevVariable>( var ) );
    
}



/** 
 * Setting the name of the internal variable. The default implementation does nothing.
 * Note: Rev types that can be used as types in the DAG should override this method.
 */
void RevObject::setName(std::string const &n)
{
}


/**
 * Subtraction operation.
 * Since we don't know the types and thus don't know the special behavior we simply throw and error.
 *
 * \param[in]   rhs     The right hand side value of the division operation.
 *
 * \return              NULL
 */
RevObject* RevObject::subtract(const RevObject &rhs) const
{
    throw RbException() << "Cannot subtract a value of type '" << rhs.getType() << "' from a value of type '" << this->getType() << "'.";
    
    return NULL;
}


/**
 * Get the info about this object as a string.
 * By default we simply call the overloaded operator<<.
 *
 */
std::string RevObject::toString( void ) const
{
    std::stringstream o;
    printValue(o,true);
    
    return o.str();
}

/** Make sure we can print the value of the object easily */
std::ostream& operator<<(std::ostream& o, const RevObject& x)
{
    
    x.printValue(o,true);
    return o;
}

/** Make sure we can print the value of the object easily */
std::ostream& RevLanguage::operator<<(std::ostream& o, const RevObject& x)
{
    
    x.printValue(o,true);
    return o;
}


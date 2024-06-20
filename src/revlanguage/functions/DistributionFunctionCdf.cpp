#include "DistributionFunctionCdf.h"

#include <cstddef>

#include "ArgumentRule.h"
#include "CumulativeDistributionFunction.h"
#include "DeterministicNode.h"
#include "Probability.h"
#include "Real.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "RevVariable.h"
#include "RlContinuousDistribution.h"
#include "RlFunction.h"
#include "RlPositiveContinuousDistribution.h"

namespace RevBayesCore { class ContinuousDistribution; }

using namespace RevLanguage;

/** Constructor */
DistributionFunctionCdf::DistributionFunctionCdf( ContinuousDistribution *d ) : TypedFunction<Probability>(),
    templateObject( d ),
    templateObjectPositive( NULL )
{
    
    argRules.push_back( new ArgumentRule("x", Real::getClassTypeSpec(), "The value for which to compute the probability.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    const ArgumentRules &memberRules = templateObject->getParameterRules();
    for (std::vector<ArgumentRule*>::const_iterator it = memberRules.begin(); it != memberRules.end(); ++it)
    {
        argRules.push_back( (*it)->clone() );
    }
    
}

/** Constructor */
DistributionFunctionCdf::DistributionFunctionCdf( PositiveContinuousDistribution *d ) : TypedFunction<Probability>(),
    templateObject( NULL ),
    templateObjectPositive( d )
{
    
    argRules.push_back( new ArgumentRule("x", Real::getClassTypeSpec(), "The value for which to compute the probability.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    const ArgumentRules &memberRules = templateObjectPositive->getParameterRules();
    for (std::vector<ArgumentRule*>::const_iterator it = memberRules.begin(); it != memberRules.end(); ++it)
    {
        argRules.push_back( (*it)->clone() );
    }
    
}


/** Constructor */
DistributionFunctionCdf::DistributionFunctionCdf(const DistributionFunctionCdf& obj) : TypedFunction<Probability>(obj),
    argRules( obj.argRules )
{
    
    if ( obj.templateObject != NULL )
    {
        templateObject = obj.templateObject->clone();
    }
    else
    {
        templateObject = NULL;
    }
    
    if ( obj.templateObjectPositive != NULL )
    {
        templateObjectPositive = obj.templateObjectPositive->clone();
    }
    else
    {
        templateObjectPositive = NULL;
    }
    
}


/** Constructor */
DistributionFunctionCdf::~DistributionFunctionCdf(void)
{
    
    delete templateObject;
    delete templateObjectPositive;
    
}


DistributionFunctionCdf& DistributionFunctionCdf::operator=(const DistributionFunctionCdf &c) {
    
    if (this != &c)
    {
        Function::operator=(c);
        delete templateObject;
        delete templateObjectPositive;
        
        if ( c.templateObject != NULL )
        {
            templateObject = c.templateObject->clone();
        }
        else
        {
            templateObject = NULL;
        }
        
        if ( c.templateObjectPositive != NULL )
        {
            templateObjectPositive = c.templateObjectPositive->clone();
        }
        else
        {
            templateObjectPositive = NULL;
        }
        
        argRules = c.argRules;
    }
    
    return *this;
}


/** Clone the object */
DistributionFunctionCdf* DistributionFunctionCdf::clone(void) const
{
    
    return new DistributionFunctionCdf(*this);
}


RevBayesCore::TypedFunction<double>* DistributionFunctionCdf::createFunction( void ) const
{
    
    RevBayesCore::ContinuousDistribution *d = NULL;
    
    if ( templateObject != NULL )
    {
    
        ContinuousDistribution* copyObject = templateObject->clone();
    
        for ( size_t i = 1; i < args.size(); i++ )
        {
        
            if ( args[i].isConstant() )
            {
                copyObject->setConstParameter( args[i].getLabel(), RevPtr<const RevVariable>( (const RevVariable*) args[i].getVariable() ) );
            }
            else
            {
                copyObject->setParameter( args[i].getLabel(), args[i].getReferenceVariable() );
            }
        }
        
        d = copyObject->createDistribution();
        
    }
    else
    {
        
        
        PositiveContinuousDistribution* copyObject = templateObjectPositive->clone();
        
        for ( size_t i = 1; i < args.size(); i++ )
        {
            
            if ( args[i].isConstant() )
            {
                copyObject->setConstParameter( args[i].getLabel(), RevPtr<const RevVariable>( (const RevVariable*) args[i].getVariable() ) );
            }
            else
            {
                copyObject->setParameter( args[i].getLabel(), args[i].getReferenceVariable() );
            }
        }
        
        d = copyObject->createDistribution();
        
    }
    
    
    
    RevBayesCore::TypedDagNode<double>* arg = static_cast<const Probability &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::CumulativeDistributionFunction* f = new RevBayesCore::CumulativeDistributionFunction( arg, d );
    
    return f;
}


/** Get argument rules */
const ArgumentRules& DistributionFunctionCdf::getArgumentRules(void) const
{
    
    return argRules;
}


/** Get Rev type of object */
const std::string& DistributionFunctionCdf::getClassType(void)
{
    
    static std::string rev_type = "DistributionFunctionCdf";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& DistributionFunctionCdf::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the aliases for the function.
 * We simple return the aliases of the distribution.
 */
std::vector<std::string> DistributionFunctionCdf::getFunctionNameAliases( void ) const
{
    
    std::vector<std::string> dist_aliases;
    if ( templateObject != NULL )
    {
        dist_aliases = templateObject->getDistributionFunctionAliases();
    }
    else if ( templateObjectPositive != NULL )
    {
        dist_aliases = templateObjectPositive->getDistributionFunctionAliases();
    }
    std::vector<std::string> aliases;

    for (size_t i = 0; i < dist_aliases.size(); ++i)
    {
        std::string f_name = "p" + dist_aliases[i];
        aliases.push_back( f_name );
    }
    
    return aliases;
}


/**
 * Get the primary Rev name for this function.
 */
std::string DistributionFunctionCdf::getFunctionName( void ) const
{
    // create a name variable that is NOT the same for all instance of this class
    std::string f_name = "p" + (templateObject != NULL ? templateObject->getDistributionFunctionName() : templateObjectPositive->getDistributionFunctionName() );
    
    return f_name;
}


/** Get type spec */
const TypeSpec& DistributionFunctionCdf::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

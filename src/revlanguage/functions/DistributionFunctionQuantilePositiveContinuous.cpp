#include "DistributionFunctionQuantilePositiveContinuous.h"

#include <cstddef>

#include "ArgumentRule.h"
#include "DeterministicNode.h"
#include "Probability.h"
#include "QuantileFunction.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlPositiveContinuousDistribution.h"

namespace RevBayesCore { class ContinuousDistribution; }

using namespace RevLanguage;


/** Constructor */
DistributionFunctionQuantilePositiveContinuous::DistributionFunctionQuantilePositiveContinuous( PositiveContinuousDistribution *d ) : TypedFunction<RealPos>(),
    templateObjectPositive( d )
{
    
    argRules.push_back( new ArgumentRule("p", Probability::getClassTypeSpec(), "The probability for this quantile.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    const ArgumentRules &memberRules = templateObjectPositive->getParameterRules();
    for (std::vector<ArgumentRule*>::const_iterator it = memberRules.begin(); it != memberRules.end(); ++it)
    {
        argRules.push_back( (*it)->clone() );
    }
    
}


/** Constructor */
DistributionFunctionQuantilePositiveContinuous::DistributionFunctionQuantilePositiveContinuous(const DistributionFunctionQuantilePositiveContinuous& obj) : TypedFunction<RealPos>(obj),
    argRules( obj.argRules )
{
    
    if ( obj.templateObjectPositive != NULL )
    {
        templateObjectPositive = obj.templateObjectPositive->clone();
    }
    else
    {
        templateObjectPositive = NULL;
    }
    
}



DistributionFunctionQuantilePositiveContinuous::~DistributionFunctionQuantilePositiveContinuous( void )
{
    delete templateObjectPositive;
}

DistributionFunctionQuantilePositiveContinuous& DistributionFunctionQuantilePositiveContinuous::operator=(const DistributionFunctionQuantilePositiveContinuous &c)
{
    
    if (this != &c)
    {
        TypedFunction<RealPos>::operator=(c);
        
        delete templateObjectPositive;
        
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
DistributionFunctionQuantilePositiveContinuous* DistributionFunctionQuantilePositiveContinuous::clone(void) const {
    
    return new DistributionFunctionQuantilePositiveContinuous(*this);
}


RevBayesCore::TypedFunction<double>* DistributionFunctionQuantilePositiveContinuous::createFunction( void ) const
{
    
    RevBayesCore::TypedFunction<double>* f = NULL;
    RevBayesCore::TypedDagNode<double>* arg = static_cast<const Probability &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
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
        
    RevBayesCore::ContinuousDistribution *d = copyObject->createDistribution();
    f = new RevBayesCore::QuantileFunction( arg, d );
    
    
    // return the value
    return f;
}


/** Get argument rules */
const ArgumentRules& DistributionFunctionQuantilePositiveContinuous::getArgumentRules(void) const {
    
    return argRules;
}


/** Get Rev type of object */
const std::string& DistributionFunctionQuantilePositiveContinuous::getClassType(void)
{
    
    static std::string rev_type = "DistributionFunctionQuantilePositiveContinuous";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& DistributionFunctionQuantilePositiveContinuous::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string DistributionFunctionQuantilePositiveContinuous::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "q" + templateObjectPositive->getDistributionFunctionName();
    
    return f_name;
}


/**
 * Get the aliases for the function.
 * We simple return the aliases of the distribution.
 */
std::vector<std::string> DistributionFunctionQuantilePositiveContinuous::getFunctionNameAliases( void ) const
{
    
    std::vector<std::string> dist_aliases = ( templateObjectPositive != NULL ? templateObjectPositive->getDistributionFunctionAliases() : std::vector<std::string>() );
    std::vector<std::string> aliases;
    
    for (size_t i = 0; i < dist_aliases.size(); ++i)
    {
        std::string f_name = "q" + dist_aliases[i];
        aliases.push_back( f_name );
    }
    
    return aliases;
}


/** Get type spec */
const TypeSpec& DistributionFunctionQuantilePositiveContinuous::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

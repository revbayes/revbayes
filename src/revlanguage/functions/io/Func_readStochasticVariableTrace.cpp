#include "Func_readStochasticVariableTrace.h"

#include <sstream>
#include <vector>

#include "ArgumentRule.h"
#include "Delimiter.h"
#include "RlString.h"
#include "RlModelTrace.h"
#include "TraceReader.h"
#include "WorkspaceVector.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "Probability.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "Trace.h"
#include "TypeSpec.h"
#include "WorkspaceToCoreWrapperObject.h"


using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_readStochasticVariableTrace* Func_readStochasticVariableTrace::clone( void ) const
{
    
    return new Func_readStochasticVariableTrace( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_readStochasticVariableTrace::execute( void )
{
    
    // get the information from the arguments for reading the file
    const std::string&     fn       = static_cast<const RlString&>( args[0].getVariable()->getRevObject() ).getValue();
    // get the column delimiter
    const std::string& delimiter    = static_cast<const RlString&>( args[1].getVariable()->getRevObject() ).getValue();
    
    RevObject& b = args[2].getVariable()->getRevObject();
    
    
    RevBayesCore::TraceReader reader;
    std::vector<RevBayesCore::ModelTrace> data = reader.readStochasticVariableTrace(fn, delimiter);
    
    WorkspaceVector<ModelTrace> *rv = new WorkspaceVector<ModelTrace>();
    for (std::vector<RevBayesCore::ModelTrace>::iterator it = data.begin(); it != data.end(); ++it)
    {
        int burnin = 0;

        if ( b.isType( Integer::getClassTypeSpec() ) )
        {
            burnin = (int)static_cast<const Integer &>(b).getValue();
        }
        else
        {
            double burninFrac = static_cast<const Probability &>(b).getValue();
            burnin = int( floor( it->size()*burninFrac ) );
        }
        
        it->setBurnin(burnin);

        rv->push_back( ModelTrace( *it ) );
    }
    
    // return the vector of traces
    return new RevVariable( rv );
}


/** Get argument rules */
const ArgumentRules& Func_readStochasticVariableTrace::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        
        argumentRules.push_back( new ArgumentRule( "file"     , RlString::getClassTypeSpec(), "The name of the file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new Delimiter() );
        std::vector<TypeSpec> burninTypes;
        burninTypes.push_back( Probability::getClassTypeSpec() );
        burninTypes.push_back( Integer::getClassTypeSpec() );
        argumentRules.push_back( new ArgumentRule( "burnin"   , burninTypes     , "The fraction/number of samples to discard as burnin.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability(0.25) ) );
        // Sebastian: currently thinning is not supported but maybe should be later
//        argumentRules.push_back( new ArgumentRule( "thinning", Natural::getClassTypeSpec(), "The frequency of samples to read, i.e., we will only used every n-th sample where n is defined by this argument.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural( 1 ) ) );

        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_readStochasticVariableTrace::getClassType(void)
{
    
    static std::string rev_type = "Func_readStochasticVariableTrace";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Func_readStochasticVariableTrace::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readStochasticVariableTrace::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readStochasticVariableTrace";
    
    return f_name;
}

/** Get type spec */
const TypeSpec& Func_readStochasticVariableTrace::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_readStochasticVariableTrace::getReturnType( void ) const
{
    static TypeSpec return_typeSpec = WorkspaceVector<ModelTrace>::getClassTypeSpec();
    return return_typeSpec;
}



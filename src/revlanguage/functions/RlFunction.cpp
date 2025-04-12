#include <cstdlib>
#include <sstream>
#include <cstddef>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "OptionRule.h"
#include "RbException.h"
#include "RevObject.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "RbHelpArgument.h"
#include "RbHelpFunction.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"

namespace RevBayesCore { class RbHelpEntry; }
namespace RevLanguage { class Environment; }

using namespace RevLanguage;

/** Basic constructor. */
Function::Function(void) : RevObject( false ),
    args( )
{

    argsProcessed = false;
}

/** Copy constructor. */
Function::Function(const Function &x) : RevObject( x ), 
    argsProcessed( x.argsProcessed ),
    args( x.args ),
    env( x.env )
{
    
}


/** Destructor. We need to free the arguments here. */
Function::~Function(void)
{
    
    // we don't own the enclosing environment -> we don't delete it.
}


/** Debug info about object */
std::string Function::callSignature(void) const
{
    
    std::ostringstream o;
    o << getType() << ": " << std::endl;
    
    if (argsProcessed)
    {

        for ( size_t i = 0;  i < args.size(); i++ )
        {
            o << " args[" << i << "] = ";
            args[i].getVariable()->getRevObject().printValue( o, true);
            o << std::endl;
        }
    }
    else
    {
        
        for ( size_t i = 0;  i < getArgumentRules().size(); i++ )
        {
            o << " args[" << i << "] = ";
            getArgumentRules()[i].printValue( o );
            o << std::endl;
        }
    }
    
    return o.str();
}


/**
 * @brief Check arguments
 *
 * This function processes arguments based on argument rules. First it deletes
 * any previously stored arguments. If the matching of the new arguments
 * succeeds, the processedArguments will be set to the new vector of processed
 * arguments and the function returns true. Any subsequent calls to execute()
 * will then use the processed arguments. You can also call the function with
 * the arguments directly, in which case processArguments will be called first
 * before the operation is actually performed.
 *
 * In matching arguments to argument rules, we use the same rules as in R with
 * the addition that types are also used in the matching process, after arguments
 * have been reordered as in R. The FunctionTable ensure that all argument rules
 * are distinct. However, several functions can nevertheless match the same
 * arguments because of the inheritance hierarchy. In these cases, the closest
 * match is chosen based on the first argument, then on the second, etc.
 *
 * The function returns a match score based on how closely the arguments match the
 * rules, with 0 being perfect match, 1 being a match to an immediate base class type,
 * 2 a match to a grandparent class, etc. A large number is used for arguments that
 * need type conversion.
 *
 * The evaluateOnce parameter MemberObject whether the function is to be evaluated once,
 * immediately after the call to processArguments, or whether the processed arguments
 * will be used in repeated function calls in a function node. Argument matching is
 * based on current values in the first case, but on the wrapper type in the second.
 *
 * @todo Labels need to be stored for ellipsis arguments.
 *
 * These are the argument matching rules:
 *
 *  1. If the last argument rule is an ellipsis, and it is the kth argument passed
 *     in, then all arguments passed in, from position k to the end, are wrapped
 *     in a single ContainerNode object. These arguments are not matched to any
 *     rules.
 *  2. The remaining arguments are matched to labels using exact matching. If the
 *     type does not match the type of the rule, it is an error.
 *  3. The remaining arguments are matched to any remaining slots using partial
 *     matching. If there is ambiguity or the types do not match, it is an error.
 *  4. The remaining arguments are used for the empty slots in the order they were
 *     passed in. If the types do not match, it is an error.
 *  5. Any remaining empty slots are filled with default values stored in the argument
 *     rules (we use copies of the values, of course).
 *  6. If there are still empty slots, the arguments do not match the rules.
 *
 * @todo Fredrik: The code and the logic has been changed without changing the comments, so these
 *       are out of date. Also note that the argument matching is problematic for unlabeled
 *       arguments (order can be changed based on argument types, which may cause unintended
 *       consequences). Furthermore, there is redundant code left from the old implementation.
 *       Finally, the ellipsis arguments no longer have to be last among the rules, but they
 *       are still the last arguments after processing.
 */
bool Function::checkArguments( const std::vector<Argument>& passed_args, std::vector<double>* match_score, std::vector<bool>& arg_mapped, bool once)
{
    
    /*********************  0. Initialization  **********************/
    
    double MAX_ARGS = 1000.0;
    
    /* Get the argument rules */
    const ArgumentRules& the_rules = getArgumentRules();
    
    /* Get the number of argument rules */
    size_t num_rules = the_rules.size();
    
    /* Keep track of which arguments we have used, and which argument slots we have filled, and with what passed arguments */
    std::vector<bool>   taken           = std::vector<bool>( passed_args.size(), false );
    std::vector<bool>   filled          = std::vector<bool>( the_rules.size(), false );
    
    
    /*********************  1. Do exact matching  **********************/
    
    /* Do exact matching of labels */
    for (size_t i=0; i<passed_args.size(); i++)
    {
        
        /* Test if swallowed by ellipsis; if so, we can quit because the remaining passed_args will also be swallowed */
        if ( taken[i] )
        {
            break;
        }
        
        /* Skip if no label */
        if ( passed_args[i].getLabel().size() == 0 )
        {
            continue;
        }
        
        /* Check for matches in all regular rules (we assume that all labels are unique; this is checked by FunctionTable) */
        for (size_t j=0; j<num_rules; j++)
        {
            std::vector<std::string> aliases = the_rules[j].getArgumentAliases();

            for (size_t k=0; k < aliases.size(); k++)
            {
                if ( passed_args[i].getLabel() == aliases[k] )
                {
                    
                    if ( filled[j] )
                    {
                        return false;
                    }
                    
                    Argument &arg = const_cast<Argument&>(passed_args[i]);
                    double penalty = the_rules[j].isArgumentValid( arg );
                    if ( penalty != -1 )
                    {
                        taken[i]          = true;
                        filled[j]         = true;
                        arg_mapped[i]     = true;

                        if ( match_score != NULL)
                        {
                            double score = computeMatchScore(passed_args[i].getVariable(), the_rules[j]);
                            score += abs(int(i)-int(j)) / MAX_ARGS;
                            score += penalty*100.0;
                            match_score->push_back(score);
                        }
                    }
                    else
                    {
                        return false;
                    }

                    // We got an exact match -> we can skip the other labels for checking
                    j = num_rules;
                    break;
                }
            }
        }
    }
    
    
    /*********************  2. Do partial matching  **********************/
    
    /* Do partial matching of labels */
    for (size_t i=0; i<passed_args.size(); i++)
    {
        
        /* Skip if already matched */
        if ( taken[i] )
        {
            continue;
        }
        
        /* Skip if no label */
        if ( passed_args[i].getLabel().size() == 0 )
        {
            continue;
        }
        
        /* Initialize match index and number of matches */
        int num_matches = 0;
        int match_rule = -1;
        
        /* Try all rules */
        for (size_t j=0; j<num_rules; j++)
        {
            std::vector<std::string> aliases = the_rules[j].getArgumentAliases();

            for (size_t k=0; k < aliases.size(); k++)
            {
                if ( filled[j] == false && aliases[k].compare(0, passed_args[i].getLabel().size(), passed_args[i].getLabel()) == 0 )
                {
                    ++num_matches;
                    match_rule = static_cast<int>( j );
                }
            }
        }
        
        if ( num_matches == 1 )
        {
            
            Argument &arg = const_cast<Argument&>(passed_args[i]);
            double penalty = the_rules[match_rule].isArgumentValid(arg );
            if ( penalty != -1 )
            {
                taken[i]            = true;
                filled[match_rule]  = true;
                arg_mapped[i]       = true;
            
                if ( match_score != NULL)
                {
                    double score = computeMatchScore(passed_args[i].getVariable(), the_rules[match_rule]);
                    score += abs(int(i)-int(match_rule)) / MAX_ARGS;
                    score += penalty*100.0;
                    match_score->push_back(score);
                }
            }
            else
            {
                return false;
            }
        }
    
    }
    
    
    /*********************  3. Fill with unused passed_args  **********************/
    
    /* Fill in empty slots using the remaining args in order */
    for (size_t i=0; i<passed_args.size(); i++)
    {
        
        /* Skip if already matched */
        if ( taken[i] )
        {
            continue;
        }
        
        /* Find first empty slot and try to fit argument there */
        for (size_t j=0; j<num_rules; j++)
        {
            
            if ( filled[j] == false )
            {
                
                Argument &arg = const_cast<Argument&>(passed_args[i]);
                double penalty = the_rules[j].isArgumentValid( arg );
                if ( penalty != -1 )
                {
                    taken[i]          = true;
                    arg_mapped[i]     = true;
                    if ( the_rules[j].isEllipsis() == false )
                    {
                        filled[j]     = true;
                    }
                    
                    if ( match_score != NULL)
                    {
                        const RevPtr<const RevVariable>& arg_var = passed_args[i].getVariable();
                        double score = computeMatchScore(arg_var, the_rules[j]);
                        score += abs(int(i)-int(j)) / MAX_ARGS;
                        score += penalty*100.0;
                        match_score->push_back(score);
                    }
                    
                    break;
                }
            }
        }
        
        /* Final test if we found a match */
        if ( taken[i] == false )
        {
            return false;
        }
        
    }
    
    
    
    /*********************  4. Fill with default values  **********************/
    
    /* Fill in empty slots using default values */
    for (size_t i=0; i<num_rules; i++)
    {
        
        if ( filled[i] == true || the_rules[i].isEllipsis() )
        {
            continue;
        }
        
        if ( the_rules[i].hasDefault() == false )
        {
            return false;
        }
        
    }
    
    return true;
}


/* Delete processed args */
void Function::clear(void)
{

    args.clear();
    
    argsProcessed = false;
    
    clearArguments();
}


/* Delete processed args */
void Function::clearArguments(void)
{
    
    // We have already done everything, but derived classes might want to overwrite this function.
    
}


/** Compute the match score between the argument and the argument rule. */
double Function::computeMatchScore(const RevVariable *var, const ArgumentRule &rule)
{
   
    double     score = 10000;   // Needs to be larger than the max depth of the class hierarchy

    const TypeSpec& arg_class = var->getRevObject().getTypeSpec();
    const std::vector<TypeSpec> &rule_arg_types = rule.getArgumentTypeSpec();
    for ( std::vector<TypeSpec>::const_iterator it = rule_arg_types.begin(); it != rule_arg_types.end(); ++it)
    {
        int j = 0;
        const TypeSpec* parent = &arg_class;
        do
        {
            if ( *parent == *it )
            {
                score = j;
                break;
            }
            parent = parent->getParentTypeSpec();
            j++;
            if ( j >= score ) 
            {
                break;
            }
        } while (parent != NULL);
    }
    
    return score;    // We needed type conversion for this argument
}


const std::vector<Argument>& Function::getArguments(void) const
{
    return args;
}


std::vector<Argument>& Function::getArguments(void)
{
    return args;
}


/** Get Rev type of object */
const std::string& Function::getClassType(void)
{
    
    static std::string rev_type = "Function";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Function::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/** Get execution environment of function */
std::shared_ptr<Environment> Function::getEnvironment(void) const
{
    
    return env;
}


RevBayesCore::RbHelpFunction* Function::constructTypeSpecificHelp( void ) const
{
    
    return new RevBayesCore::RbHelpFunction();
}


/**
 * Get the help entry for this class
 */
void Function::addSpecificHelpFields(RevBayesCore::RbHelpEntry *e) const
{
    // create the help function entry that we will fill with some values
    RevBayesCore::RbHelpFunction *help = static_cast<RevBayesCore::RbHelpFunction*>( e );
    RevBayesCore::RbHelpFunction &help_entry = *help;
    
    // usage
    help_entry.setUsage( getHelpUsage() );
    
    // arguments
    const MemberRules& rules = getArgumentRules();
    std::vector<RevBayesCore::RbHelpArgument> arguments = std::vector<RevBayesCore::RbHelpArgument>();
    
    for ( size_t i=0; i<rules.size(); ++i )
    {
        const ArgumentRule &the_rule = rules[i];
        
        RevBayesCore::RbHelpArgument argument = RevBayesCore::RbHelpArgument();
        
        std::string label;
        const std::vector<std::string> aliases = the_rule.getArgumentAliases();
        for ( std::vector<std::string>::const_iterator it = aliases.begin(); it != aliases.end(); it++ )
        {
            if ( it != aliases.begin() )
            {
                label += "/";
            }
            label += *it;
        }
        argument.setLabel( label );
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
    
    help_entry.setArguments( arguments );
    
    // return value
    help_entry.setReturnType( getReturnType().getType() );
    
}


/**
 * Get the name for this procedure.
 */
std::string Function::getConstructorFunctionName( void ) const
{
    return getFunctionName();
}


/**
 * Get the name for this procedure.
 */
std::vector<std::string> Function::getConstructorFunctionAliases( void ) const
{
    
    return getFunctionNameAliases();
}


/**
 * Get the usage for this function.
 */
std::string Function::getHelpUsage( void ) const
{
    std::string usage = getRevDeclaration();
    
    return usage;
}


/**
 * Get the name for this procedure.
 */
std::vector<std::string> Function::getFunctionNameAliases( void ) const
{
    std::vector<std::string> aliases;
    
    return aliases;
}


/** 
 * Get Rev declaration of the function, formatted for output
 */
std::string Function::getRevDeclaration(void) const
{
    
    std::ostringstream o;
 
    // Sebastian: We don't want to print the return type in the usage.
    // It only confuses.
//    o << getReturnType().getType();
    if ( getFunctionName() == "" )
    {
        o << "<unnamed>(";
    }
    else
    {
        o << "" << getFunctionName() << "(";
    }
    
    const ArgumentRules& argRules = getArgumentRules();
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


/** Print value for user */
void Function::printValue(std::ostream& o, bool user) const
{
    const ArgumentRules& argRules = getArgumentRules();

    size_t n_space = getFunctionName().size() + 2;
    std::string space(n_space, ' ');
    
    o << getFunctionName() << " (";
    for (size_t i=0; i<argRules.size(); i++)
    {
        if (i != 0)
        {
            o << ",\n";
            o << space;
        }
        argRules[i].printValue(o);
    }
    o << ")";
}


/**
 * @brief Process arguments
 *
 * This function processes arguments based on argument rules. First it deletes
 * any previously stored arguments. If the matching of the new arguments
 * succeeds, the processedArguments will be set to the new vector of processed
 * arguments and the function returns true. Any subsequent calls to execute()
 * will then use the processed arguments. You can also call the function with
 * the arguments directly, in which case processArguments will be called first
 * before the operation is actually performed.
 *
 * In matching arguments to argument rules, we use the same rules as in R with
 * the addition that types are also used in the matching process, after arguments
 * have been reordered as in R. The FunctionTable ensure that all argument rules
 * are distinct. However, several functions can nevertheless match the same
 * arguments because of the inheritance hierarchy. In these clases, the closest
 * match is chosen based on the first argument, then on the second, etc.
 *
 * The function returns a match score based on how closely the arguments match the
 * rules, with 0 being perfect match, 1 being a match to an immediate base class type,
 * 2 a match to a grandparent class, etc. A large number is used for arguments that
 * need type conversion.
 *
 * The evaluateOnce parameter MemberObject whether the function is to be evaluated once,
 * immediately after the call to processArguments, or whether the processed arguments
 * will be used in repeated function calls in a function node. Argument matching is
 * based on current values in the first case, but on the wrapper type in the second.
 *
 * @todo Labels need to be stored for ellipsis arguments.
 *
 * These are the argument matching rules:
 *
 *  1. If the last argument rule is an ellipsis, and it is the kth argument passed
 *     in, then all arguments passed in, from position k to the end, are wrapped
 *     in a single ContainerNode object. These arguments are not matched to any
 *     rules.
 *  2. The remaining arguments are matched to labels using exact matching. If the
 *     type does not match the type of the rule, it is an error.
 *  3. The remaining arguments are matched to any remaining slots using partial
 *     matching. If there is ambiguity or the types do not match, it is an error.
 *  4. The remaining arguments are used for the empty slots in the order they were
 *     passed in. If the types do not match, it is an error.
 *  5. Any remaining empty slots are filled with default values stored in the argument
 *     rules (we use copies of the values, of course).
 *  6. If there are still empty slots, the arguments do not match the rules.
 *
 * @todo Fredrik: The code and the logic has been changed without changing the comments, so these
 *       are out of date. Also note that the argument matching is problematic for unlabeled
 *       arguments (order can be changed based on argument types, which may cause unintended
 *       consequences). Furthermore, there is redundant code left from the old implementation.
 *       Finally, the ellipsis arguments no longer have to be last among the rules, but they
 *       are still the last arguments after processing.
 *
 * @todo Fredrik: Static and dynamic type conversion added, but partly hack-ish, so the implementation
 *       needs to be revised
 */
void Function::processArguments( const std::vector<Argument>& passed_args )
{

    /*********************  0. Initialization  **********************/
    /* Get my own copy of the argument vector */
    std::vector<Argument> p_args = passed_args;
    
    /* Get the argument rules */
    const ArgumentRules& the_rules = getArgumentRules();

    /* Get the number of argument rules */
    size_t num_rules = the_rules.size();

    /* Clear previously processed arguments */
    args.clear();
    
    /* Keep track of which arguments we have used, and which argument slots we have filled, and with what passed arguments */
    std::vector<bool>   taken            = std::vector<bool>( passed_args.size(), false );
    std::vector<bool>   filled           = std::vector<bool>( num_rules, false );
    std::vector<size_t> passed_arg_index = std::vector<size_t>( num_rules, 1000 );
    std::vector<Argument> ellipsis_args;

    /*********************  1. Do exact matching  **********************/

    /* Do exact matching of labels */
    for (size_t i=0; i<passed_args.size(); i++)
    {

        /* Test if swallowed by ellipsis; if so, we can quit because the remaining passed_args will also be swallowed */
        if ( taken[i] )
        {
            break;
        }
        
        /* Skip if no label */
        if ( passed_args[i].getLabel().size() == 0 )
        {
            continue;
        }
        
        /* Check for matches in all regular rules (we assume that all labels are unique; this is checked by FunctionTable) */
        for (size_t j=0; j<num_rules; j++)
        {
            std::vector<std::string> aliases = the_rules[j].getArgumentAliases();

            for (size_t k=0; k < aliases.size(); k++)
            {
                if ( passed_args[i].getLabel() == aliases[k] )
                {

                    if ( filled[j] )
                    {
                        throw RbException() << "Duplicate argument labels '" << passed_args[i].getLabel() ; 
                    }

                    p_args[i]               = the_rules[j].fitArgument( p_args[i] );
                    taken[i]                = true;
                    filled[j]               = true;
                    passed_arg_index[j]     = static_cast<int>( i );

                    // We got an exact match -> we can skip the other labels for checking
                    j = num_rules;
                    break;
                }
            }
        }
    }

 
    /*********************  2. Do partial matching  **********************/

    /* Do partial matching of labels */
    for (size_t i=0; i<passed_args.size(); i++)
    {

        /* Skip if already matched */
        if ( taken[i] )
        {
            continue;
        }
        
        /* Skip if no label */
        if ( passed_args[i].getLabel().size() == 0 )
        {
            continue;
        }
        
        /* Initialize match index and number of matches */
        int num_matches = 0;
        int match_rule = -1;
        std::string label;

        /* Try all rules */
        for (size_t j=0; j<num_rules; j++)
        {
            std::vector<std::string> aliases = the_rules[j].getArgumentAliases();

            for (size_t k=0; k < aliases.size(); k++)
            {
                if ( !filled[j] && aliases[k].compare(0, passed_args[i].getLabel().size(), passed_args[i].getLabel()) == 0 )
                {
                    ++num_matches;
                    match_rule = static_cast<int>( j );
                    label = aliases[k];
                }
            }
        }

        if (num_matches > 1)
        {
            throw RbException() << "Argument label '" << passed_args[i].getLabel() << "' matches mutliple parameter labels." ; 
        }
        else if (num_matches < 1)
        {
            throw RbException() << "Argument label '" << passed_args[i].getLabel() << "' matches no untaken parameter labels." ; 
        }
        
        if ( num_matches == 1)
        {
            p_args[i].setLabel(label);
            p_args[i]                       = the_rules[match_rule].fitArgument( p_args[i] );
            taken[i]                        = true;
            filled[match_rule]              = true;
            passed_arg_index[match_rule]    = static_cast<int>( i );
        }
        
    }


    /*********************  3. Fill with unused passed_args  **********************/

    /* Fill in empty slots using the remaining args in order */
    for (size_t i=0; i<passed_args.size(); i++)
    {

        /* Skip if already matched */
        if ( taken[i] )
        {
            continue;
        }
        
        /* Find first empty slot and try to fit argument there */
        for (size_t j=0; j<num_rules; j++)
        {

            if ( filled[j] == false &&
                 ( (!the_rules[j].isEllipsis() && passed_args[i].getLabel().size() == 0) || (the_rules[j].isEllipsis()) ) )
            {
                
                Argument &arg = const_cast<Argument&>(passed_args[i]);
                double penalty = the_rules[j].isArgumentValid( arg );
                if ( penalty != -1 )
                {
                    p_args[i].setLabel( the_rules[j].getArgumentAliases().front() );
                    p_args[i]           = the_rules[j].fitArgument( p_args[i] );
                    taken[i]            = true;
                    if ( the_rules[j].isEllipsis() == false )
                    {
                        filled[j]       = true;
                        passed_arg_index[j] = static_cast<int>( i );
                    }
                    else 
                    {
                        ellipsis_args.push_back( p_args[i] );
                    }
                    
                    break;
                }
            }
        }
        
        /* Final test if we found a match */
        if ( taken[i] == false )
        {
            throw RbException() << "Superfluous argument of type '" << passed_args[i].getVariable()->getRevObject().getType() << "' and name '" << passed_args[i].getLabel() << "' passed to function '" << getType() << "'.";
        }
    }

    /*********************  4. Fill with default values  **********************/

    /* Fill in empty slots using default values */
    for (size_t i=0; i<num_rules; i++)
    {

        if ( filled[i] == true || the_rules[i].isEllipsis())
        {
            continue;
        }
        
        // we just leave the optional arguments empty
        if ( !the_rules[i].hasDefault() )
        {
            throw RbException() << "No argument found for parameter '" << the_rules[i].getArgumentLabel() << "'.";
        }
        
        const ArgumentRule& theRule = the_rules[i];
        RevPtr<RevVariable> the_var = theRule.getDefaultVariable().clone();
        the_var->setName( "." + theRule.getArgumentLabel() );
        the_var->setRequiredTypeSpec( theRule.getDefaultVariable().getRequiredTypeSpec() );
        size_t idx = p_args.size();
        passed_arg_index[i] = idx;
        p_args.push_back( Argument( the_var, theRule.getArgumentLabel(), theRule.getEvaluationType() != ArgumentRule::BY_CONSTANT_REFERENCE ) );
    }

    argsProcessed = true;
    
    /*********************  5. Insert arguments into argument list  **********************/
    for (size_t j=0; j<num_rules; j++)
    {
        if ( passed_arg_index[j] < 1000 )
        {
            args.push_back( p_args[ passed_arg_index[j] ] );
        }
        
    }
    
    /*********************  6. Insert ellipsis arguments  **********************/
    for (std::vector<Argument>::iterator i = ellipsis_args.begin(); i != ellipsis_args.end(); i++) 
    {
        args.push_back( *i );
    }

}


void Function::setExecutionEnviroment(const std::shared_ptr<Environment>& e)
{
    
    env = e;

}





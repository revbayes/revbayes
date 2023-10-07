#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "RbException.h"
#include "RlFunction.h"
#include "TypeSpec.h"
#include "Workspace.h"
#include "DagNode.h"
#include "Environment.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"

using namespace RevLanguage;


/**
 * Construct rule with single type;
 * use "" for no label.
 */
ArgumentRule::ArgumentRule(const std::string& argName, const TypeSpec& argTypeSp, const std::string& argDesc, EvaluationType et, DagNodeType dt, RevObject *defVal) :
    argTypeSpecs( 1, argTypeSp ),
    defaultVar( new RevVariable( defVal ) ),
    evalType( et ),
    nodeType( dt ),
    aliases( std::vector<std::string>(1, argName) ),
    label( argName ),
    description( argDesc ),
    hasDefaultVal( true )
{
    
}


/**
 * Construct rule with multiple types;
 * use "" for no label.
 */
ArgumentRule::ArgumentRule(const std::string& argName, const std::vector<TypeSpec>& argTypeSp, const std::string& argDesc, EvaluationType et, DagNodeType dt, RevObject *defVal) :
    argTypeSpecs( argTypeSp ),
    defaultVar( new RevVariable( defVal ) ),
    evalType( et ),
    nodeType( dt ),
    aliases( std::vector<std::string>(1, argName) ),
    label( argName ),
    description( argDesc ),
    hasDefaultVal( true )
{
    
}


/**
 * Construct rule with single type;
 * use "" for no label.
 */
ArgumentRule::ArgumentRule(const std::string& argName, const TypeSpec& argTypeSp, const std::string& argDesc, EvaluationType et, DagNodeType dt) :
    argTypeSpecs( 1, argTypeSp ),
    defaultVar( NULL ),
    evalType( et ),
    nodeType( dt ),
    aliases( std::vector<std::string>(1, argName) ),
    label( argName ),
    description( argDesc ),
    hasDefaultVal( false )
{

}


/**
 * Construct rule with multiple types;
 * use "" for no label.
 */
ArgumentRule::ArgumentRule(const std::string& argName, const std::vector<TypeSpec>& argTypeSp, const std::string& argDesc, EvaluationType et, DagNodeType dt) :
    argTypeSpecs( argTypeSp ),
    defaultVar( NULL ),
    evalType( et ),
    nodeType( dt ),
    aliases( std::vector<std::string>(1, argName) ),
    label( argName ),
    description( argDesc ),
    hasDefaultVal( false )
{

}


/**
 * Construct rule with single type and multiple names;
 * use "" for no label.
 */
ArgumentRule::ArgumentRule(const std::vector<std::string>& argNames, const TypeSpec& argTypeSp, const std::string& argDesc, EvaluationType et, DagNodeType dt) :
    argTypeSpecs( 1, argTypeSp ),
    defaultVar( NULL ),
    evalType( et ),
    nodeType( dt ),
    aliases( argNames ),
    label( argNames.front() ),
    description( argDesc ),
    hasDefaultVal( false )
{

}


/**
 * Construct rule with multiple types and multiple names;
 * use "" for no label.
 */
ArgumentRule::ArgumentRule(const std::vector<std::string>& argNames, const std::vector<TypeSpec>& argTypeSp, const std::string& argDesc, EvaluationType et, DagNodeType dt) :
    argTypeSpecs( argTypeSp ),
    defaultVar( NULL ),
    evalType( et ),
    nodeType( dt ),
    aliases( argNames ),
    label( argNames.front() ),
    description( argDesc ),
    hasDefaultVal( false )
{

}


/**
 * Construct rule with single type;
 * use "" for no label.
 */
ArgumentRule::ArgumentRule(const std::vector<std::string>& argNames, const TypeSpec& argTypeSp, const std::string& argDesc, EvaluationType et, DagNodeType dt, RevObject *defVal) :
    argTypeSpecs( 1, argTypeSp ),
    defaultVar( new RevVariable( defVal ) ),
    evalType( et ),
    nodeType( dt ),
    aliases( argNames ),
    label( argNames.front() ),
    description( argDesc ),
    hasDefaultVal( true )
{

}


/**
 * Construct rule with multiple types;
 * use "" for no label.
 */
ArgumentRule::ArgumentRule(const std::vector<std::string>& argNames, const std::vector<TypeSpec>& argTypeSp, const std::string& argDesc, EvaluationType et, DagNodeType dt, RevObject *defVal) :
    argTypeSpecs( argTypeSp ),
    defaultVar( new RevVariable( defVal ) ),
    evalType( et ),
    nodeType( dt ),
    aliases( argNames ),
    label( argNames.front() ),
    description( argDesc ),
    hasDefaultVal( true )
{

}


ArgumentRule* RevLanguage::ArgumentRule::clone( void ) const
{

    return new ArgumentRule( *this );
}


/**
 * Fit a variable into an argument according to the argument rule. If necessary and
 * appropriate, we do type conversion or type promotion.
 *
 *
 *
 * @todo To conform to the old code we change the required type of the incoming
 *       variable wrapper here. We need to change this so that we do not change
 *       the wrapper here, but make sure that if the argument variable is inserted
 *       in a member variable or container element slot, that the slot variable
 *       wrapper, which should be unique (not the same as the incoming variable
 *       wrapper), has the right required type.
 */
Argument ArgumentRule::fitArgument( Argument& arg, bool once ) const
{

    RevPtr<RevVariable> the_var = arg.getVariable();
    if ( evalType == BY_VALUE || the_var->isWorkspaceVariable() || the_var->getRevObject().isConstant() )
    {
        once = true;
    }
    
    
    for ( std::vector<TypeSpec>::const_iterator it = argTypeSpecs.begin(); it != argTypeSpecs.end(); ++it )
    {
        if ( evalType == BY_VALUE || the_var->isWorkspaceVariable() == true )
        {
            if ( the_var->getRevObject().isType( *it ) )
            {
                RevPtr<RevVariable> valueVar = RevPtr<RevVariable>( new RevVariable(the_var->getRevObject().clone(),the_var->getName() ) );
                return Argument( valueVar, arg.getLabel(), false );
            }
            else if ( the_var->getRevObject().isConvertibleTo( *it, once ) != -1 )
            {
                // Fit by type conversion. For now, we also modify the type of the incoming variable wrapper.
                RevObject* convertedObject = the_var->getRevObject().convertTo( *it );

                RevPtr<RevVariable> valueVar = RevPtr<RevVariable>( new RevVariable(convertedObject,the_var->getName() ) );
                return Argument( valueVar, arg.getLabel(), false );
                
            }
        } // if (by-value)
        else
        {
            if ( the_var->getRevObject().isType( *it ) )
            {
                // For now, change the required type of the incoming variable wrapper
                the_var->setRequiredTypeSpec( *it );
            
                if ( isEllipsis() == false )
                {
                    return Argument( the_var, arg.getLabel(), evalType == BY_CONSTANT_REFERENCE );
                }
                else
                {
                    return Argument( the_var, arg.getLabel(), true );
                }
            
            }
            else if ( the_var->getRevObject().isConvertibleTo( *it, once ) != -1  && (*it).isDerivedOf( the_var->getRequiredTypeSpec() ) )
            {
                // Fit by type conversion. For now, we also modify the type of the incoming variable wrapper.
                RevObject* converted_object = the_var->getRevObject().convertTo( *it );
                
                RevPtr<RevVariable> the_new_var = NULL;
                if ( the_var->getRevObject().isConstant() == true )
                {
                    the_new_var = the_var;
                    the_new_var->replaceRevObject( converted_object );
                    the_new_var->setRequiredTypeSpec( *it );
                }
                else
                {
                    the_new_var = RevPtr<RevVariable>( new RevVariable(converted_object, the_var->getName() ) );
                }
                
                if ( !isEllipsis() )
                {
                    return Argument( the_new_var, arg.getLabel(), false );
                }
                else
                {
                    return Argument( the_new_var, arg.getLabel(), false );
                }
            }
            else
            {
                // Fit by type conversion function
            
                const TypeSpec& typeFrom = the_var->getRevObject().getTypeSpec();
                const TypeSpec& typeTo   = *it;
            
                // create the function name
                std::string function_name = "_" + typeFrom.getType() + "2" + typeTo.getType();
                
                // Package arguments
                std::vector<Argument> args;
                Argument theArg = Argument( the_var, "arg" );
                args.push_back( the_var );
                
                Environment& env = Workspace::globalWorkspace();
            
                try
                {
                    Function* func = env.getFunction(function_name, args, once).clone();

                    // Allow the function to process the arguments
                    func->processArguments( args, once );
            
                    // Set the execution environment of the function
                    func->setExecutionEnviroment( &env );
                
                    // Evaluate the function
                    RevPtr<RevVariable> conversionVar = func->execute();
                
                    // free the memory
                    delete func;
                
                    conversionVar->setHiddenVariableState( true );
                    conversionVar->setRequiredTypeSpec( *it );
                
                    return Argument( conversionVar, arg.getLabel(), evalType == BY_CONSTANT_REFERENCE );
                
                }
                catch (RbException& e)
                {
                // we do nothing here
                }
                
            } // else (type conversion function)

            
        } // else (not by-value)

    }
        
    throw RbException( "Argument type mismatch while fitting variable with name \"" + the_var->getName() + "\" of type " + the_var->getRevObject().getType() + " to the argument with name \"" + getArgumentLabel() + "\" and type " +
                        getArgumentTypeSpec()[0].getType()  );
}


const std::vector<std::string>& ArgumentRule::getArgumentAliases( void ) const
{
    return aliases;
}


ArgumentRule::DagNodeType ArgumentRule::getArgumentDagNodeType( void ) const
{
    // return the internal value
    return nodeType;
}


const std::string& ArgumentRule::getArgumentLabel( void ) const
{
    return label;
}


const std::vector<TypeSpec>& ArgumentRule::getArgumentTypeSpec(void) const
{
    return argTypeSpecs;
}



const RevVariable& ArgumentRule::getDefaultVariable( void ) const
{
    
    if ( defaultVar == NULL ) 
    {
        throw RbException("Cannot get default variable \"" + label + "\"");
    }
    
    return *defaultVar;
}


const std::string& ArgumentRule::getArgumentDescription( void ) const
{
    return description;
}


ArgumentRule::EvaluationType ArgumentRule::getEvaluationType( void ) const
{
    
    return evalType;
}


bool ArgumentRule::hasDefault(void) const
{
    
    return hasDefaultVal;
}


/**
 * Test if argument is valid. The boolean flag 'once' is used to signal whether the argument matching
 * is done in a static or a dynamic context. If the rule is constant, then the argument matching
 * is done in a static context (evaluate-once context) regardless of the setting of the once flag.
 * If the argument is constant, we try type promotion if permitted by the variable required type.
 *
 * @todo See the TODOs for fitArgument(...)
 */
double ArgumentRule::isArgumentValid( Argument &arg, bool once) const
{
    
    RevPtr<RevVariable> the_var = arg.getVariable();
    if ( the_var == NULL )
    {
        return -1;
    }
    
    if ( evalType == BY_VALUE || the_var->isWorkspaceVariable() || ( the_var->getRevObject().isModelObject() && the_var->getRevObject().getDagNode()->getDagNodeType() == RevBayesCore::DagNode::CONSTANT) )
    {
        once = true;
    }
    
    if ( nodeType == STOCHASTIC && the_var->getRevObject().getDagNode()->getDagNodeType() != RevBayesCore::DagNode::STOCHASTIC )
    {
        return -1;
    }
    else if ( nodeType == DETERMINISTIC && the_var->getRevObject().getDagNode()->getDagNodeType() != RevBayesCore::DagNode::DETERMINISTIC )
    {
        return -1;
    }

    for ( std::vector<TypeSpec>::const_iterator it = argTypeSpecs.begin(); it != argTypeSpecs.end(); ++it )
    {
        if ( the_var->getRevObject().isType( *it ) )
        {
            return 0.0;
        }
        
        double penalty = -1;
        // make sure that we only perform type casting when the variable will not be part of a model graph
        if ( once == true || the_var->getRevObject().isConstant() == true )
        {
            penalty = the_var->getRevObject().isConvertibleTo( *it, once );
        }
        
        if ( penalty != -1 && (*it).isDerivedOf( the_var->getRequiredTypeSpec() ) )
        {
            return penalty;
        }
        else if ( penalty != -1 && evalType == BY_VALUE )
        {
            return penalty;
        }

//        else if ( once == true &&
////                 !var->isAssignable() &&
//                  the_var->getRevObject().isConvertibleTo( *it, true ) != -1 &&
//                  (*it).isDerivedOf( the_var->getRequiredTypeSpec() )
//                )
//        {
//            return the_var->getRevObject().isConvertibleTo( *it, true );
//        }
        else if ( nodeType != STOCHASTIC )
        {
            
            const TypeSpec& typeFrom = the_var->getRevObject().getTypeSpec();
            const TypeSpec& typeTo   = *it;
            
            // create the function name
            std::string function_name = "_" + typeFrom.getType() + "2" + typeTo.getType();
            
            // Package arguments
            std::vector<Argument> args;
            Argument theArg = Argument( the_var, "arg" );
            args.push_back( the_var );
            
            Environment& env = Workspace::globalWorkspace();
            try
            {
                // we just want to check if the function exists and can be found
                env.getFunction(function_name, args, once);
                return 0.1;
            }
            catch (RbException& e)
            {
                // we do nothing here
            }

        }
        
    }
    
    return -1;
}



bool RevLanguage::ArgumentRule::isEllipsis( void ) const
{
    
    return false;
}
 


/**
 * Print value for user (in descriptions of functions, for instance). We apparently do
 * not use the isConst flag to denote whether an argument is supposed to be passed as
 * a constant currently, so the printing of this modifier is suspended for now.
 *
 */
void ArgumentRule::printValue(std::ostream &o) const
{

    for ( std::vector<TypeSpec>::const_iterator it = argTypeSpecs.begin(); it != argTypeSpecs.end(); ++it )
    {
        if ( it != argTypeSpecs.begin() )
        {
            o << "|";
        }
        
        o << (*it).getType();
    }
    
    // create the default DAG type of the passed-in argument
    std::string dagtype = "";
    // get the type if the variable wasn't NULL
    if ( nodeType == ArgumentRule::DETERMINISTIC )
    {
        dagtype = "<deterministic>";
    }
    else if ( nodeType == ArgumentRule::STOCHASTIC )
    {
        dagtype = "<stochastic>";
    }
    else if ( nodeType == ArgumentRule::CONSTANT )
    {
        dagtype = "<constant>";
    }
    else if ( nodeType == ArgumentRule::DYNAMIC )
    {
        dagtype = "<dynamic>";
    }
    else if ( nodeType == ArgumentRule::ANY )
    {
        dagtype = "<any>";
    }
    else
    {
        dagtype = "<?>";
    }
    o << dagtype;
    
    o << " " << label;
}


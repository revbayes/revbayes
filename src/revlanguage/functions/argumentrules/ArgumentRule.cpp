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
Argument ArgumentRule::fitArgument( Argument& arg ) const
{
    bool convert_by_value = false;
    RevPtr<RevVariable> the_var = arg.getVariable();
    if (not the_var->getRevObject().hasDagNode() or the_var->getRevObject().getDagNode()->isConstant())
    {
        convert_by_value = true;
    }
    else if (evalType == BY_VALUE || the_var->isWorkspaceVariable() || the_var->getRevObject().isConstant())
    {
        convert_by_value = true;
    }
    
    for ( auto& argTypeSpec: argTypeSpecs )
    {
        bool by_value = evalType == BY_VALUE || the_var->isWorkspaceVariable() == true;

        if ( the_var->getRevObject().isType( argTypeSpec ) )
        {
            if (by_value)
            {
                RevPtr<RevVariable> valueVar = RevPtr<RevVariable>( new RevVariable(the_var->getRevObject().clone(),the_var->getName() ) );
                return Argument( valueVar, arg.getLabel(), true );
            }
//<<<<<<< HEAD
//            else if ( the_var->getRevObject().isConvertibleTo( *it, convert_by_value ) != -1 )
//            {
//                // Fit by type conversion. For now, we also modify the type of the incoming variable wrapper.
//                RevObject* convertedObject = the_var->getRevObject().convertTo( *it );
//
//                RevPtr<RevVariable> valueVar = RevPtr<RevVariable>( new RevVariable(convertedObject,the_var->getName() ) );
//                return Argument( valueVar, arg.getLabel(), true );
//                
//            }
//        } // if (by-value)
//        else
//        {
//            if ( the_var->getRevObject().isType( *it ) )
//            {
//                // For now, change the required type of the incoming variable wrapper
//                the_var->setRequiredTypeSpec( *it );
//            
//                if ( isEllipsis() == false )
//                {
//                    return Argument( the_var, arg.getLabel(), evalType == BY_CONSTANT_REFERENCE );
//                }
//                else
//                {
//                    return Argument( the_var, arg.getLabel(), true );
//                }
//            
//            }
//            else if ( the_var->getRevObject().isConvertibleTo( *it, convert_by_value ) != -1  && (*it).isDerivedOf( the_var->getRequiredTypeSpec() ) )
//            {
//                // Fit by type conversion. For now, we also modify the type of the incoming variable wrapper.
//                RevObject* converted_object = the_var->getRevObject().convertTo( *it );
//                
//                RevPtr<RevVariable> the_new_var = NULL;
//                if ( the_var->getRevObject().isConstant() == true )
//                {
//                    the_new_var = the_var;
//                    the_new_var->replaceRevObject( converted_object );
//                    the_new_var->setRequiredTypeSpec( *it );
//                }
//                else
//                {
//                    the_new_var = RevPtr<RevVariable>( new RevVariable(converted_object, the_var->getName() ) );
//                }
//                
//                if ( !isEllipsis() )
//                {
//                    return Argument( the_new_var, arg.getLabel(), false );
//                }
//                else
//                {
//                    return Argument( the_new_var, arg.getLabel(), false );
//                }
//            }
//=======
//>>>>>>> bec51c4a9f3581046b7c8af121197e7ef0dae9d4
            else
            {
                // For now, change the required type of the incoming variable wrapper
                the_var->setRequiredTypeSpec( argTypeSpec );
            
                return Argument( the_var, arg.getLabel(), isEllipsis() or evalType == BY_CONSTANT_REFERENCE );
            }
        }
        else if ( the_var->getRevObject().isConvertibleTo( argTypeSpec, convert_by_value ) != -1 )
        {
            // Fit by type conversion. For now, we also modify the type of the incoming variable wrapper.
            RevObject* convertedObject = the_var->getRevObject().convertTo( argTypeSpec );

            RevPtr<RevVariable> the_new_var = nullptr;
            if ( the_var->getRevObject().isConstant() and not by_value )
            {
                the_new_var = the_var;
                the_new_var->replaceRevObject( convertedObject );
                the_new_var->setRequiredTypeSpec( argTypeSpec );
            }
            else
                the_new_var = RevPtr<RevVariable>( new RevVariable(convertedObject, the_var->getName() ) );
                
            return Argument( the_new_var, arg.getLabel(), false );
        }
        else
        {
            // Fit by type conversion function
            
            const TypeSpec& typeFrom = the_var->getRevObject().getTypeSpec();
            const TypeSpec& typeTo   = argTypeSpec;
            
            // create the function name
            std::string function_name = "_" + typeFrom.getType() + "2" + typeTo.getType();
                
            // Package arguments
            std::vector<Argument> args;
            Argument theArg = Argument( the_var, "arg" );
            args.push_back( the_var );
                
            auto env = Workspace::globalWorkspacePtr();
            
            if (auto orig_func = env->findFunction(function_name, args))
            {
                Function* func = orig_func->clone();

                // Allow the function to process the arguments
                func->processArguments( args );
            
                // Set the execution environment of the function
                func->setExecutionEnviroment( env );
                
                // Evaluate the function
                RevPtr<RevVariable> conversionVar = func->execute();
                
                // free the memory
                delete func;
                
                conversionVar->setHiddenVariableState( true );
                conversionVar->setRequiredTypeSpec( argTypeSpec );
                
                return Argument( conversionVar, arg.getLabel(), evalType == BY_CONSTANT_REFERENCE );
            }
        } 
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
        throw RbException() << "Cannot get default variable \"" << label << "\"";
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


double ArgumentRule::isArgumentValid( Argument &arg) const
{
    RevPtr<RevVariable> the_var = arg.getVariable();
    if ( the_var == NULL )
    {
        return -1;
    }
    
    bool convert_by_value = false;
    if (not the_var->getRevObject().hasDagNode() or the_var->getRevObject().getDagNode()->isConstant())
    {
        convert_by_value = true;
    }
    else if ( evalType == BY_VALUE || the_var->isWorkspaceVariable() || the_var->getRevObject().isConstant())
    {
        convert_by_value = true;
    }

    if ( nodeType == STOCHASTIC || nodeType == DETERMINISTIC )
    {
        convert_by_value = false;
    }
    
    if ( nodeType == STOCHASTIC && the_var->getRevObject().getDagNode()->getDagNodeType() != RevBayesCore::DagNode::STOCHASTIC )
    {
        return -1;
    }
    else if ( nodeType == DETERMINISTIC && the_var->getRevObject().getDagNode()->getDagNodeType() != RevBayesCore::DagNode::DETERMINISTIC )
    {
        return -1;
    }
    
    // we need to store and check all arg types
    std::vector<double> penalties;
    for ( auto& req_arg_type_spec: argTypeSpecs )
    {
        if ( the_var->getRevObject().isType( req_arg_type_spec ) )
        {
            return 0.0;
        }
            
        double penalty = -1;
        // make sure that we only perform type casting when the variable will not be part of a model graph
        if ( the_var->getRevObject().isConstant() == true )
        {
            penalty = the_var->getRevObject().isConvertibleTo( req_arg_type_spec, convert_by_value );
        }
        
        if ( penalty != -1 && req_arg_type_spec.isDerivedOf( the_var->getRequiredTypeSpec() ) )
        {
            penalties.push_back( penalty );
        }
        else if ( penalty != -1 && evalType == BY_VALUE )
        {
            penalties.push_back( penalty );
        }
        else if ( nodeType != STOCHASTIC )
        {
            
            const TypeSpec& typeFrom = the_var->getRevObject().getTypeSpec();
            const TypeSpec& typeTo   = req_arg_type_spec;
            
            // create the function name
            std::string function_name = "_" + typeFrom.getType() + "2" + typeTo.getType();
                
            // Package arguments
            std::vector<Argument> args;
            Argument theArg = Argument( the_var, "arg" );
            args.push_back( the_var );
                
            Environment& env = Workspace::globalWorkspace();
	    if (env.findFunction(function_name, args) != nullptr)
		return 0.1;
        }
    }
        
    // check which one was the best penalty
    double best_penalty = -1;
    for (size_t i=0; i<penalties.size(); ++i)
    {
        if ( penalties[i] != -1 )
        {
            if (best_penalty == -1 || best_penalty > penalties[i])
            {
                best_penalty = penalties[i];
            }
        }
    }
    
    
    return best_penalty;
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


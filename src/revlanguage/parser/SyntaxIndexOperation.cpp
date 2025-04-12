#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "RlContainer.h"
#include "Integer.h"
#include "IntegerPos.h"
#include "RevNullObject.h"
#include "RlMemberMethod.h"
#include "SyntaxIndexOperation.h"
#include "Workspace.h"
#include "Argument.h"
#include "DagNode.h"
#include "Environment.h"
#include "MethodTable.h"
#include "RbException.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "StringUtilities.h"
#include "SyntaxElement.h"

using namespace RevLanguage;


/** Construct from base variable (member object), identifier and index */
SyntaxIndexOperation::SyntaxIndexOperation(SyntaxElement* var, SyntaxElement* indx) : SyntaxElement(),
    index(indx),
    base_variable(var)
{

}


/** Deep copy constructor */
SyntaxIndexOperation::SyntaxIndexOperation(const SyntaxIndexOperation& x) : SyntaxElement(x)
{

    if ( x.base_variable != NULL )
    {
        base_variable = x.base_variable->clone();
    }
    else
    {
        base_variable = NULL;
    }

    if ( x.index != NULL )
    {
        index = x.index->clone();
    }
    else
    {
        index = NULL;
    }

}


/** Destructor deletes variable, identifier and index */
SyntaxIndexOperation::~SyntaxIndexOperation()
{

    delete base_variable;
    delete index;
}


/** Assignment operator */
SyntaxIndexOperation& SyntaxIndexOperation::operator=(const SyntaxIndexOperation& x)
{

    if (&x != this)
    {
        delete base_variable;
        delete index;


        SyntaxElement::operator=(x);


        if ( x.base_variable != NULL )
        {
            base_variable = x.base_variable->clone();
        }
        else
        {
            base_variable = NULL;
        }

        if ( x.index != NULL )
        {
            index = x.index->clone();
        }
        else
        {
            index = NULL;
        }

    }

    return (*this);
}



/** Clone syntax element */
SyntaxIndexOperation* SyntaxIndexOperation::clone () const
{

    return new SyntaxIndexOperation(*this);
}


/**
 * @brief Get semantic value (l-value)
 *
 * This function is similar to getValue(), but we only get the
 * slot and index of the variable. Another difference is that we
 * do not throw an error if the variable does not exist in the
 * frame; instead, we return a NULL pointer and set theSlot pointer
 * to NULL as well.
 */
RevPtr<RevVariable> SyntaxIndexOperation::evaluateLHSContent( const std::shared_ptr<Environment>& env, const std::string &varType)
{

    RevPtr<RevVariable> indexVar     = index->evaluateContent(env);
    RevPtr<RevVariable> theParentVar = base_variable->evaluateLHSContent(env, varType);
    
    // first we need to check if the parent variable is a deterministic vector
    if ( theParentVar->isVectorVariable() )
    {
        // everything is fine and we can add this element
    }
    else if ( theParentVar->getRevObject() == RevNullObject::getInstance() )
    {

    }
    else
    {
        RevObject& theParentObj = theParentVar->getRevObject();
        
        Container* c = dynamic_cast<Container*>( &theParentObj );
        if ( c != NULL )
        {
            if ( theParentObj.getDagNode()->isConstant() == false )
            {
                throw RbException("We cannot create a composite container from a non-constant container");
            }
            else if ( c->allowsModificationToCompositeContainer() == false )
            {
                throw RbException() << "An object of type '" << theParentObj.getType() << "' does not allow transformation into a composite container.";
            }
            else
            {
                // We need to set this to a vector variable before we call addIndex on it.
                theParentVar->setToVectorVariable();

                for ( size_t i=1; i<=c->size(); ++i)
                {
                    std::string elementIdentifier = theParentVar->getName() + "[" + std::to_string(i) + "]";
                
                    // first ensure that the environment has a slot named `elementIdentifier`
                    if ( !env->existsVariable( elementIdentifier ) )
                    {
                        // create a new variable for environment slot `elementIdentifier`
                        RevPtr<RevVariable> theElementVar = RevPtr<RevVariable>( new RevVariable( c->getElement(i-1) ) );
                        env->addVariable( elementIdentifier, theElementVar );
                        theElementVar->setName( elementIdentifier );
                    }

                    // then look up the element variable in the environment by its identifier
                    RevPtr<RevVariable> theElementVar  = env->getVariable( elementIdentifier );
                
                    // set the element variable as a hidden variable so that it doesn't show in ls()
                    theElementVar->setElementVariableState( true );
                    
                    theParentVar->addIndex( int(i) , theElementVar );
                }
            }
        }
        else
        {
            throw RbException() << "We cannot make a composite container from variable of type '" << theParentObj.getType() << "'.";
        }
    }

    // compute the index and internal name for this variable
    int idx = (int)static_cast<Integer&>( indexVar->getRevObject() ).getValue();
    std::string identifier = theParentVar->getName() + "[" + std::to_string(idx) + "]";

    // first ensure that we've got an entry with name `identifier` in the environment
    if ( not env->existsVariable( identifier ) )
    {
        // create a new variable called `identifier` if the environment doesn't have one.
        RevPtr<RevVariable> the_var = RevPtr<RevVariable>( new RevVariable( NULL ) );
        env->addVariable( identifier, the_var );
        the_var->setName( identifier );
    }

    // then look up the variable in the environment by its identifier
    RevPtr<RevVariable> the_var  = env->getVariable( identifier );
    
    // set this variable as an element variable; which is by default a hidden variable so that it doesn't show in ls()
    the_var->setElementVariableState( true );

    // mark the parent variable as a vector variable
    theParentVar->setToVectorVariable();
    theParentVar->addIndex( idx, the_var );

    return the_var;
}


/**
 * @brief Get semantic value (r-value)
 *
 * The variable can either be a member or a base variable. In the latter
 * case, its "base_variable" member is NULL. If the element is a base variable,
 * we get the semantic value of the element by looking it up in the frame.
 * If it is a member variable, we try to find it in the member variable
 * frame of a composite variable found by another SyntaxIndexOperation element.
 *
 * If dynamic == true, then we cannot compute a static index for indexed variables.
 * Instead, we need to deliver an index consisting of variables resulting
 * from dynamic evaluation of the index variables. These need to be put
 * in a dynamic lookup variable.
 */
RevPtr<RevVariable> SyntaxIndexOperation::evaluateContent( const std::shared_ptr<Environment>& env, bool dynamic)
{

    RevPtr<RevVariable> indexVar     = index->evaluateContent(env,dynamic);
    RevPtr<RevVariable> theParentVar = base_variable->evaluateContent( env );
    std::string identifier = theParentVar->getName() + "[" + indexVar->getRevObject().toString() + "]";

    RevPtr<RevVariable> the_var = NULL;

    // test whether we want to directly assess the variable or if we want to assess subelement of this container
    // if this variable already exists in the workspace
    if ( env->existsVariable( identifier ) )
    {
        // get the from the workspace
        the_var = env->getVariable( identifier );

    }
    else
    {


        try
        {
            // convert the value into a member object
            std::vector<Argument> args;
            args.push_back( Argument( theParentVar, "v" ) );
            args.push_back( Argument( indexVar, "index" ) );
            
            Function* f = Workspace::userWorkspace().getFunction("[]", args, false).clone();
            f->processArguments( args );
            the_var = f->execute();
            the_var->setName( identifier );
            the_var->setElementVariableState(true);
            
            delete f;
        }
        catch (RbException &e)
        {
            // We are trying to find a member function
        
            // Now we get a reference to the member object inside
            RevObject &theMemberObject = theParentVar->getRevObject();
            
            const MethodTable& mt = theMemberObject.getMethods();
            
            // convert the value into a member object
            std::vector<Argument> args;
            args.push_back( Argument( indexVar, "index" ) );
            
            Function* the_function = mt.getFunction( "[]", args, !dynamic ).clone();
            the_function->processArguments(args);
            
            MemberMethod* theMemberMethod = dynamic_cast<MemberMethod*>( the_function );
            if ( theMemberMethod != NULL )
            {
                theMemberMethod->setMemberObject( theParentVar );

                the_var = the_function->execute();
                the_var->setName( identifier );
                the_var->setElementVariableState(true);
            }
            else
            {
                delete the_function;
                throw e;
            }
            
        }

    }

    return the_var;

}


/**
 * Get the pointer to the internal base variable for this syntax element.
 */
SyntaxElement* SyntaxIndexOperation::getBaseVariable( void )
{
    // return the internal pointer
    return base_variable;
}


/**
 * Update the variable.
 * We need to refresh the composite variables so that the relationships are properly set.
 */
void SyntaxIndexOperation::updateVariable( Environment& env, const std::string &n )
{
    
    std::string varName = n;
    size_t pos = varName.find_last_of('[');
    if ( pos != std::string::npos)
    {
        std::string parentName = varName.substr(0,pos);
        
        if ( env.existsVariable(parentName) )
        {
            RevPtr<RevVariable> &parentVariable = env.getVariable(parentName);
            
//            const std::set<int>& indices = parentVariable->getElementIndices();
//            if ( indices.empty() )
//            {
//                throw RbException() << "Cannot create a vector variable with name '" << parentName << "' because it doesn't have elements.";
//            }
            size_t max_index = parentVariable->getMaxElementIndex();
            std::vector<Argument> args;
            for (size_t i = 1; i <= max_index; ++i)
            {
                std::string element_identifier = parentName + "[" + std::to_string(i) + "]";
                RevPtr<RevVariable>& elementVar = env.getVariable( element_identifier );
                // check that the element is not NULL
                if ( elementVar == NULL || elementVar->getRevObject() == RevNullObject::getInstance() )
                {
                    throw RbException() << "Cannot create vector variable with name '" << parentName << "' because element with name '" << element_identifier << "' is NULL." ;
                }
                args.push_back( Argument( elementVar ) );
            }
            Function* func = Workspace::userWorkspace().getFunction("v",args,false).clone();
            func->processArguments(args);
            
            // Evaluate the function (call the static evaluation function)
            RevPtr<RevVariable> funcReturnValue = func->execute();
            
            // free the memory of our copy
            delete func;
            
            parentVariable->replaceRevObject( funcReturnValue->getRevObject().clone() );
            
            SyntaxIndexOperation *parentExpression = dynamic_cast< SyntaxIndexOperation *>( base_variable );
            if ( parentExpression != NULL )
            {
                parentExpression->updateVariable(env, parentName);
            }
            
        }
        
    }
    
}


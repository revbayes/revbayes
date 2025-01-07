#include <cstddef>
#include <string>
#include <sstream>
#include <vector>
#include <memory>

#include "Argument.h"
#include "RbException.h"
#include "RevNullObject.h"
#include "RlFunction.h"
#include "TypeSpec.h"
#include "RevVariable.h"
#include "Workspace.h"
#include "DagNode.h"
#include "Environment.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "StringUtilities.h"

using namespace RevLanguage;

/** Constructor of empty RevVariable with specified type. */
RevVariable::RevVariable( const TypeSpec& ts, const std::string& n ) :
    name( n ),
    required_type_spec( ts )
{
    
}

/** Constructor of filled RevVariable (no type restrictions). */
RevVariable::RevVariable(RevObject *v, const std::string &n) :
    name( n ),
    required_type_spec( RevObject::getClassTypeSpec() )
{
    replaceRevObject( v );
}


/** Constructor of reference RevVariable (no type restrictions). */
RevVariable::RevVariable(const RevPtr<RevVariable>& refVar, const std::string &n) :
    name( n ),
    referenced_variable( refVar ),
    required_type_spec( RevObject::getClassTypeSpec() )
{
    
}


/** Copy constructor */
RevVariable::RevVariable(const RevVariable &v) :
    needs_building( v.needs_building ),
    is_element_var( v.is_element_var ),
    is_hidden_var( v.is_hidden_var ),
    vector_var_elements( v.vector_var_elements ),
    is_workspace_var( v.is_workspace_var ),
    name( v.name ),
    ref_count( 0 ),
    referenced_variable( v.referenced_variable ),
    rev_object( nullptr ),
    required_type_spec( v.required_type_spec )
{
    
    if ( v.rev_object )
    {
        replaceRevObject( v.rev_object->clone() );
    }
    else
    {
        rev_object = nullptr;
    }
    
}


RevVariable::~RevVariable( void )
{
    
    if ( not isReferenceVariable() and rev_object != NULL )
    {
        delete rev_object;
    }
    
}


RevVariable& RevVariable::operator=(const RevVariable &v)
{
    if ( this != &v )
    {
        
        name                = v.name;
        required_type_spec  = v.required_type_spec;
        is_element_var      = v.is_element_var;
        is_hidden_var       = v.is_hidden_var;
        vector_var_elements = v.vector_var_elements;
        is_workspace_var    = v.is_workspace_var;
        needs_building      = v.needs_building;
        referenced_variable = v.referenced_variable;
        
        if ( isReferenceVariable() )
        {
            rev_object = NULL;
        }
        else
        {
            if ( rev_object != NULL )
            {
                delete rev_object;
                rev_object = NULL;
            }
            
            if ( v.rev_object != NULL )
            {
                replaceRevObject( v.rev_object->clone() );
            }
            
        }
    }
    
    return *this;
}


/** Resize the vector to include this index. */
void RevVariable::addIndex(int idx, const RevPtr<RevVariable>& element)
{
    assert( isVectorVariable() );

    if (idx < 1)
        throw RbException()<<"Index "<<idx<<" not allowed in '"<<name<<"["<<idx<<"]'.  The first element should be '"<<name<<"[1]'";

    if (idx > vector_var_elements->size())
        vector_var_elements->resize(idx);

    (*vector_var_elements)[idx-1] = element;
    
    needs_building = true;
}


/* Clone RevVariable and RevVariable */
RevVariable* RevVariable::clone( void ) const
{
    
    return new RevVariable( *this );
}


/* Decrement the reference count. */
size_t RevVariable::decrementReferenceCount( void ) const
{
    
    ref_count--;
    
    return ref_count;
}


size_t RevVariable::getMaxElementIndex( void ) const
{
    assert(isVectorVariable());
    return vector_var_elements->size();
}


const std::string& RevVariable::getName( void ) const
{
    
    return name;
}



/* Get the reference count for this instance. */
size_t RevVariable::getReferenceCount(void) const
{
    return ref_count;
}


/* Get the value of the RevVariable */
RevObject& RevVariable::getRevObject(void) const
{
    
    if ( isReferenceVariable() )
    {
        return referenced_variable->getRevObject();
    }
    
    if ( isVectorVariable() && needs_building )
    {
        needs_building = false;

        std::vector<Argument> args;
        int i=1;
        for (auto& element_var: *vector_var_elements)
        {
            // check that the element is not NULL
            if ( element_var == NULL || element_var->getRevObject() == RevNullObject::getInstance() )
            {
                std::string element_identifier = name + "[" + std::to_string(i) + "]";
                throw RbException()<<"Cannot create vector variable with name '"<<name
                                   <<"' because element with name '"<<element_identifier<<"' is NULL.";

                /* NOTE: This scenario can occur if we define x[4] but not x[3] and then try to print x.
                 *       When x[4] is defined, vector_var_elements is padded with nullptrs to length 4.
                 *       But only x[4] is then set to a non-NULL value.
                 */
            }

            args.push_back( Argument( element_var ) );
            i++;
        }

        // Sebastian: We absolutely must use dynamic functions because we may want to use this vector in a dynamic, i.e., determistic function
        // if this is not dynamic, then the vector is always converted into a constant variable and does not allow it's elements to change.
        bool dynamic = true;
        std::unique_ptr<Function> func( Workspace::userWorkspace().getFunction("v", args, dynamic).clone() );
        func->processArguments(args, dynamic);
        
        // Evaluate the function (call the static evaluation function)
        RevPtr<RevVariable> func_return_value = func->execute();
        
        const_cast<RevVariable*>(this)->replaceRevObject( func_return_value->getRevObject().clone() );
    }
    
    if ( rev_object == NULL )
    {
        return RevNullObject::getInstance();
    }
    
    return *rev_object;
}


/**
 * Get the required type specs for values stored inside this RevVariable.
 * We return our own type specification even if we reference another
 * RevVariable. By reassignment, we receive the new value, so it is
 * more important what our type spec is.
 */
const TypeSpec& RevVariable::getRequiredTypeSpec(void) const
{
    return required_type_spec;
}

RevPtr<RevVariable> RevVariable::getElementVariable(int i) const
{
    assert(isVectorVariable());
    return (*vector_var_elements)[i];
}

/* Increment the reference count for this instance. */
void RevVariable::incrementReferenceCount( void ) const
{
    ref_count++;
}



/**
 * Is the RevVariable or any of its members (upstream DAG nodes) assignable, that is,
 * modifiable by the user? For them to be assignable, they have to be named, otherwise
 * there is no chance for the user to change them.
 */
bool RevVariable::isAssignable( void ) const
{
    // Check if we are assignable
    if ( name != "" )
        return true;
    
    // Ask our object for assignable upstream RevVariables
    if ( rev_object != NULL && rev_object->isAssignable() )
        return true;
    
    // No possibility left to modify us
    return false;
}


/**
 * Return the internal flag signalling whether the RevVariable is an element of a vector, e.g., x[1] would be.
 */
bool RevVariable::isElementVariable( void ) const
{
    return is_element_var;
}


/** Return the internal flag signalling whether the RevVariable is currently a hidden RevVariable. Hidden RevVariables will not show in the ls() function. */
bool RevVariable::isHiddenVariable( void ) const
{
    return is_hidden_var;
}


/** Return the internal flag signalling whether the RevVariable is currently a reference RevVariable */
bool RevVariable::isReferenceVariable( void ) const
{
    return referenced_variable;
}


/** Return the internal flag signalling whether the RevVariable is currently a vector RevVariable, that is, should be computed by x := v(x[1],...) */
bool RevVariable::isVectorVariable( void ) const
{
    return bool(vector_var_elements);
}


/** Return the internal flag signalling whether the RevVariable is currently a workspace (control) RevVariable */
bool RevVariable::isWorkspaceVariable( void ) const
{
    if ( isReferenceVariable() )
    {
        return referenced_variable->isWorkspaceVariable();
    }
    else
    {
        return is_workspace_var;
    }
}


/** Make this RevVariable a reference to another RevVariable. Make sure we delete any object we held before. */
void RevVariable::makeReference(const RevPtr<RevVariable>& refVar)
{
    if ( not isReferenceVariable ())
    {
        if ( rev_object != NULL )
        {
            delete rev_object;
        }
        
        rev_object = NULL;
        is_workspace_var = false;
    }
    
    referenced_variable = refVar;
}


/* Print value of the RevVariable RevVariable */
void RevVariable::printValue(std::ostream& o, bool toScreen) const
{
    
    if ( rev_object == NULL )
    {
        o << "NULL";
    }
    else
    {
        rev_object->printValue( o, toScreen );
    }
    
}


/** Replace the RevVariable and update the DAG structure. */
void RevVariable::replaceRevObject( RevObject *newValue )
{
    if ( isReferenceVariable() )
    {
        referenced_variable = NULL;
    }
//
//    // Make sure default assignment is not a workspace (control) RevVariable assignment
//    is_workspace_var = false;
    
    if ( rev_object != NULL )
    {
        
        if ( rev_object->isModelObject() && rev_object->getDagNode() != NULL && newValue->isModelObject() )
        {
            rev_object->getDagNode()->replace( newValue->getDagNode() );
        }
        
    }
    
    delete rev_object;
    rev_object = newValue;
    
    if ( rev_object != NULL )
    {
        rev_object->setName( name );
    }
    
    try
    {
        if ( rev_object != NULL && rev_object->isModelObject() == true && rev_object->getDagNode() != NULL )
        {
            rev_object->getDagNode()->setHidden( is_hidden_var );
            rev_object->getDagNode()->setElementVariable( is_element_var );
        }
    }
    catch (RbException &e)
    {
        // do nothing
    }
    
}




/**
 * Set whether this RevVariable is an element of a vector RevVariable.
 * All element RevVariable are also hidden.
 * Throw an error if the RevVariable is a reference RevVariable.
 * If so, you need to set the Rev object first, and then set the hidden RevVariable flag.
 */
void RevVariable::setElementVariableState(bool flag)
{
    
    if ( isReferenceVariable() )
    {
        throw "A reference RevVariable cannot be made a hidden RevVariable";
    }
    
    is_element_var = flag;
    
    try
    {
        
        if ( rev_object != NULL && rev_object->getDagNode() != NULL )
        {
            rev_object->getDagNode()->setElementVariable( flag );
        }
        
    }
    catch (RbException &e)
    {
        // do nothing
    }
}


/**
 * Set whether this RevVariable is a hidden RevVariable. Throw an error if the RevVariable
 * is a reference RevVariable. If so, you need to set the Rev object first, and then set
 * the hidden RevVariable flag.
 */
void RevVariable::setHiddenVariableState(bool flag)
{
    if ( isReferenceVariable() )
    {
        throw "A reference RevVariable cannot be made a hidden RevVariable";
    }
    
    is_hidden_var = flag;
    
    try
    {
        
        if ( rev_object != NULL && rev_object->getDagNode() != NULL )
        {
            rev_object->getDagNode()->setHidden( flag );
        }
        
    }
    catch (RbException &e)
    {
        // do nothing
    }

}


/**
 * Set whether this RevVariable is a vector RevVariable. Throw an error if the RevVariable
 * is a reference RevVariable. If so, you need to set the Rev object first, and then set
 * the vector RevVariable flag.
 */
void RevVariable::setToVectorVariable()
{
    if ( isReferenceVariable() )
    {
        throw "A reference RevVariable cannot be made a vector RevVariable";
    }
    
    if (not vector_var_elements)
    {
        vector_var_elements = std::vector<RevPtr<RevVariable>>();
        needs_building = true;
    }
}


/**
 * Check whether this RevVariable is a workspace (control) RevVariable. Throw an error if the RevVariable
 * is a reference RevVariable. If so, you need to set the Rev object first, and then set
 * the workspace (control) RevVariable flag.
 */
void RevVariable::setWorkspaceVariableState(bool flag)
{
    if ( isReferenceVariable() and flag == true)
    {
        throw "A reference RevVariable cannot be made a workspace (control) RevVariable";
    }
    
    is_workspace_var = flag;
}


/** Set the name of the RevVariable */
void RevVariable::setName(std::string const &n)
{
    
    name = n;
    if ( rev_object != NULL )
    {
        rev_object->setName( n );
    }
    
}



/**
 * We set here the required value type spec. An error is thrown if the
 * current Rev object of the variable, if any, is not of the specified type.
 */
void RevVariable::setRequiredTypeSpec(const TypeSpec &ts)
{
    
    const RevObject& theObject = this->getRevObject();
    if ( theObject != RevNullObject::getInstance() )
    {
        if ( theObject.isType( ts ) == false )
        {
            throw RbException()<<"setRequiredTypeSpec: failed setting object of type "<<theObject.getType()<<" to "<<ts.getType();
        }
    }
    
    required_type_spec = ts;
}

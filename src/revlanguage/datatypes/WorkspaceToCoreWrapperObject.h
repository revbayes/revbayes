#ifndef WorkspaceToCoreWrapperObject_H
#define WorkspaceToCoreWrapperObject_H

#include "WorkspaceObject.h"

namespace RevLanguage {
    
    /**
     * @file
     * This file contains the declaration of the Rev language WorkspaceToCoreWrapperObject, which is
     * the Rev abstract base class for all language objects wrapping core datatypes
     * that can be stored inside DAG nodes and hence used in Workspace graphs.
     *
     * @brief Declaration of WorkspaceToCoreWrapperObject
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2014-10-21
     */
    template <typename rbType>
    class WorkspaceToCoreWrapperObject : public WorkspaceObject {
        
    public:
        virtual                                            ~WorkspaceToCoreWrapperObject(void);                             //!< Virtual destructor
        
        WorkspaceToCoreWrapperObject&                       operator=(const WorkspaceToCoreWrapperObject& x);               //!< Assignment operator
        
        // the value type definition
        typedef rbType valueType;
        
        // Basic utility functions you have to override (also getClassTypeSpec()!)
        virtual WorkspaceToCoreWrapperObject<rbType>*       clone(void) const = 0;                                          //!< Clone object
//        virtual void                                        constructInternalObject(void) = 0;                              //!< We construct the a new internal object.
        static const std::string&                           getClassType(void);                                             //!< Get Rev type
        static const TypeSpec&                              getClassTypeSpec(void);                                         //!< Get class type spec
        virtual const TypeSpec&                             getTypeSpec(void) const = 0;                                    //!< Get the type spec of the instance
        
        
        // Getters and setters
        rbType&                                             getValue(void) const;                                           //!< Get value
        void                                                setValue(const rbType& x);                                      //!< Set value
        virtual void                                        printValue(std::ostream& o, bool user) const;                                                      //!< Print value (for user)
        virtual void                                        printValue(std::ostream& o) const { printValue(o, true); };

    protected:
        WorkspaceToCoreWrapperObject(void);
        WorkspaceToCoreWrapperObject(rbType *v);
        WorkspaceToCoreWrapperObject(const WorkspaceToCoreWrapperObject &v);
        
        rbType*                                             value;
    };
    
}


#include "TypeSpec.h"

template <typename rbType>
RevLanguage::WorkspaceToCoreWrapperObject<rbType>::WorkspaceToCoreWrapperObject() : WorkspaceObject(),
    value( NULL )
{
    
}



template <typename rbType>
RevLanguage::WorkspaceToCoreWrapperObject<rbType>::WorkspaceToCoreWrapperObject(rbType *v) : WorkspaceObject(),
    value( v )
{
    
}



template <typename rbType>
RevLanguage::WorkspaceToCoreWrapperObject<rbType>::WorkspaceToCoreWrapperObject(const WorkspaceToCoreWrapperObject &v) : WorkspaceObject( v ),
    value( NULL )
{
    
    if ( v.value != NULL )
    {
        value = v.value->clone();
    }
    
}



template <typename rbType>
RevLanguage::WorkspaceToCoreWrapperObject<rbType>::~WorkspaceToCoreWrapperObject()
{
    
    delete value;
}


template <typename rbType>
RevLanguage::WorkspaceToCoreWrapperObject<rbType>& RevLanguage::WorkspaceToCoreWrapperObject<rbType>::operator=(const WorkspaceToCoreWrapperObject &v) {
    
    if ( this != &v )
    {
        WorkspaceObject::operator=( v );
        
        // free the memory
        delete value;
        value = NULL;
        
        // create own copy
        if ( v.value != NULL ) {
            value = v.value->clone();
        }
    }
    
    return *this;
}


template <typename rbType>
const std::string& RevLanguage::WorkspaceToCoreWrapperObject<rbType>::getClassType(void)
{
    
    static std::string rev_type = "WorkspaceToCoreWrapperObject";
    
    return rev_type;
}


/** Get class type spec describing type of object */
template <typename rlType>
const RevLanguage::TypeSpec& RevLanguage::WorkspaceToCoreWrapperObject<rlType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceObject::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * This function allows you to access the internal value of the
 * workspace object. Because we do not need to guard our value
 * against inappropriate modification, unlike model objects, we
 * give out a non-const reference to it allowing the caller to
 * modify the object.
 */
template <typename rbType>
rbType& RevLanguage::WorkspaceToCoreWrapperObject<rbType>::getValue( void ) const
{
    return *value;
}


template <typename rbType>
void RevLanguage::WorkspaceToCoreWrapperObject<rbType>::setValue(const rbType &x) {
    
    // free memory
    delete value;
    value = new rbType( x );
}


template <typename rbType>
void RevLanguage::WorkspaceToCoreWrapperObject<rbType>::printValue(std::ostream &o, bool user) const
{
    
    o << value;
}



#endif

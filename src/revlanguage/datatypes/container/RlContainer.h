#ifndef Container_H
#define Container_H

#include "RevObject.h"
#include "TypeSpec.h"

#include <iostream>
#include <vector>

namespace RevLanguage {
    
    /**
     * @brief Generic container base class
     *
     * Container is the interface for all generic Rev containers, that is, Rev objects with a type of the
     * form element type followed by one or more pairs of square brackets marking the dimension of the
     * container. For instance, a generic Rev vector of Real elements will have type "Real[]", and a
     * vector of probabilities will have type "Probability[]". The actual generic vector types used
     * in Rev are classes deriving directly from either ModelVector (for ModelObject elements) or
     * from WorkspaceVector (for WorkspaceObject elements).
     *
     * More specialized containers derive from the generic container classes, and thus inherit all
     * of the container properties.
     *
     * Generic containers support assignment to individual elements by default, but derived classes
     * defining specialized classes typically override this behavior by appropriate modification of
     * the element access functions.
     */
    class Container {
        
    public:
        virtual                                    ~Container(void) {}
        
        // Container functions you have to override
        virtual void                                push_back(const RevObject &x) = 0;                                  //!< Append element to end
        virtual RevObject*                          getElement(size_t idx) const = 0;                                   //!< Get element variable (vector of indices)
        virtual size_t                              size(void) const = 0;                                               //!< Get number of elements in container
        
        virtual bool                                allowsModificationToCompositeContainer(void) const { return true;}  //!< Does an object of this type allow transformation into a composite container?

    protected:
        Container(void) {}                                                                                              //!< Default constructor
      //Container&                                  operator=(const Container& x);                                      //!< Assignment operator
    };
    
}

#endif


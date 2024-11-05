#ifndef RbContainer_H
#define RbContainer_H

#include <iostream>

namespace RevBayesCore {
    
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
        virtual                                     ~Container(void) {}
        
        // Container functions you have to override
        virtual size_t                              size(void) const = 0;                                                                                   //!< Get number of elements in container
        virtual void                                printElement(std::ostream &o, size_t i, std::string sep="\t", std::int64_t l=-1, bool left=true) const = 0;                  //!< Print the i-th element
        
        
    protected:
        Container(void) {}                                                                                                                                  //!< Default constructor
        
        //        Container&                                  operator=(const Container& x);                                                                //!< Assignment operator
        
    };
    
}

#endif


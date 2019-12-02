#ifndef RevAbstractType_H
#define RevAbstractType_H

#include "WorkspaceObject.h"

#include <ostream>
#include <string>


namespace RevLanguage {
    
    class TypeSpec;
    
    /**
     * @brief RevAbstractType: Object representing abstract Rev types
     *
     * This class is used to represent abstract Rev types in contexts where an object is
     * needed, such as in type checking.
     *
     */
    class RevAbstractType : public WorkspaceObject {
        
    public:
                                            RevAbstractType(const TypeSpec& TypeSpec, RevObject* obj);              //!< Standard constructor
                                            RevAbstractType(const RevAbstractType& n);                              //!< Copy constructor
        virtual                            ~RevAbstractType(void) {}                                                //!< Virtual destructor just in case
        
        // Basic utility functions you have to override (also getClassTypeSpec()!)
        RevAbstractType*                    clone(void) const;                                                      //!< Clone object
        const TypeSpec&                     getTypeSpec(void) const;                                                //!< Get the type spec for this instance
        bool                                isAbstract(void) const;                                                 //!< Return true because we are abstract
        
        // RevAbstractType function
        RevObject*                          makeExampleObject(void) const;                                          //!< Make new example object instance
        void                                printValue(std::ostream& o, bool user) const;                           //!< Print value
        
    private:
        const TypeSpec&                     myTypeSpec;                                                             //!< Abstract type specification
        RevObject*                          exampleObject;                                                          //!< Example of non-abstract derived class
    };

    // Global functions using the class
    std::ostream&                           operator<<(std::ostream& o, const RevAbstractType& x);                    //!< Overloaded output operator
    
}


#endif

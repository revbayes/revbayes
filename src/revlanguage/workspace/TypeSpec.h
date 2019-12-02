#ifndef TypeSpec_H
#define TypeSpec_H

#include <string>


namespace RevLanguage {

    /**
     * @brief TypeSpec: holding type specifications
     *
     * A type specification is a vector of strings describing
     * the inheritance for a Rev type. For instance, the Rev
     * type "Natural" has the type specification:
     *
     * [ "RevObject", "Integer", "Natural" ]
     *
     * and the type "Simplex" has the type specification:
     *
     * [ "RevObject", "Container", "ModelContainer", "Real[]", "Simplex" ]
     *
     * The TypeSpec class is simply a way of organizing the type
     * specifications, and providing some related type checking
     * functionality.
     *
     * The Rev type specification is inspired by the "class" attribute
     * in R. The type specification can be accessed using the "class"
     * function in Rev, while the type can be accessed using the "type"
     * function.
     *
     * We do not need independent copies of the parent or the element
     * type spec pointers because of the way type specifications are used
     * in Rev. Rev types and their type specifications are guaranteed to
     * be created in a top-down manner, and remain in place throughout the
     * duration of the run because of the use of the static type spec
     * functions in all Rev classes.
     *
     * We give out pointers to parentTypeSpec and elementTypeSpec here
     * without any checking because this is performance-critical code.
     * We assume that the caller uses the class appropriately.
     */
    class TypeSpec {
    
    public:
                                    TypeSpec(const std::string& objType, const TypeSpec* par);                              //!< Standard constructor
                                    TypeSpec(const std::string& objType, const TypeSpec* par, const TypeSpec* elemType);    //!< Constructor for container object
        
                                    // Operators
        bool                        operator==(const TypeSpec& x) const;                                            //!< Equals operator
        bool                        operator!=(const TypeSpec& x) const;                                            //!< Not equals operator

                                    // Regular functions
        const std::string&          getElementType(void) const;                                                     //!< Get element type, if container
        const TypeSpec*             getElementTypeSpec(void) const;                                                 //!< Get element type spec, if container
        const std::string&          getParentType(void) const;                                                      //!< Get parent type
        const TypeSpec*             getParentTypeSpec(void) const;                                                  //!< Get parent type spec
        const std::string&          getType(void) const;                                                            //!< Get object type
        bool                        isDerivedOf(const TypeSpec& x) const;                                           //!< Test whether the type is the same or derived from x

    private:
        
        // Helper functions
        const TypeSpec*             getElementTypeSpecFromClass(void) const;                                        //!< Make element type spec, if container, using type spec hierarchy

        // Member variables
        const TypeSpec*             element;                                                                        //!< Pointer to element type spec, if any
        const TypeSpec*             parent;                                                                         //!< Pointer to the parent type
        std::string                 type;                                                                           //!< The type
    };

    // Global functions using the class
    std::ostream&                   operator<<(std::ostream& o, const TypeSpec& x);                                 //!< Overloaded output operator

}


#endif

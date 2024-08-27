#ifndef Integer_H
#define Integer_H

#include "ModelObject.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class Real;

    class Integer : public ModelObject<long> {

    public:
        Integer(void);                                                                                          //!< Default constructor
        Integer(RevBayesCore::TypedDagNode<long> *v);                                                           //!< Constructor from DAG node
        Integer(long v);                                                                                        //!< Constructor from long
//        Integer(unsigned long v);                                                                               //!< Constructor from unsigned long

        // Basic operator functions
        virtual RevObject*              add(const RevObject &rhs) const;                                        //!< Addition operator used for example in '+=' statements
        Integer*                        add(const Integer &rhs) const;                                          //!< Addition operator used for example in '+=' statements
        Real*                           add(const Real &rhs) const;                                             //!< Addition operator used for example in '+=' statements
        void                            decrement(void);                                                        //!< Decrement operator used for example in 'a--' statements
        virtual RevObject*              divide(const RevObject &rhs) const;                                     //!< Division operator used for example in '/=' statements
        Real*                           divide(const Integer &rhs) const;                                       //!< Division operator used for example in '/=' statements
        Real*                           divide(const Real &rhs) const;                                          //!< Division operator used for example in '/=' statements
        void                            increment(void);                                                        //!< Increment operator used for example in 'a++' statements
        virtual RevObject*              multiply(const RevObject &rhs) const;                                   //!< Multiplication operator used for example in '*=' statements
        Integer*                        multiply(const Integer &rhs) const;                                     //!< Multiplication operator used for example in '*=' statements
        Real*                           multiply(const Real &rhs) const;                                        //!< Multiplication operator used for example in '*=' statements
        virtual RevObject*              subtract(const RevObject &rhs) const;                                   //!< Subtraction operator used for example in '-=' statements
        Integer*                        subtract(const Integer &rhs) const;                                     //!< Subtraction operator used for example in '-=' statements
        Real*                           subtract(const Real &rhs) const;                                        //!< Subtraction operator used for example in '-=' statements

        // Basic utility functions
        virtual Integer*                clone(void) const;                                                      //!< Clone object
        virtual RevObject*              convertTo(const TypeSpec& type) const;                                  //!< Convert to type
        static const std::string&       getClassType(void);                                                     //!< Get Rev type
        static const TypeSpec&          getClassTypeSpec(void);                                                 //!< Get class type spec
        virtual const TypeSpec&         getTypeSpec(void) const;                                                //!< Get language type of the object
        virtual double                  isConvertibleTo(const TypeSpec& type, bool once) const;                 //!< Is convertible to type?

        template<class rbTypeTo> RevObject* convertTo() const;
        
        std::string                     getGuiName(void) { return "Integer"; }
        std::string                     getGuiUnicodeSymbol(void) { return "Z"; }
        std::string                     getGuiInfo(void) { return ""; }
    };

    template<class rbTypeTo>
    RevObject* Integer::convertTo() const
    {   
        typedef typename rbTypeTo::valueType valueType;
        Func__conversion<Integer,rbTypeTo>* rlFunc = new Func__conversion<Integer,rbTypeTo>();
        RevBayesCore::TypeConversionFunction<Integer::valueType,valueType>* func = new RevBayesCore::TypeConversionFunction<Integer::valueType,valueType>(dag_node);
        DeterministicNode<valueType>* newnode = new DeterministicNode<valueType>(dag_node->getName() + "2" + rbTypeTo::getClassType(), func, rlFunc);
        return new rbTypeTo(newnode);
    }
    
}

#endif


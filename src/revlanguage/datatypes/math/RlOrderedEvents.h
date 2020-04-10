/*
 * RlOrderedEvents.h
 *
 *  Created on: Apr 9, 2020
 *      Author: mrmay
 */

#ifndef MATH_RlOrderedEvents_H_
#define MATH_RlOrderedEvents_H_

#include "RlOrderedEvents.h"

#include <string>

#include "ModelObject.h"
#include "OrderedEvents.h"

namespace RevBayesCore { class RbHelpReference; }


namespace RevLanguage {
class Argument;
class RevVariable;
class TypeSpec;

	template<class rlType>
	class RlOrderedEvents : public ModelObject<RevBayesCore::OrderedEvents<typename rlType::valueType> > {

	public:

		RlOrderedEvents(void);                                                                                    //!< Constructor requires character type
		RlOrderedEvents(RevBayesCore::OrderedEvents<typename rlType::valueType> *v);                              //!< Constructor requires character type
		RlOrderedEvents(const RevBayesCore::OrderedEvents<typename rlType::valueType> &v);                        //!< Constructor requires character type
		RlOrderedEvents(RevBayesCore::TypedDagNode<RevBayesCore::OrderedEvents<typename rlType::valueType> > *n); //!< Constructor requires character type

        typedef RevBayesCore::OrderedEvents<typename rlType::valueType> valueType;

        RlOrderedEvents<rlType>*                     clone(void) const;                                                                  //!< Clone object
        static const std::string&                    getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&                       getClassTypeSpec(void);                                                             //!< Get class type spec
        const TypeSpec&                              getTypeSpec(void) const;                                                            //!< Get language type of the object

        // Member method inits
        virtual RevPtr<RevVariable>                  executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f); //!< Map member methods to internal functions

        std::string                                  getGuiName(void) { return "OrderedEvents"; }
        std::string                                  getGuiUnicodeSymbol(void) { return ""; }
        std::string                                  getGuiInfo(void) { return ""; }

    protected:

        void                                        initMethods(void);


	};

} /* namespace RevLanguage */


#include "TypeSpec.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevVariable.h"


/** Default constructor */
template <class rlType>
RevLanguage::RlOrderedEvents<rlType>::RlOrderedEvents(void) : ModelObject<RevBayesCore::OrderedEvents<typename rlType::valueType> >()
{
    initMethods();
}

/** Construct from core OrderedEvents */
template <class rlType>
RevLanguage::RlOrderedEvents<rlType>::RlOrderedEvents(RevBayesCore::OrderedEvents<typename rlType::valueType> *c) : ModelObject<RevBayesCore::OrderedEvents<typename rlType::valueType> >( c )
{
    initMethods();
}

/** Construct from core OrderedEvents */
template <class rlType>
RevLanguage::RlOrderedEvents<rlType>::RlOrderedEvents(const RevBayesCore::OrderedEvents<typename rlType::valueType> &t) : ModelObject<RevBayesCore::OrderedEvents<typename rlType::valueType> >( new RevBayesCore::OrderedEvents<typename rlType::valueType>( t ) )
{
    initMethods();
}

/** Construct from DAG node */
template <class rlType>
RevLanguage::RlOrderedEvents<rlType>::RlOrderedEvents(RevBayesCore::TypedDagNode<RevBayesCore::OrderedEvents<typename rlType::valueType> > *n) : ModelObject<RevBayesCore::OrderedEvents<typename rlType::valueType> >( n )
{
    initMethods();
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <class rlType>
RevLanguage::RlOrderedEvents<rlType>* RevLanguage::RlOrderedEvents<rlType>::clone(void) const
{
    return new RevLanguage::RlOrderedEvents<rlType>(*this);
}

/* Map calls to member methods */
template <class rlType>
RevLanguage::RevPtr<RevLanguage::RevVariable> RevLanguage::RlOrderedEvents<rlType>::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    return ModelObject<RevBayesCore::OrderedEvents<typename rlType::valueType> >::executeMethod( name, args, found );
}

/** Get Rev type of object */
template <class rlType>
const std::string& RevLanguage::RlOrderedEvents<rlType>::getClassType(void)
{
    static std::string rev_type = "OrderedEvents";

    return rev_type;
}

/** Get class type spec describing type of object */
template <class rlType>
const TypeSpec& RevLanguage::RlOrderedEvents<rlType>::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );

    return rev_type_spec;
}

/** Get type spec */
template <class rlType>
const TypeSpec& RevLanguage::RlOrderedEvents<rlType>::getTypeSpec( void ) const
{
    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}

/**
 * Initialize the member methods.
 */
template <class rlType>
void RevLanguage::RlOrderedEvents<rlType>::initMethods( void )
{
}











#endif /* RlOrderedEvents_H_ */

/*
 * RlOrderedEventTimes.h
 *
 *  Created on: Apr 9, 2020
 *      Author: mrmay
 */

#ifndef MATH_RLORDEREDEVENTTIMES_H_
#define MATH_RLORDEREDEVENTTIMES_H_

#include "ModelObject.h"
#include "OrderedEventTimes.h"

namespace RevBayesCore { class RbHelpReference; }


namespace RevLanguage {
class Argument;
class RevVariable;
class TypeSpec;

	class RlOrderedEventTimes : public ModelObject<RevBayesCore::OrderedEventTimes> {

	public:

		RlOrderedEventTimes(void);                                                                                    //!< Constructor requires character type
		RlOrderedEventTimes(RevBayesCore::OrderedEventTimes *v);                                                      //!< Constructor requires character type
		RlOrderedEventTimes(const RevBayesCore::OrderedEventTimes &v);                                                //!< Constructor requires character type
		RlOrderedEventTimes(RevBayesCore::TypedDagNode<RevBayesCore::OrderedEventTimes> *n);                          //!< Constructor requires character type

        typedef RevBayesCore::OrderedEventTimes valueType;

        RlOrderedEventTimes*                        clone(void) const;                                                                  //!< Clone object
        static const std::string&                   getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                             //!< Get class type spec
        const TypeSpec&                             getTypeSpec(void) const;                                                            //!< Get language type of the object

        // Member method inits
        virtual RevPtr<RevVariable>                 executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f); //!< Map member methods to internal functions

        std::string                                 getGuiName(void) { return "OrderedEventTimes"; }
        std::string                                 getGuiUnicodeSymbol(void) { return ""; }
        std::string                                 getGuiInfo(void) { return ""; }

    protected:

        void                                        initMethods(void);


	};

} /* namespace RevLanguage */

#endif /* MATH_RLORDEREDEVENTTIMES_H_ */

#ifndef RlTreeVector_H
#define RlTreeVector_H

#include <string>
#include <vector>
#include <iosfwd>

#include "ModelObject.h"
#include "TreeVector.h"
#include "Tree.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "RevPtr.h"
#include "RlConstantNode.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class RbHelpReference; }


namespace RevLanguage {
class Argument;
class RevVariable;
class TypeSpec;

    /**
     * @file
     * This file contains the declaration of a TreeVector, which is
     * the class that contains the name of a TreeVector with associated species and date.
     *
     * @brief Declaration of TreeVector
     *
     * (c) Copyright 2009-
     * @date Last modified: $Date: $
     * @author The RevBayes Development Core Team
     * @license GPL version 3
     *
     *
     */

    class TreeVector : public ModelObject<RevBayesCore::TreeVector> {

    public:
        TreeVector(void);                                                                        //!< Constructor requires character type
        TreeVector(RevBayesCore::TreeVector *v);                                                      //!< Constructor requires character type
        TreeVector(const RevBayesCore::TreeVector &v);                                                //!< Constructor requires character type
        TreeVector(RevBayesCore::TypedDagNode<RevBayesCore::TreeVector> *n);                          //!< Constructor requires character type

        typedef RevBayesCore::TreeVector         valueType;

        // Basic utility functions
        TreeVector*                            clone(void) const;                                                                  //!< Clone object
        static const std::string&                   getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                             //!< Get class type spec
        const TypeSpec&                             getTypeSpec(void) const;                                                            //!< Get language type of the object

        // Member method inits
        virtual RevPtr<RevVariable>                 executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f); //!< Map member methods to internal functions

        std::string                                 getGuiName(void) { return "TreeVector"; }
        std::string                                 getGuiUnicodeSymbol(void) { return "OTU"; }
        std::string                                 getGuiInfo(void) { return ""; }

    protected:

        void                                        initMethods(void);


    };
}

#endif

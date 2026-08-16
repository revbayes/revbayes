/*!
 * @file This file contains utilities for the RevLanguage, including type conversion functions.
 */

#ifndef RlUtil_H
#define RlUtil_H

#include "RevNullObject.h"
#include "RevObject.h"
#include "ModelObject.h"

#include <string>

namespace RevLanguage {
    
    namespace RlUtils {

        // Empty return type spec
        static const TypeSpec& Void    = TypeSpec( "void", NULL );

        class RlTypeConverter {
            public:
            static RevObject*                toReal(double x);
            static RevObject*                toString(const std::string &x);
            template<class rbTypeFrom, class rbTypeTo> static RevObject* convertTo(const ModelObject<typename rbTypeFrom::valueType>* input);
        };

        /**
         * Auxiliary function (templated) to convert to another type.
         * We create a deterministic node tied to the source node by a type conversion function, which will
         * handle updates. This also applies to constant source nodes so that the conversion does not change the
         * source variable or discard its dependency.
         * NB: ArgumentRule::fitArgument must keep this object separate from the source variable; replacing the
         * source node with its own conversion child would create a cycle.
         * NB: we do *not* check whether a similar conversion node already exists, as the performance cost of
         * checking seems likely to be higher than the cost of duplicating in most circumstances.
         * COSTS HAVE NOT BEEN CHECKED and in theory this could end up with many duplicated conversion nodes, so
         * this may need fixing in the future.
         *
         * \return the type-converted object
         */
        template<class rbTypeFrom, class rbTypeTo>
        RevObject* RlTypeConverter::convertTo(const ModelObject<typename rbTypeFrom::valueType>* input) 
        {
            typedef typename rbTypeFrom::valueType fromValueType;
            typedef typename rbTypeTo::valueType toValueType;
            Func__conversion<rbTypeFrom,rbTypeTo>* rlFunc = new Func__conversion<rbTypeFrom,rbTypeTo>();
            RevBayesCore::TypeConversionFunction<fromValueType,toValueType>* func = new RevBayesCore::TypeConversionFunction<fromValueType,toValueType>(input->getDagNode());
            DeterministicNode<toValueType>* newnode = new DeterministicNode<toValueType>(input->getDagNode()->getName() + "2" + rbTypeTo::getClassType(), func, rlFunc);
            return new rbTypeTo(newnode);
        }
    }
}

#endif

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
         * Auxiliary function (templated) to convert to another type
         * If the current node is constant, we just create another constant node with the correct type to replace this one, as this is faster
         * NB: changing this behaviour for constant nodes also requires changing the corresponding code in ArgumentRule::fitArgument, otherwise infinite loops are created
         * Otherwise, the node value may change, so we create a deterministic node tied to this one by a type conversion function, which will handle the updates
         * NB: we do *not* check whether a similar conversion node already exists, as the perfomance cost of checking seems likely to be higher than the cost of duplicating in most circumstances
         * COSTS HAVE NOT BEEN CHECKED and in theory this could end up with many duplicated conversion nodes, so this may need fixing in the future
         * 
         * \return the type-converted object
         */
        template<class rbTypeFrom, class rbTypeTo>
        RevObject* RlTypeConverter::convertTo(const ModelObject<typename rbTypeFrom::valueType>* input) 
        {   
            if(!input->isConstant()) {
                typedef typename rbTypeFrom::valueType fromValueType;
                typedef typename rbTypeTo::valueType toValueType;
                Func__conversion<rbTypeFrom,rbTypeTo>* rlFunc = new Func__conversion<rbTypeFrom,rbTypeTo>();
                RevBayesCore::TypeConversionFunction<fromValueType,toValueType>* func = new RevBayesCore::TypeConversionFunction<fromValueType,toValueType>(input->getDagNode());
                DeterministicNode<toValueType>* newnode = new DeterministicNode<toValueType>(input->getDagNode()->getName() + "2" + rbTypeTo::getClassType(), func, rlFunc);
                return new rbTypeTo(newnode);
            } else {
                return new rbTypeTo(input->getDagNode()->getValue());
            }
        }
    }
}

#endif

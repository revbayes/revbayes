/**
 * @file
 * This file contains the declaration of the probability density function of some distribution.
 * The distribution is given as a parameter so this file is the wrapper to call the pdf of the distribution.
 * Hence, this function can be used inside deterministic nodes.
 *
 * @brief Declaration of the probability density function.
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-06-20 22:57:09 +0200 (Wed, 20 Jun 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-08-13, version 1.0
 *
 * $Id: RandomNumberFactory.h 1643 2012-06-20 20:57:09Z hoehna $
 */

#ifndef DiscretizeBetaQuadratureFunction_H
#define DiscretizeBetaQuadratureFunction_H

#include <cstdint>

#include "TypedFunction.h"
#include "MatrixReal.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    class DiscretizeBetaQuadratureFunction : public TypedFunction< MatrixReal >{
        
    public:
        DiscretizeBetaQuadratureFunction(const TypedDagNode<double> *a, const TypedDagNode<double> *b, const TypedDagNode<std::int64_t> *nc);
        
        DiscretizeBetaQuadratureFunction*   clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value
        
    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        const TypedDagNode<double>*         alpha;
        const TypedDagNode<double>*         beta;
        const TypedDagNode<std::int64_t>*           numCats;
    };
}


#endif

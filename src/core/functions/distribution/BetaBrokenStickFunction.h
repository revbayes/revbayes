#ifndef BetaBrokenStickFunction_H
#define BetaBrokenStickFunction_H

#include "TypedFunction.h"
#include "Simplex.h"
#include "RbVector.h"
#include "TypedDagNode.h"

/**
 *
 * @brief Declaration of the discretized Beta function.
 *
 */

namespace RevBayesCore {

    class BetaBrokenStickFunction : public TypedFunction< Simplex >{

    public:
        BetaBrokenStickFunction(const TypedDagNode<double> *s, const TypedDagNode<double> *r, const TypedDagNode<long> *nc);

        BetaBrokenStickFunction*            clone(void) const;                                                  //!< Create a clon.
        void                                update(void);                                                       //!< Recompute the value

    protected:
        void                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters

    private:
        const TypedDagNode<double>*         alpha;
        const TypedDagNode<double>*         beta;
        const TypedDagNode<long>*           numCats;
    };
}


#endif

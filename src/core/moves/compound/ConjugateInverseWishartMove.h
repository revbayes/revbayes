#ifndef ConjugateInverseWishartMove_H
#define ConjugateInverseWishartMove_H

#include <map>
#include <ostream>
#include <set>
#include <string>

#include "AbstractGibbsMove.h"
#include "MatrixReal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {
    
    class ConjugateInverseWishartMove : public AbstractGibbsMove {
        
    public:
        ConjugateInverseWishartMove(StochasticNode<MatrixReal>* s, double w);                                                         //!<  constructor
        virtual                                ~ConjugateInverseWishartMove(void);                                                             //!< Destructor
        
        // Basic utility functions
        ConjugateInverseWishartMove*            clone(void) const;                                                                  //!< Clone object
        const std::string&                      getMoveName(void) const;                                                            //!< Get the name of the move for summary printing
        
    protected:
        void                                    performGibbsMove(void);                                 //!< Perform move
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                                             //!< Swap the pointers to the variable on which the move works on.
        
    private:
        
        // member variables
        StochasticNode<MatrixReal>*             sigma;
        
    };
    
}

#endif

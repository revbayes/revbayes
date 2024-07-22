/**
 * @file
 * This file contains the declaration of the CAFERateMatrixFunction class.
 * This class is derived from the function class and is used to
 * create the CAFE rate matrix.
 *
 * @brief Declaration of the CAFERateMatrixFunction.
 *
 * (c) Copyright 2014- under GPL version 3
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.2.5
 */


#ifndef CAFERateMatrixFunction_H
#define CAFERateMatrixFunction_H

#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {

    class DagNode;
    template <class valueType> class TypedDagNode;
    
    class CAFERateMatrixFunction : public TypedFunction<RateGenerator> {

    public:
        CAFERateMatrixFunction( size_t n, const TypedDagNode<double> *b, const TypedDagNode<double> *d);
        
        virtual                                     ~CAFERateMatrixFunction(void);                                      //!< Virtual destructor
        
        // public member functions
        CAFERateMatrixFunction*                     clone(void) const;                                                  //!< Create an independent clone
        void                                        update(void);
        
    protected:
        void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);    //!< Implementation of swaping parameters
        
    private:
        
        // members

        const TypedDagNode<double>*                 birth;
        const TypedDagNode<double>*                 death;

    };
    
}

#endif

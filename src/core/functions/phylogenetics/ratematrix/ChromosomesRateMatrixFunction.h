/**
 * @file
 * This file contains the declaration of the ChromosomesRateMatrixFunction class.
 * This class is derived from the function class and is used to
 * create the chromosomes rate matrix.
 *
 * @brief Declaration of the ChromosomesRateMatrixFunction.
 *
 * (c) Copyright 2014- under GPL version 3
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 */


#ifndef ChromosomesRateMatrixFunction_H
#define ChromosomesRateMatrixFunction_H

#include <cstdint>

#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    class ChromosomesRateMatrixFunction : public TypedFunction<RateGenerator> {

    public:
        ChromosomesRateMatrixFunction(const TypedDagNode<std::int64_t> *n, const TypedDagNode<double> *g, const TypedDagNode<double> *d, const TypedDagNode<double> *r, const TypedDagNode<double> *e, const TypedDagNode<double> *g_l, const TypedDagNode<double> *d_l);
        
        virtual                                     ~ChromosomesRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        ChromosomesRateMatrixFunction*              clone(void) const;                                                              //!< Create an independent clone
        void                                        update(void);
        
    protected:
        void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members

        const TypedDagNode<std::int64_t>*                    n;
        const TypedDagNode<double>*                 gamma;
        const TypedDagNode<double>*                 delta;
        const TypedDagNode<double>*                 rho;
        const TypedDagNode<double>*                 eta;
        const TypedDagNode<double>*                 gamma_l;
        const TypedDagNode<double>*                 delta_l;

    };
    
}

#endif

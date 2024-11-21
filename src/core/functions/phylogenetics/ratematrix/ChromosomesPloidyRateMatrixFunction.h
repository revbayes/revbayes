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


#ifndef ChromosomesPloidyRateMatrixFunction_H
#define ChromosomesPloidyRateMatrixFunction_H

#include <cstdint>

#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    class ChromosomesPloidyRateMatrixFunction : public TypedFunction<RateGenerator> {

    public:
        ChromosomesPloidyRateMatrixFunction(const TypedDagNode<std::int64_t> *n, const TypedDagNode<double> *g_d, const TypedDagNode<double> *d_d, const TypedDagNode<double> *r_d, const TypedDagNode<double> *e_d, const TypedDagNode<double> *g_dl, const TypedDagNode<double> *d_dl, const TypedDagNode<double> *g_p, const TypedDagNode<double> *d_p, const TypedDagNode<double> *r_p, const TypedDagNode<double> *e_p, const TypedDagNode<double> *g_pl, const TypedDagNode<double> *d_pl, const TypedDagNode<double> *b);
        
        virtual                                     ~ChromosomesPloidyRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        ChromosomesPloidyRateMatrixFunction*        clone(void) const;                                                              //!< Create an independent clone
        void                                        update(void);
        
    protected:
        void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<std::int64_t>*                    n;
        const TypedDagNode<double>*                 gamma_d;
        const TypedDagNode<double>*                 delta_d;
        const TypedDagNode<double>*                 rho_d;
        const TypedDagNode<double>*                 eta_d;
        const TypedDagNode<double>*                 gamma_dl;
        const TypedDagNode<double>*                 delta_dl;
        const TypedDagNode<double>*                 gamma_p;
        const TypedDagNode<double>*                 delta_p;
        const TypedDagNode<double>*                 rho_p;
        const TypedDagNode<double>*                 eta_p;
        const TypedDagNode<double>*                 gamma_pl;
        const TypedDagNode<double>*                 delta_pl;
        const TypedDagNode<double>*                 beta;

    };
    
}

#endif

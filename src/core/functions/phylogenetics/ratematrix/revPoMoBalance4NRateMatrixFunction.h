#ifndef revPoMoBalance4NRateMatrixFunction_H
#define revPoMoBalance4NRateMatrixFunction_H

#include "RateMatrix_revPoMoBalance4N.h"
#include "RbVector.h"
#include "Simplex.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>
#include <cstdint>

namespace RevBayesCore {
    
    /**
     * @brief HKY rate matrix function.
     *
     * This function creates the HKY rates matrix object by setting the transition-transversion parameter kappa
     * and the base frequencies. The rate matrix takes care of the setting of the actual rates and transition probabilities.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0.7, 2017-10-16
     *
     */
    class revPoMoBalance4NRateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:

        revPoMoBalance4NRateMatrixFunction( const TypedDagNode< std::int64_t > *n,const TypedDagNode< Simplex  > *p, const TypedDagNode< RbVector<double> > *r, const TypedDagNode< RbVector<double> > *s, const TypedDagNode< RbVector<double> > *b, const TypedDagNode< RbVector<std::int64_t> > *Bf  );

        virtual                                            ~revPoMoBalance4NRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        revPoMoBalance4NRateMatrixFunction*      clone(void) const;                                                              //!< Create an independent clone

        void                                                    update(void);
        
    protected:
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<std::int64_t>*                            N;
        const TypedDagNode< Simplex >*                      pi;
        const TypedDagNode< RbVector<double> >*             rho;
        const TypedDagNode< RbVector<double> >*             phi;
        const TypedDagNode< RbVector<double> >*             beta;
        const TypedDagNode< RbVector<std::int64_t> >*                B;


    };
    
}

#endif



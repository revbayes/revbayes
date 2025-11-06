#ifndef AssembleOrder2ContinuousMRFFunction_H
#define AssembleOrder2ContinuousMRFFunction_H

#include <cstddef>

#include "RbVector.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;

    /**
     * @brief Make a MRF of order 2 on a continuous vector using increments.
     *
     * This function takes theta[1] on the log or non-log scale, and the second order differences ((theta[i+1] - theta[i]) - (theta[i] - theta[i-1])).
     * It then computes the entire vector theta, optionally exponentiating the output first.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team
     * @since Version 1.0, 2014-07-04
     *
     */
    class AssembleOrder2ContinuousMRFFunction : public TypedFunction<RbVector<double> > {

    public:
        AssembleOrder2ContinuousMRFFunction(const TypedDagNode< double > *theta1, const TypedDagNode< RbVector<double> > *increments, bool log_theta1);
        virtual                                            ~AssembleOrder2ContinuousMRFFunction(void);                                                    //!< Virtual destructor

        // public member functions
        AssembleOrder2ContinuousMRFFunction*                              clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);

    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters

    private:

        // members
        const TypedDagNode< RbVector<double> >*             increments;
        const TypedDagNode< double >*                       theta1;

        size_t                                              field_size;
        bool                                                theta_1_is_log;

    };

}

#endif

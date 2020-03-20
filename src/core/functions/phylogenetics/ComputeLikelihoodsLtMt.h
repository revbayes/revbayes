#ifndef ComputeLikelihoodsLtMt_H
#define ComputeLikelihoodsLtMt_H

#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "MatrixReal.h"

#include <string>
#include <vector>

namespace RevBayesCore {
        
        // std::vector<Event>      poolTimes(  const TypedDagNode<double> *start_age,
        //                                     const TypedDagNode< RbVector<double> > *dn_time_points,
        //                                     const TypedDagNode<Tree> *timeTree );
        
        MatrixReal              ComputeLikelihoodsBackwardsLt(  const TypedDagNode<double> *start_age,
                                                                const TypedDagNode<double> *lambda,
                                                                const TypedDagNode<double> *mu,
                                                                const TypedDagNode<double> *psi,
                                                                const TypedDagNode<double> *omega,
                                                                const TypedDagNode<double> *rho,
                                                                const TypedDagNode<double> *removalPr,
                                                                const TypedDagNode<long> *maxHiddenLin,

                                                                const std::string& cond,
                                                                const std::vector<double> &dn_time_points,
                                                                bool useOrigin,
                                                                const Tree &timeTree);

        MatrixReal              ComputeLikelihoodsForwardsMt(   const TypedDagNode<double> *start_age,
                                                                const TypedDagNode<double> *lambda,
                                                                const TypedDagNode<double> *mu,
                                                                const TypedDagNode<double> *psi,
                                                                const TypedDagNode<double> *omega,
                                                                const TypedDagNode<double> *rho,
                                                                const TypedDagNode<double> *removalPr,
                                                                const TypedDagNode<long> *maxHiddenLin,

                                                                const std::string& cond,
                                                                const std::vector<double> &dn_time_points,
                                                                bool useOrigin,
                                                                const Tree &timeTree);
    };
    
#endif

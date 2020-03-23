#ifndef ComputeLikelihoodsLtMt_H
#define ComputeLikelihoodsLtMt_H

#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "MatrixReal.h"

#include <string>
#include <vector>

namespace RevBayesCore {

        MatrixReal              ComputeLikelihoodsBackwardsLt(  const TypedDagNode<double> *start_age,
                                                                const TypedDagNode<double> *lambda,
                                                                const TypedDagNode<double> *mu,
                                                                const TypedDagNode<double> *psi,
                                                                const TypedDagNode<double> *omega,
                                                                const TypedDagNode<double> *rho,
                                                                const TypedDagNode<double> *removalPr,
                                                                const TypedDagNode<long> *maxHiddenLin,

                                                                const std::string& cond,
                                                                const std::vector<double> &time_points,
                                                                bool useOrigin,
                                                                const std::vector<double> &occurrence_ages,
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
                                                                const std::vector<double> &time_points,
                                                                bool useOrigin,
                                                                const std::vector<double> &occurrence_ages,
                                                                const Tree &timeTree);
    };

#endif

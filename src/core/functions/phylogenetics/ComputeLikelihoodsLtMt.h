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

        MatrixReal              ComputeLikelihoodsBackwardsLtPiecewise(   const TypedDagNode<double> *start_age,
                                                                          const std::vector<double> &timeline,
                                                                          const std::vector<double> &lambda,
                                                                          const std::vector<double> &mu,
                                                                          const std::vector<double> &psi,
                                                                          const std::vector<double> &omega,
                                                                          const TypedDagNode<double> *rho,
                                                                          const std::vector<double> &removalPr,
                                                                          const TypedDagNode<long> *maxHiddenLin,
                                                                          const std::string& cond,
                                                                          const std::vector<double> &time_points,
                                                                          bool useOrigin,
                                                                          const std::vector<double> &occurrence_ages,
                                                                          const Tree &timeTree);

        MatrixReal              ComputeLikelihoodsForwardsMtPiecewise(  const TypedDagNode<double> *start_age,
                                                                        const std::vector<double> &timeline,
                                                                        const std::vector<double> &lambda,
                                                                        const std::vector<double> &mu,
                                                                        const std::vector<double> &psi,
                                                                        const std::vector<double> &omega,
                                                                        const TypedDagNode<double> *rho,
                                                                        const std::vector<double> &removalPr,
                                                                        const TypedDagNode<long> *maxHiddenLin,
                                                                        const std::string& cond,
                                                                        const std::vector<double> &time_points,
                                                                        bool useOrigin,
                                                                        const std::vector<double> &occurrence_ages,
                                                                        const Tree &timeTree);

        size_t                  LocateTimeSliceIndex(const double &t, const std::vector<double> &timeline);

    };

#endif

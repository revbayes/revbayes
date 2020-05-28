#ifndef ComputeLikelihoodsLtMt_H
#define ComputeLikelihoodsLtMt_H

#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "MatrixReal.h"

#include <string>
#include <vector>

struct Event {
    Event(double d, std::string s) : time(d), type(s) {};

    std::string type;
    double time;

    std::string getEventType(void){ return type; }
    double getEventTime(void){ return time; }
};

namespace RevBayesCore {

    MatrixReal              ComputeLnProbabilityDensitiesOBDP(  const TypedDagNode<double> *start_age,
                                                                const TypedDagNode<double> *lambda,
                                                                const TypedDagNode<double> *mu,
                                                                const TypedDagNode<double> *psi,
                                                                const TypedDagNode<double> *omega,
                                                                const TypedDagNode<double> *rho,
                                                                const TypedDagNode<double> *removalPr,
                                                                const TypedDagNode<long> *maxHiddenLin,
                                                                const std::string &cond,
                                                                const std::vector<double> &time_points,
                                                                bool useOrigin,
                                                                bool useMt,
                                                                bool verbose,
                                                                const std::vector<double> &occurrence_ages,
                                                                const Tree &timeTree);

MatrixReal          ComputeLnProbabilityDensitiesOBDPPiecewise( const TypedDagNode<double> *start_age,
                                                                const std::vector<double> &timeline,
                                                                const std::vector<double> &lambda,
                                                                const std::vector<double> &mu,
                                                                const std::vector<double> &psi,
                                                                const std::vector<double> &omega,
                                                                const TypedDagNode<double> *rho,
                                                                const std::vector<double> &removalPr,
                                                                const TypedDagNode<long> *maxHiddenLin,
                                                                const std::string &cond,
                                                                const std::vector<double> &time_points,
                                                                bool useOrigin,
                                                                bool useMt,
                                                                bool verbose,
                                                                const std::vector<double> &occurrence_ages,
                                                                const Tree &timeTree);

    double                  ComputeLnLikelihoodOBDP(    const TypedDagNode<double> *start_age,
                                                        const TypedDagNode<double> *lambda,
                                                        const TypedDagNode<double> *mu,
                                                        const TypedDagNode<double> *psi,
                                                        const TypedDagNode<double> *omega,
                                                        const TypedDagNode<double> *rho,
                                                        const TypedDagNode<double> *removalPr,
                                                        const TypedDagNode<long> *maxHiddenLin,
                                                        const std::string &cond,
                                                        bool useOrigin,
                                                        bool useMt,
                                                        bool verbose,
                                                        const std::vector<double> &occurrence_ages,
                                                        const Tree &timeTree);

double               ComputeLnLikelihoodOBDPPiecewise(    const TypedDagNode<double> *start_age,
                                                          const std::vector<double> &timeline,
                                                          const std::vector<double> &lambda,
                                                          const std::vector<double> &mu,
                                                          const std::vector<double> &psi,
                                                          const std::vector<double> &omega,
                                                          const TypedDagNode<double> *rho,
                                                          const std::vector<double> &removalPr,
                                                          const TypedDagNode<long> *maxHiddenLin,
                                                          const std::string &cond,
                                                          bool useOrigin,
                                                          bool useMt,
                                                          bool verbose,
                                                          const std::vector<double> &occurrence_ages,
                                                          const Tree &timeTree);

    std::vector<Event>      PoolEvents( const TypedDagNode<double> *start_age,
                                        const std::vector<double> &time_points,
                                        bool verbose,
                                        const std::vector<double> &occurrence_ages,
                                        const Tree &timeTree);

    std::vector<Event>      PoolEventsPiecewise(  const TypedDagNode<double> *start_age,
                                                  const std::vector<double> &time_points,
                                                  const std::vector<double> &occurrence_ages,
                                                  const Tree &timeTree,
                                                  const std::vector<double> &timeline);

    MatrixReal              ForwardsTraversalMt(    const TypedDagNode<double> *start_age,
                                                    const TypedDagNode<double> *lambda,
                                                    const TypedDagNode<double> *mu,
                                                    const TypedDagNode<double> *psi,
                                                    const TypedDagNode<double> *omega,
                                                    const TypedDagNode<double> *rho,
                                                    const TypedDagNode<double> *removalPr,
                                                    const TypedDagNode<long> *maxHiddenLin,
                                                    const std::string &cond,
                                                    const std::vector<double> &time_points,
                                                    bool useOrigin,
                                                    bool returnLogLikelihood,
                                                    bool verbose,
                                                    const std::vector<double> &occurrence_ages,
                                                    const Tree &timeTree);

    MatrixReal              BackwardsTraversalLt(   const TypedDagNode<double> *start_age,
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
                                                    bool verbose,
                                                    const std::vector<double> &occurrence_ages,
                                                    const Tree &timeTree);


    MatrixReal              ForwardsTraversalMtPiecewise(  const TypedDagNode<double> *start_age,
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
                                                                    bool returnLogLikelihood,
                                                                    bool verbose,
                                                                    const std::vector<double> &occurrence_ages,
                                                                    const Tree &timeTree);

    MatrixReal              BackwardsTraversalLtPiecewise(   const TypedDagNode<double> *start_age,
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
                                                                      bool verbose,
                                                                      const std::vector<double> &occurrence_ages,
                                                                      const Tree &timeTree);

    size_t                  LocateTimeSliceIndex(const double &t, const std::vector<double> &timeline);

    /////leftovers
    MatrixReal ComputeLikelihoodsForwardsMtPiecewise(  const TypedDagNode<double> *start_age,
                                                        const std::vector< double > &timeline,
                                                        const std::vector< double > &lambda,
                                                        const std::vector< double > &mu,
                                                        const std::vector< double > &psi,
                                                        const std::vector< double > &omega,
                                                        const TypedDagNode<double> *rho,
                                                        const std::vector< double > &removalPr,
                                                        const TypedDagNode<long> *maxHiddenLin,
                                                        const std::string& cond,
                                                        const std::vector<double> &time_points,
                                                        bool useOrigin,
                                                        const std::vector<double> &occurrence_ages,
                                                        const Tree &timeTree  );

    MatrixReal ComputeLikelihoodsBackwardsLtPiecewise(  const TypedDagNode<double> *start_age,
                                                        const std::vector< double > &timeline,
                                                        const std::vector< double > &lambda,
                                                        const std::vector< double > &mu,
                                                        const std::vector< double > &psi,
                                                        const std::vector< double > &omega,
                                                        const TypedDagNode<double> *rho,
                                                        const std::vector< double > &removalPr,
                                                        const TypedDagNode<long> *maxHiddenLin,
                                                        const std::string& cond,
                                                        const std::vector<double> &time_points,
                                                        bool useOrigin,
                                                        const std::vector<double> &occurrence_ages,
                                                        const Tree &timeTree  );

    };

#endif

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


MatrixReal          ComputeLnProbabilityDensitiesOBDP(          const TypedDagNode<double> *start_age,
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


double               ComputeLnLikelihoodOBDP(    const TypedDagNode<double> *start_age,
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


    std::vector<Event>      PoolEvents         (  const TypedDagNode<double> *start_age,
                                                  const std::vector<double> &time_points,
                                                  const std::vector<double> &occurrence_ages,
                                                  bool verbose,
                                                  const Tree &timeTree,
                                                  const std::vector<double> &timeline);




    MatrixReal              ForwardsTraversalMt(  const TypedDagNode<double> *start_age,
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

    MatrixReal              BackwardsTraversalLt(   const TypedDagNode<double> *start_age,
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


    double              likelihoodWithAllSamplesRemoved(   const TypedDagNode<double> *start_age,
                                                                      const std::vector<double> &timeline,
                                                                      const std::vector<double> &lambda,
                                                                      const std::vector<double> &mu,
                                                                      const std::vector<double> &psi,
                                                                      const std::vector<double> &omega,
                                                                      const TypedDagNode<double> *rho,
                                                                      const std::vector<double> &removalPr,
                                                                      const std::string& cond,
                                                                      const std::vector<double> &time_points,
                                                                      bool useOrigin,
                                                                      bool verbose,
                                                                      const std::vector<double> &occurrence_ages,
                                                                      const Tree &timeTree);
    
    double GetQ(  const double t, const double beta, const double rhoc, const double mu, const double psi, const double omega );
    
    double GetDerivativeQ(  const double t, const double beta, const double rhoc, const double mu, const double psi, const double omega, const unsigned n );

    unsigned nChoosek( unsigned n, unsigned k );

    unsigned factorial( unsigned n );

    double GetMultiDerivativeRecQ(  const double t, const double beta, const double rhoc, const double mu, const double psi, const double omega, const unsigned n, const unsigned NumObservedLineages );

    double GetP0(  const double t, const double beta, const double rhoc, const double mu, const double psi, const double omega);

    double GetDerP0(  const double t, const double beta, const double rhoc, const double mu, const double psi, const double omega, const unsigned n);

    std::vector<double> TransformDerivativeContrVec( const double t, const double beta, const double rhoc, const double mu, const double psi, const double omega, const unsigned NumObservedLineages, const std::vector<double> v );

    MatrixReal IncompleteBellPolynomial(unsigned N, unsigned K, std::vector<double> Vector);



    };
#endif

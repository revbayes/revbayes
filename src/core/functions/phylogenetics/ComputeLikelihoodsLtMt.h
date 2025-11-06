/**
 * @file
 * This file contains the functions to compute Mt and Lt matrices of log-probability densities
 * of the tree and occurrences under the Occurrence Birth-Death Process with a breadth-first
 * forward or backward traversal algorithm, and to deduce the log phylogenetic likelihood.
 *
 * @brief Compute joint log-probability density of the tree and occurrences under the Occurrence Birth-Death Process.
 *
 * @author Antoine Zwaans & Jérémy Andréoletti
 * @version 1.0
 * @since 2020-03, version 1.0
 *
 */

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


    std::vector<Event>  PoolEvents                        (   const double &start_age,
                                                              const std::vector<double> &time_points,
                                                              const std::vector<double> &occurrence_ages,
                                                              bool verbose,
                                                              const Tree &timeTree,
                                                              const std::vector<double> &timeline);           //!< Construct the vector containig all branching and sampling times + time points at which we want to compute the density.


    MatrixReal          ComputeLnProbabilityDensitiesOBDP (   const double &start_age,
                                                              const std::vector<double> &timeline,
                                                              const std::vector<double> &lambda,
                                                              const std::vector<double> &mu,
                                                              const std::vector<double> &psi,
                                                              const std::vector<double> &omega,
                                                              const TypedDagNode<double> *rho,
                                                              const std::vector<double> &removalPr,
                                                              const TypedDagNode<std::int64_t> *maxHiddenLin,
                                                              const std::string &cond,
                                                              const std::vector<double> &time_points,
                                                              bool useMt,
                                                              bool verbose,
                                                              const std::vector<double> &occurrence_ages,
                                                              const Tree &timeTree);                          //!< Compute the joint log-probability density of the observations made up to any time t (direction depending on the algorithm).


    double              ComputeLnLikelihoodOBDP           (   const double &start_age,
                                                              const std::vector<double> &timeline,
                                                              const std::vector<double> &lambda,
                                                              const std::vector<double> &mu,
                                                              const std::vector<double> &psi,
                                                              const std::vector<double> &omega,
                                                              const TypedDagNode<double> *rho,
                                                              const std::vector<double> &removalPr,
                                                              const TypedDagNode<std::int64_t> *maxHiddenLin,
                                                              const std::string &cond,
                                                              bool useMt,
                                                              bool verbose,
                                                              const std::vector<double> &occurrence_ages,
                                                              const Tree &timeTree);                          //!< Compute the joint log-probability density of the tree and occurrences.


    MatrixReal          ForwardsTraversalMt               (   const double &start_age,
                                                              const std::vector<double> &timeline,
                                                              const std::vector<double> &lambda,
                                                              const std::vector<double> &mu,
                                                              const std::vector<double> &psi,
                                                              const std::vector<double> &omega,
                                                              const TypedDagNode<double> *rho,
                                                              const std::vector<double> &removalPr,
                                                              const TypedDagNode<std::int64_t> *maxHiddenLin,
                                                              const std::string& cond,
                                                              const std::vector<double> &time_points,
                                                              bool returnLogLikelihood,
                                                              bool verbose,
                                                              const std::vector<double> &occurrence_ages,
                                                              const Tree &timeTree);                          //!< Compute the joint probability density of the observations made up to any time t and the population size at that time, as time decreases towards present (breadth-first forward traversal algorithm), with piecewise constant rates.


    MatrixReal          BackwardsTraversalLt              (   const double &start_age,
                                                              const std::vector<double> &timeline,
                                                              const std::vector<double> &lambda,
                                                              const std::vector<double> &mu,
                                                              const std::vector<double> &psi,
                                                              const std::vector<double> &omega,
                                                              const TypedDagNode<double> *rho,
                                                              const std::vector<double> &removalPr,
                                                              const TypedDagNode<std::int64_t> *maxHiddenLin,
                                                              const std::string& cond,
                                                              const std::vector<double> &time_points,
                                                              bool verbose,
                                                              const std::vector<double> &occurrence_ages,
                                                              const Tree &timeTree);                          //!< Compute the joint log probability density of observations made down any time t, conditioned on the population size at that time, as time increases towards the past (breadth-first backward traversal algorithm), with piecewise constant rates.


    double              likelihoodWithAllSamplesRemoved   (   const double &start_age,
                                                              const std::vector<double> &timeline,
                                                              const std::vector<double> &lambda,
                                                              const std::vector<double> &mu,
                                                              const std::vector<double> &psi,
                                                              const std::vector<double> &omega,
                                                              const TypedDagNode<double> *rho,
                                                              const std::vector<double> &removalPr,
                                                              const std::string& cond,
                                                              const std::vector<double> &time_points,
                                                              bool verbose,
                                                              const std::vector<double> &occurrence_ages,
                                                              const Tree &timeTree);                          //!< Compute the joint log probability density of observations made down any time t, conditioned on the population size at that time, as time increases towards the past (breadth-first backward traversal algorithm), with piecewise constant rates. This function is borrowed from Ankit Gupta, see his paper: "The probability distribution of the reconstructed phylogenetic tree with occurrence data".

    double     GetQ( const double t, const double beta, const double rhoc, const double mu, const double psi, const double omega );

    double     GetDerivativeQ( const double t, const double beta, const double rhoc, const double mu, const double psi, const double omega, const unsigned n );

    unsigned   nChoosek( unsigned n, unsigned k );

    unsigned   factorial( unsigned n );

    double     GetMultiDerivativeRecQ( const double t, const double beta, const double rhoc, const double mu, const double psi, const double omega, const unsigned n, const unsigned NumObservedLineages );

    double     GetP0( const double t, const double beta, const double rhoc, const double mu, const double psi, const double omega );

    double     GetDerP0( const double t, const double beta, const double rhoc, const double mu, const double psi, const double omega, const unsigned n );

    void       TransformDerivativeContrVec( const double t, const double beta, const double rhoc, const double mu, const double psi, const double omega, const unsigned NumObservedLineages, std::vector<double>& v );

    MatrixReal IncompleteBellPolynomial( unsigned N, unsigned K, const std::vector<double> Vector );

    std::vector<double> GetFunctionUandP( const double &start_age, const std::vector<double> &timeline, const std::vector<double> &lambda, const std::vector<double> &mu,
                                          const std::vector<double> &psi, const std::vector<double> &omega, const TypedDagNode<double> *rho, const std::vector<double> &removalPr );


    };
#endif

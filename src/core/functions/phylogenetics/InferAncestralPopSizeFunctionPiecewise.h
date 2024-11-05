#ifndef InferAncestralPopSizeFunctionPiecewise_h
#define InferAncestralPopSizeFunctionPiecewise_h

#include "RlTypedFunction.h"
#include "MatrixReal.h"

#include <string>
#include <vector>

namespace RevBayesCore {
    class DagNode;
    class Tree;
    class MatrixReal;
    class Clade;
    class Taxon;

    template <class valueType> class RbVector;
    template <class valueType> class TypedDagNode;

    /**
    * @file
    * This file contains a function to estimate ancestral population size given a tree, occurrences, and piecewise-constant phylodynamic
    * parameters under the Occurrence Birth-Death Process, as a density matrix of the number of hidden lineages through time, using the
    * algorithm introduced in Manceau & al. 2020 (http://dx.doi.org/10.1101/755561).
    *
    * @brief Computes the probability density of the number of hidden lineages through time, under the piecewise-constant Occurrence Birth-Death Process)
    *
    * @return A density matrix of the number of hidden lineages through time.
    *
    * @author Antoine Zwaans, Jérémy Andréoletti, Rachel Warnock & Marc Manceau
    * @since 2020-03, version 1.0
    */

    class InferAncestralPopSizeFunctionPiecewise : public TypedFunction<MatrixReal> {

    public:
    InferAncestralPopSizeFunctionPiecewise(             const TypedDagNode<double> *sa,
                                                        const DagNode *inspeciation,
                                                        const DagNode *inextinction,
                                                        const DagNode *inserialsampling,
                                                        const DagNode *inoccurrence,
                                                        const DagNode *ineventsampling,
                                                        const DagNode *intreatment,
                                                        const TypedDagNode<std::int64_t> *n,
                                                        const std::string& cdt,
                                                        const TypedDagNode< RevBayesCore::RbVector<double> > *O,
                                                        const std::vector<double> &tau,
                                                        TypedDagNode<Tree> *tr,
                                                        const TypedDagNode< RbVector<double> > *ht); //!< Compute the density matrix of the number of hidden lineages through time.

    virtual                                             ~InferAncestralPopSizeFunctionPiecewise(void);

        // public member functions
        InferAncestralPopSizeFunctionPiecewise*         clone(void) const;                           //!< Create an independent clone
        void                                            update(void);                                //!< Update the value of the function
        void                                            updateVectorParameters(void) const;

    protected:
          // Parameter management functions
        void                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);     //!< Swap a parameter

    private:
        // Members
        const TypedDagNode<double >*                    homogeneous_lambda;                          //!< The homogeneous birth rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_lambda;                        //!< The heterogeneous birth rates.
        const TypedDagNode<double >*                    homogeneous_mu;                              //!< The homogeneous death rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_mu;                            //!< The heterogeneous death rates.
        const TypedDagNode<double >*                    homogeneous_r;                               //!< The homogeneous conditional probability of death upon treatment.
        const TypedDagNode<RbVector<double> >*          heterogeneous_r;                             //!< The heterogeneous conditional probability of death upon treatment.
        const TypedDagNode<double >*                    homogeneous_o;                               //!< The homogeneous conditional probability of death upon treatment.
        const TypedDagNode<RbVector<double> >*          heterogeneous_o;                             //!< The heterogeneous conditional probability of death upon treatment.
        const TypedDagNode<double >*                    homogeneous_psi;                             //!< The homogeneous sampling rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_psi;                           //!< The heterogeneous sampling rates.
        const TypedDagNode<double >*                    homogeneous_rho;                             //!< The probability of sampling a tip at the present.
        const TypedDagNode<RbVector<double> >*          interval_times;


        const TypedDagNode< double > *                  start_age;                                   //!< Start age of the process.
        const TypedDagNode< std::int64_t > *                    maxHiddenLin;                                //!< The maximal number of hidden lineages.
        const std::string&                              cond;                                        //!< Condition of the process ("time" or "survival")
        const std::vector<double>                       time_points;                                 //!< Times at which density is computed
        const TypedDagNode< RbVector<double> > *        occurrences;                                 //!< Occurrence ages of incomplete fossils
        const TypedDagNode< Tree > *                    timeTree;                                    //!< Tree for which ancestral pop. size has to be computed.

        mutable std::vector<double>                     lambda;                                      //!< The speciation rate.
        mutable std::vector<double>                     mu;                                          //!< The extinction rate.
        mutable std::vector<double>                     psi;                                         //!< The fossil sampling rate.
        mutable std::vector<double>                     omega;                                       //!< The occurrence sampling rate.
        mutable std::vector<double>                     r;                                           //!< The sampling probability of extant taxa.
        mutable std::vector<double>                     psi_event;
        mutable std::vector<double>                     timeline;                                    //!< The times of the instantaneous events and rate shifts.


        };
}

#endif

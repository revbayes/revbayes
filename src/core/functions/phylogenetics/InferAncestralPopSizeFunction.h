#ifndef InferAncestralPopSizeFunction_h
#define InferAncestralPopSizeFunction_h

#include "RlTypedFunction.h"
#include "MatrixReal.h"

#include <string>
#include <cstdint>
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
    * This file contains a function to estimate ancestral population size given a tree, occurrences, and constant phylodynamic
    * parameters under the Occurrence Birth-Death Process, as a density matrix of the number of hidden lineages through time,
    * using the algorithm introduced in Manceau & al. 2020 (http://dx.doi.org/10.1101/755561).
    *
    * @brief Computes the probability density of the number of hidden lineages through time, under the Occurrence Birth-Death Process)
    *
    * @return A density matrix of the number of hidden lineages through time.
    *
    * @author Antoine Zwaans, Jérémy Andréoletti, Rachel Warnock & Marc Manceau
    * @version 1.0
    * @since 2020-03, version 1.0
    *
    */

    class InferAncestralPopSizeFunction : public TypedFunction<MatrixReal> {

    public:
        InferAncestralPopSizeFunction(                      const TypedDagNode<double> *sa,
                                                            const TypedDagNode<double> *inspeciation,
                                                            const TypedDagNode<double> *inextinction,
                                                            const TypedDagNode<double> *inserialsampling,
                                                            const TypedDagNode<double> *inoccurrence,
                                                            const TypedDagNode<double> *ineventsampling,
                                                            const TypedDagNode<double> *intreatment,
                                                            const TypedDagNode<std::int64_t>   *n,

                                                            const std::string& cdt,
                                                            const TypedDagNode< RevBayesCore::RbVector<double> > *O,
                                                            const std::vector<double> &tau,
                                                            bool vb,
                                                            TypedDagNode<Tree> *tr);               //!< Compute the density matrix of the number of hidden lineages through time

    virtual                                                 ~InferAncestralPopSizeFunction(void);

      // public member functions
      InferAncestralPopSizeFunction*                        clone(void) const;                     //!< Create an independent clone
      void                                                  update(void);                          //!< Update the value of the function


    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);     //!< Swap a parameter

    private:
        // Members
        const TypedDagNode< double > *                      start_age;                             //!< Start age of the process.
        const TypedDagNode< double > *                      lambda;                                //!< The speciation rate.
        const TypedDagNode< double > *                      mu;                                    //!< The extinction rate.
        const TypedDagNode< double > *                      psi;                                   //!< The fossil sampling rate.
        const TypedDagNode< double > *                      omega;                                 //!< The occurrence sampling rate.
        const TypedDagNode< double > *                      rho;                                   //!< The sampling probability of extant taxa.
        const TypedDagNode< double > *                      removalPr;                             //!< The removal probability after sampling.
        const TypedDagNode< std::int64_t > *                        maxHiddenLin;                          //!< The maximal number of hidden lineages.
        const std::string&                                  cond;                                  //!< Condition of the process ("survival" or "survival2")
        const std::vector<double>                           time_points;                           //!< Times at which density is computed
        const TypedDagNode< RbVector<double> > *            occurrences;                           //!< Occurrence ages of fossils not included in the tree
        const bool                                          verbose;                               //!< Display warnings and information messages.
        const TypedDagNode< Tree > *                        timeTree;                              //!< Facultative initial tree

      };
}

#endif

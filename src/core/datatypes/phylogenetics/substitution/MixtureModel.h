#ifndef MIXTURE_MODEL_H
#define MIXTURE_MODEL_H

#include <stddef.h>
#include <vector>

//#include "MatrixReal.h"
#include "Cloneable.h"
#include "TransitionProbabilityMatrix.h"

namespace RevBayesCore {

    class Tree;

    /**
     * @brief Abstract mixture model class.
     *
     * Derived models that depend on the branch lengths will additionally need to
     * incorporate a reference to the tree.
     *
     * @copyright Copyright 2022-
     * @author Benjamin Redelings
     * @since 2022-09-22, version 1.0
     */
    class MixtureModel: public Cloneable
    {
        int n_mixture_components = 0;
        int n_states = 0;

    public:
        virtual                          ~MixtureModel() = default;                                                  //!< Destructor
        virtual MixtureModel*            clone() const = 0;

        int                              getNumberOfComponents() const;
        int                              getNumberOfStates() const;
        virtual void                     calculateTransitionProbabilities(const Tree& t,
                                                                             int node,
                                                                             int mixture_component,
                                                                             TransitionProbabilityMatrix&) const = 0;

        TransitionProbabilityMatrix      calculateTransitionProbabilities(const Tree& t, int node, int mixture_component) const;
        std::vector<TransitionProbabilityMatrix> calculateTransitionProbabilities(const Tree& t, int node) const;

        virtual std::vector<double>      getRootFrequencies(int mixture_component) const = 0;
        virtual std::vector<double>      componentProbs() const = 0;

        // MatrixReal                          weightedFrequencyMatrix() const;
        // virtual std::vector<int>            get_emitted_letters() const;                                          //!<Find out what alphet letter each state emits, for markov modulated models.

    protected:

        // prevent instantiation
        MixtureModel(int m, int n);

        // virtual void                        update(void) = 0;                                                         //!< Update the rate entries of the matrix (is needed if stationarity freqs or similar have changed)
        // virtual bool simulateStochasticMapping(int node, std::vector<size_t>& transition_states, std::vector<double>& transition_times);
        // bool checkTimeReversibility(double tolerance);
        // virtual void multiplyMatrices(TransitionProbabilityMatrix& p,  TransitionProbabilityMatrix& q,  TransitionProbabilityMatrix& r) const;
        // bool                                needs_update;
    };

    // We need this for TypedDagNode<MixtureModel> for some reason...
    std::ostream&                                       operator<<(std::ostream& o, const MixtureModel& x);
}

#endif

#ifndef SUBSTITUTION_MIXTURE_MODEL_H
#define SUBSTITUTION_MIXTURE_MODEL_H

#include <stddef.h>
#include <vector>
#include <tuple>
#include <optional>

//#include "MatrixReal.h"
#include "MemberObject.h"  // For member functions.
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

    class SubstitutionMixtureModel: public Cloneable, public MemberObject< RbVector<RbVector<RbVector<double>>> >
    {
        int n_mixture_components = 0;
        int n_states = 0;

    public:
        virtual                          ~SubstitutionMixtureModel() = default;                                                  //!< Destructor
        virtual SubstitutionMixtureModel*            clone() const = 0;

        int                              getNumberOfComponents() const;
        int                              getNumberOfStates() const;

        virtual TransitionProbabilityMatrix
                                         calculateTransitionProbabilities(const Tree& t, int node, int mixture_component, double rate) const = 0;
        std::vector<TransitionProbabilityMatrix>
                                         calculateTransitionProbabilities(const Tree& t, int node, double rate) const;

        virtual bool                     simulateStochasticMapping(const Tree&, int node, int mixture_component, int rate, std::vector<size_t>&, std::vector<double>&) = 0;

        virtual std::optional<double>    rate() const = 0;
        virtual void                     scale(double f) = 0;
        void                             setRate(double r);

        virtual void                     executeMethod( const std::string &n, const std::vector<const DagNode*> &args, RbVector<RbVector<RbVector<double>>> &retValue) const;       //!< Execute the member-method

        virtual std::vector<double>      getRootFrequencies(int mixture_component) const = 0;
        virtual std::vector<double>      componentProbs() const = 0;

        // This is a hack to satisfy ModelVector<T>, which incorrectly assumes that these exist for all T.
        bool                             operator==(const SubstitutionMixtureModel&) const {return false;}
        bool                             operator!=(const SubstitutionMixtureModel&) const {return true;}
        bool                             operator<=(const SubstitutionMixtureModel&) const {return false;}

        // virtual std::vector<int>            get_emitted_letters() const;                                          //!<Find out what alphet letter each state emits, for markov modulated models.

    protected:

        // prevent instantiation
        SubstitutionMixtureModel(int m, int n);
    };

    // We need this for TypedDagNode<SubstitutionMixtureModel> for some reason...
    std::ostream&                                       operator<<(std::ostream& o, const SubstitutionMixtureModel& x);
}

#endif

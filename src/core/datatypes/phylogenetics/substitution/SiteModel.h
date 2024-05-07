#ifndef SITE_MODEL_H
#define SITE_MODEL_H

#include <cstddef>
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

    class SiteModel: public Cloneable, public MemberObject< RbVector<RbVector<double>> >
    {
    public:
        virtual                          ~SiteModel() = default;                                                  //!< Destructor
        virtual SiteModel*               clone() const = 0;

        virtual int                      getNumberOfStates() const = 0;

        virtual TransitionProbabilityMatrix
                                         calculateTransitionProbabilities(const Tree& t, int node, double rate) const = 0;

        virtual bool                     simulateStochasticMapping(const Tree&, int node, int rate, std::vector<size_t>&, std::vector<double>&) const = 0;

        virtual std::optional<double>    rate() const = 0;
        virtual void                     scale(double f) = 0;
        void                             setRate(double r);

        virtual void                     executeMethod( const std::string &n, const std::vector<const DagNode*> &args, RbVector<RbVector<double>> &retValue) const;       //!< Execute the member-method

        virtual std::vector<double>      getRootFrequencies() const = 0;

        // This is a hack to satisfy ModelVector<T>, which incorrectly assumes that these exist for all T.
        bool                             operator==(const SiteModel&) const {return false;}
        bool                             operator!=(const SiteModel&) const {return true;}
        bool                             operator<=(const SiteModel&) const {return false;}

        // virtual std::vector<int>            get_emitted_letters() const;                                          //!<Find out what alphet letter each state emits, for markov modulated models.
    };

    // We need this for TypedDagNode<SiteModel> for some reason...
    std::ostream&                                       operator<<(std::ostream& o, const SiteModel& x);
}

#endif

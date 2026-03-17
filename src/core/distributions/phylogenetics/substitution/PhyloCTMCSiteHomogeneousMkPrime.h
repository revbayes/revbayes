#ifndef PhyloCTMCSiteHomogeneousMkPrime_H
#define PhyloCTMCSiteHomogeneousMkPrime_H

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include <boost/math/special_functions/zeta.hpp>

#include "ConcreteTimeReversibleRateMatrix.h"
#include "ConstantNode.h"
#include "PhyloCTMCSiteHomogeneousConditional.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RateGenerator.h"
#include "StandardState.h"

namespace RevBayesCore {

class DagNode;
class Tree;
template <class valueType> class TypedDagNode;

class PhyloCTMCSiteHomogeneousMkPrime : public PhyloCTMCSiteHomogeneousConditional<StandardState> {
public:
    PhyloCTMCSiteHomogeneousMkPrime(const TypedDagNode<Tree>* t,
                                    size_t kObs,
                                    const TypedDagNode<std::int64_t>* u,
                                    const std::string& priorFamily,
                                    const std::vector<double>& priorParams,
                                    size_t kMax = 18,
                                    AscertainmentBias::Coding coding = AscertainmentBias::VARIABLE);
    PhyloCTMCSiteHomogeneousMkPrime(const PhyloCTMCSiteHomogeneousMkPrime& other);
    ~PhyloCTMCSiteHomogeneousMkPrime(void) override;

    PhyloCTMCSiteHomogeneousMkPrime* clone(void) const override;
    double computeLnProbability(void) override;

protected:
    void swapParameterInternal(const DagNode* oldP, const DagNode* newP) override;

private:
    RateGenerator* getOrCreateRateMatrix(size_t k);
    void updateRateMatrixIfNeeded(void);
    double computeLogRelabellingCorrection(size_t k, size_t kObs) const;
    double computeLogPrior(size_t u) const;

    size_t kObs_;
    size_t kMax_;
    size_t currentK_;
    std::string priorFamily_;
    std::vector<double> priorParams_;
    const TypedDagNode<std::int64_t>* u_;

    ConstantNode<RateGenerator>* activeRateMatrixNode_;
    std::vector<ConstantNode<RateGenerator>*> ownedRateMatrixNodes_;
    std::map<size_t, ConstantNode<RateGenerator>*> ownedNodeByK_;
    std::map<size_t, RateGenerator*> rateMatrixCache_;
    std::vector<double> logFactorial_;
};

}

inline RevBayesCore::PhyloCTMCSiteHomogeneousMkPrime::PhyloCTMCSiteHomogeneousMkPrime(
    const TypedDagNode<Tree>* t,
    size_t kObs,
    const TypedDagNode<std::int64_t>* u,
    const std::string& priorFamily,
    const std::vector<double>& priorParams,
    size_t kMax,
    AscertainmentBias::Coding coding)
    : PhyloCTMCSiteHomogeneousConditional<StandardState>(
          t,
          kMax,
          true,
          0,
          false,
          coding),
      kObs_(kObs),
      kMax_(kMax),
      currentK_(u->getValue() < 0 ? 0 : kObs + static_cast<size_t>(u->getValue())),
      priorFamily_(priorFamily),
      priorParams_(priorParams),
            u_(u),
            activeRateMatrixNode_(NULL) {
    this->addParameter(u_);

        logFactorial_.resize(kMax_ + 1, 0.0);
        for (size_t i = 1; i <= kMax_; ++i) {
                logFactorial_[i] = logFactorial_[i - 1] + std::log(static_cast<double>(i));
        }

        updateRateMatrixIfNeeded();
}

inline RevBayesCore::PhyloCTMCSiteHomogeneousMkPrime::PhyloCTMCSiteHomogeneousMkPrime(
    const PhyloCTMCSiteHomogeneousMkPrime& other)
    : PhyloCTMCSiteHomogeneousConditional<StandardState>(other),
      kObs_(other.kObs_),
      kMax_(other.kMax_),
      currentK_(other.currentK_),
      priorFamily_(other.priorFamily_),
      priorParams_(other.priorParams_),
      u_(other.u_),
    activeRateMatrixNode_(NULL),
    logFactorial_(other.logFactorial_) {
    for (std::map<size_t, RateGenerator*>::const_iterator it = other.rateMatrixCache_.begin(); it != other.rateMatrixCache_.end(); ++it) {
        rateMatrixCache_[it->first] = it->second->clone();
    }

    updateRateMatrixIfNeeded();
}

inline RevBayesCore::PhyloCTMCSiteHomogeneousMkPrime::~PhyloCTMCSiteHomogeneousMkPrime(void) {
    for (std::vector<ConstantNode<RateGenerator>*>::iterator it = ownedRateMatrixNodes_.begin(); it != ownedRateMatrixNodes_.end(); ++it) {
        delete *it;
    }

    for (std::map<size_t, RateGenerator*>::iterator it = rateMatrixCache_.begin(); it != rateMatrixCache_.end(); ++it) {
        delete it->second;
    }
}

inline RevBayesCore::PhyloCTMCSiteHomogeneousMkPrime*
RevBayesCore::PhyloCTMCSiteHomogeneousMkPrime::clone(void) const {
    return new PhyloCTMCSiteHomogeneousMkPrime(*this);
}

inline double RevBayesCore::PhyloCTMCSiteHomogeneousMkPrime::computeLnProbability(void) {
    const std::int64_t uValue = u_->getValue();
    if (uValue < 0) {
        return RbConstants::Double::neginf;
    }

    updateRateMatrixIfNeeded();

    if (currentK_ < kObs_ || currentK_ > kMax_) {
        return RbConstants::Double::neginf;
    }

    const double lnL = PhyloCTMCSiteHomogeneousConditional<StandardState>::computeLnProbability();

    size_t numChars = 0;
    for (std::vector<size_t>::const_iterator it = this->pattern_counts.begin(); it != this->pattern_counts.end(); ++it) {
        numChars += *it;
    }

    const size_t u = static_cast<size_t>(uValue);
    const double relabelTerm = static_cast<double>(numChars) * computeLogRelabellingCorrection(currentK_, kObs_);

    // Include p(u) only when u is an unclamped stochastic node.
    double priorTerm = 0.0;
    if (u_->getDagNodeType() == DagNode::STOCHASTIC && !u_->isClamped()) {
        priorTerm = computeLogPrior(u);
    }

    return lnL + relabelTerm + priorTerm;
}

inline void RevBayesCore::PhyloCTMCSiteHomogeneousMkPrime::swapParameterInternal(const DagNode* oldP,
                                                                                  const DagNode* newP) {
    if (oldP == u_) {
        u_ = static_cast<const TypedDagNode<std::int64_t>*>(newP);
    } else if (oldP->getName() == "mkPrimeRateMatrix" || newP->getName() == "mkPrimeRateMatrix") {
        // Internal mkPrime rate-matrix nodes can be replaced between cached k values.
        // Ignore explicit swap bookkeeping here; updateRateMatrixIfNeeded() maintains
        // the active pointer and dirty flags during likelihood evaluation.
        return;
    } else {
        PhyloCTMCSiteHomogeneousConditional<StandardState>::swapParameterInternal(oldP, newP);
    }
}

inline RevBayesCore::RateGenerator* RevBayesCore::PhyloCTMCSiteHomogeneousMkPrime::getOrCreateRateMatrix(size_t k) {
    std::map<size_t, RateGenerator*>::iterator cached = rateMatrixCache_.find(k);
    if (cached != rateMatrixCache_.end()) {
        return cached->second;
    }

    if (k < 2 || k > kMax_) {
        throw RbException() << "dnMkPrime: invalid k=" << k << " (expected 2.." << kMax_ << ")";
    }

    const double activeExchangeability = static_cast<double>(k) / static_cast<double>(k - 1);
    std::vector<double> exchangeability(kMax_ * (kMax_ - 1) / 2, 0.0);
    std::vector<double> stationary(kMax_, 0.0);

    for (size_t i = 0; i < k; ++i) {
        stationary[i] = 1.0 / static_cast<double>(k);
    }

    size_t idx = 0;
    for (size_t i = 0; i < kMax_; ++i) {
        for (size_t j = i + 1; j < kMax_; ++j) {
            exchangeability[idx] = (i < k && j < k) ? activeExchangeability : 0.0;
            ++idx;
        }
    }

    RateGenerator* matrix = new ConcreteTimeReversibleRateMatrix(exchangeability, stationary, boost::optional<double>());
    rateMatrixCache_[k] = matrix;
    return matrix;
}

inline void RevBayesCore::PhyloCTMCSiteHomogeneousMkPrime::updateRateMatrixIfNeeded(void) {
    const std::int64_t uValue = u_->getValue();
    if (uValue < 0) {
        currentK_ = 0;
        return;
    }

    const size_t u = static_cast<size_t>(uValue);
    const size_t proposedK = kObs_ + u;
    if (proposedK > kMax_) {
        currentK_ = proposedK;
        return;
    }

    if (activeRateMatrixNode_ != NULL && proposedK == currentK_) {
        return;
    }

    RateGenerator* matrix = getOrCreateRateMatrix(proposedK);

    // Find or create an owned ConstantNode for proposedK (one per distinct k).
    ConstantNode<RateGenerator>* newRateMatrixNode = NULL;
    std::map<size_t, ConstantNode<RateGenerator>*>::iterator it = ownedNodeByK_.find(proposedK);
    if (it != ownedNodeByK_.end()) {
        newRateMatrixNode = it->second;
    } else {
        newRateMatrixNode = new ConstantNode<RateGenerator>("mkPrimeRateMatrix", matrix->clone());
        ownedRateMatrixNodes_.push_back(newRateMatrixNode);
        ownedNodeByK_[proposedK] = newRateMatrixNode;
    }

    // Assign directly WITHOUT going through Distribution::addParameter or setRateMatrix.
    // These rate-matrix nodes are privately owned by MkPrime and must never enter the
    // Distribution::parameters vector.  If they do, StochasticNode::StochasticNode()
    // registers them as DAG parents; subsequent k-switches call removeParameter()
    // (which deletes the node when refCount reaches 0, invalidating ownedNodeByK_)
    // and desynchronise the StochasticNode parent list, causing
    // "Could not find the distribution parameter to be swapped" errors on swapParent.
    this->homogeneous_rate_matrix = newRateMatrixNode;

    // A k-change modifies Q, so both transition probabilities and partials must be recomputed.
    for (size_t i = 0; i < this->dirty_nodes.size(); ++i) {
        this->dirty_nodes[i] = true;
        this->changed_nodes[i] = true;
    }
    for (size_t i = 0; i < this->pmat_dirty_nodes.size(); ++i) {
        this->pmat_dirty_nodes[i] = true;
        this->pmat_changed_nodes[i] = true;
    }

    activeRateMatrixNode_ = newRateMatrixNode;
    currentK_ = proposedK;
}

inline double RevBayesCore::PhyloCTMCSiteHomogeneousMkPrime::computeLogRelabellingCorrection(size_t k, size_t kObs) const {
    if (k == 0 || k > kMax_ || kObs > kMax_ || kObs > k) {
        return RbConstants::Double::neginf;
    }

    return logFactorial_[kObs] -
           static_cast<double>(kObs) * std::log(static_cast<double>(kObs)) -
           logFactorial_[k] +
           logFactorial_[k - kObs] +
           static_cast<double>(kObs) * std::log(static_cast<double>(k));
}

inline double RevBayesCore::PhyloCTMCSiteHomogeneousMkPrime::computeLogPrior(size_t u) const {
    if (priorParams_.empty()) {
        throw RbException() << "dnMkPrime: prior parameter vector is empty";
    }

    if (priorFamily_ == "geometric") {
        const double p = priorParams_[0];
        if (p <= 0.0 || p >= 1.0) {
            return RbConstants::Double::neginf;
        }
        return static_cast<double>(u) * std::log(1.0 - p) + std::log(p);
    }

    if (priorFamily_ == "logseries") {
        const double c = priorParams_[0];
        if (c <= 0.0 || c >= 1.0) {
            return RbConstants::Double::neginf;
        }
        return static_cast<double>(u) * std::log(c) -
               std::log(static_cast<double>(u + 1)) -
               std::log(-std::log(1.0 - c));
    }

    if (priorFamily_ == "powerlaw") {
        const double s = priorParams_[0];
        if (s <= 1.0) {
            return RbConstants::Double::neginf;
        }
        return -s * std::log(static_cast<double>(u + 1)) - std::log(boost::math::zeta(s));
    }

    throw RbException() << "dnMkPrime: unknown prior family '" << priorFamily_ << "'";
}

#endif

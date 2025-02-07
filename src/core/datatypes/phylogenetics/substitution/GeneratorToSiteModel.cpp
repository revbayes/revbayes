#include "GeneratorToSiteModel.h"

#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstddef>

#include "MatrixReal.h"
#include "RbException.h"
#include "RbMathFunctions.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "RbMathLogic.h"

using namespace RevBayesCore;

using std::vector;
using std::unique_ptr;
using std::optional;
using std::tuple;

GeneratorToSiteModel* GeneratorToSiteModel::clone() const
{
    return new GeneratorToSiteModel(*this);
}

int GeneratorToSiteModel::getNumberOfStates() const
{
    return generator->getNumberOfStates();
}

TransitionProbabilityMatrix
GeneratorToSiteModel::calculateTransitionProbabilities(const Tree& tau, int node_index, double rate) const
{
    const TopologyNode* node = tau.getNodes()[node_index];

    if ( node->isRoot() == true )
    {
        throw RbException("GeneratorToSiteModel called updateTransitionProbabilities for the root node\n");
    }

    auto [start_age, end_age] = getStartEndAge(tau, *node);
    
    TransitionProbabilityMatrix P(getNumberOfStates());
    generator->calculateTransitionProbabilities( start_age, end_age, _scale * rate, P);
    return P;
}

bool GeneratorToSiteModel::simulateStochasticMapping(const Tree& tau, int node_index, int rate, vector<size_t>& states, vector<double>& times) const
{
    const TopologyNode* node = tau.getNodes()[node_index];

    auto [start_age, end_age] = getStartEndAge(tau, *node);

    std::vector<size_t> transition_states;
    std::vector<double> transition_times;
    return generator->simulateStochasticMapping(start_age, end_age, _scale * rate, states, times);
}

vector<double> GeneratorToSiteModel::getRootFrequencies() const
{
    return frequencies;
}

void GeneratorToSiteModel::scale(double factor)
{
    _scale *= factor;
}

optional<double> GeneratorToSiteModel::rate() const
{
    if (auto matrix = dynamic_cast<const RateMatrix*>(generator.get()))
        return _scale * matrix->averageRate();
    else
        return {};
}
        

// constructors
GeneratorToSiteModel::GeneratorToSiteModel(const RateMatrix& m, double s)
    :GeneratorToSiteModel(m, m.getStationaryFrequencies(), s)
{
}

GeneratorToSiteModel::GeneratorToSiteModel(const RateGenerator& g, const vector<double>& freqs, double s)
    :GeneratorToSiteModel(std::shared_ptr<const RateGenerator>(g.clone()), freqs, s)
{
}

GeneratorToSiteModel::GeneratorToSiteModel(const std::shared_ptr<const RateMatrix>& m, double s)
    :GeneratorToSiteModel(m, m->getStationaryFrequencies(), s)
{
}

GeneratorToSiteModel::GeneratorToSiteModel(const std::shared_ptr<const RateGenerator>& g, const vector<double>& freqs, double s)
    :generator( g ),
     frequencies(freqs),
     _scale(s)
{
}


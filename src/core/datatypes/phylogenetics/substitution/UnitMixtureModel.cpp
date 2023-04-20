#include "UnitMixtureModel.h"

#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstddef>

#include "MatrixReal.h"
#include "MixtureModel.h"
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

UnitMixtureModel* UnitMixtureModel::clone() const
{
    return new UnitMixtureModel(*this);
}

TransitionProbabilityMatrix
UnitMixtureModel::calculateTransitionProbabilities(const Tree& tau,
                                                   int node_index,
                                                   int m,
                                                   double rate) const
{
    assert(m == 0);

    const TopologyNode* node = tau.getNodes()[node_index];

    if ( node->isRoot() == true )
    {
        throw RbException("UnitMixtureModel called updateTransitionProbabilities for the root node\n");
    }

    auto [start_age, end_age] = getStartEndAge(*node);
    
    TransitionProbabilityMatrix P(getNumberOfStates());
    generator->calculateTransitionProbabilities( start_age, end_age, _scale * rate, P);
    return P;
}

bool UnitMixtureModel::simulateStochasticMapping(const Tree& tau, int node_index, int m, int rate, vector<size_t>& states, vector<double>& times)
{
    assert(m == 0);

    const TopologyNode* node = tau.getNodes()[node_index];

    auto [start_age, end_age] = getStartEndAge(*node);

    std::vector<size_t> transition_states;
    std::vector<double> transition_times;
    return generator->simulateStochasticMapping(start_age, end_age, _scale * rate, states, times);
}

vector<double> UnitMixtureModel::getRootFrequencies(int) const
{
    return frequencies;
}

vector<double> UnitMixtureModel::componentProbs() const { return {1.0};}

void UnitMixtureModel::scale(double factor)
{
    _scale *= factor;
}

optional<double> UnitMixtureModel::rate() const
{
    if (auto matrix = dynamic_cast<const RateMatrix*>(generator.get()))
        return _scale * matrix->averageRate();
    else
        return {};
}
        

// assignment operator / copy constructor
UnitMixtureModel& UnitMixtureModel::operator=(const UnitMixtureModel& m)
{
    MixtureModel::operator=(m);
    generator = unique_ptr<RateGenerator>(m.generator->clone());
    frequencies = m.frequencies;
    _scale = m._scale;
    return *this;
}

UnitMixtureModel::UnitMixtureModel(const UnitMixtureModel& m)
    :MixtureModel(m)
{
    operator=(m);
}


// move assignment operator / move  constructor
UnitMixtureModel& UnitMixtureModel::operator=(UnitMixtureModel&& m)
{
    MixtureModel::operator=(std::move(m));
    generator = std::move(m.generator);
    frequencies = m.frequencies;
    _scale = m._scale;
    return *this;
}

UnitMixtureModel::UnitMixtureModel(UnitMixtureModel&& m)
    :MixtureModel(m)
{
    operator=(std::move(m));
}


// constructors
UnitMixtureModel::UnitMixtureModel(const RateMatrix& m, double s)
    :MixtureModel(1, m.getNumberOfStates()),
     generator( unique_ptr<RateGenerator>(m.clone()) ),
     frequencies (m.getStationaryFrequencies()),
     _scale(s)
{

}

UnitMixtureModel::UnitMixtureModel(const RateMatrix& m, const vector<double>& freqs, double s)
    :MixtureModel(1, m.getNumberOfStates()),
     generator( unique_ptr<RateGenerator>(m.clone()) ),
     frequencies(freqs),
     _scale(s)
{

}

UnitMixtureModel::UnitMixtureModel(const RateGenerator& g, const vector<double>& freqs, double s)
    :MixtureModel(1, g.getNumberOfStates()),
     generator( unique_ptr<RateGenerator>(g.clone()) ),
     frequencies(freqs),
     _scale(s)
{

}


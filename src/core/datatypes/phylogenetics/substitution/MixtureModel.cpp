#include "MixtureModel.h"

#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstddef>

#include "MatrixReal.h"
#include "MixtureModel.h"
#include "RbException.h"

using namespace RevBayesCore;

using std::vector;

MixtureModel::MixtureModel(int m, int n)
    :n_mixture_components(m),
     n_states(n)
{
    if (m <= 0)
        throw RbException()<<"Mixture model: number of components is "<<m<<", but must be at least 1";
}

int MixtureModel::getNumberOfComponents() const
{
    return n_mixture_components;
}

int MixtureModel::getNumberOfStates() const
{
    return n_states;
}

TransitionProbabilityMatrix MixtureModel::calculateTransitionProbabilities(const Tree& t, int node, int mixture_component) const
{
    TransitionProbabilityMatrix P(getNumberOfStates());
    calculateTransitionProbabilities(t, node, mixture_component, P);
    return P;
}

vector<TransitionProbabilityMatrix> MixtureModel::calculateTransitionProbabilities(const Tree& t, int node) const
{
    vector<TransitionProbabilityMatrix> Ps;
    for(int c=0; c < getNumberOfComponents(); c++)
        Ps.push_back( calculateTransitionProbabilities(t, node, c) );
    return Ps;
}

namespace RevBayesCore
{

std::ostream& operator<<(std::ostream& o, const MixtureModel& m)
{
    o<<"MixtureModel with "<<m.getNumberOfComponents()<<" components and "<<m.getNumberOfStates()<<" states.";
    return o;
}

}

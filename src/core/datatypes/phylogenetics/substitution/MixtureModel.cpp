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

vector<TransitionProbabilityMatrix> MixtureModel::calculateTransitionProbabilities(const Tree& t, int node, double rate) const
{
    vector<TransitionProbabilityMatrix> Ps;
    for(int c=0; c < getNumberOfComponents(); c++)
        Ps.push_back( calculateTransitionProbabilities(t, node, c, rate) );
    return Ps;
}

void MixtureModel::setRate(double r)
{
    if (auto R = rate())
    {
        if (*R == 0 and r > 0)
            throw RbException()<<"Cannot set rate to "<<r<<" because the rate 0.";

        scale(r / *R);
    }
    else
        throw RbException()<<"Cannot set rate to "<<r<<" because the rate is not defined for this model.";
}

namespace RevBayesCore
{

std::ostream& operator<<(std::ostream& o, const MixtureModel& m)
{
    o<<"MixtureModel with "<<m.getNumberOfComponents()<<" components and "<<m.getNumberOfStates()<<" states.";
    return o;
}

}

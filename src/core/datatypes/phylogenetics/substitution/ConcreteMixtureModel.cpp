#include "ConcreteMixtureModel.h"

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

#include "UnitMixtureModel.h"
using namespace RevBayesCore;

using std::vector;
using std::unique_ptr;
using std::optional;
using std::tuple;

int sum_submodels( const vector<unique_ptr<MixtureModel>>& ms)
{
    int total = 0;
    for(auto& m: ms)
        total += m->getNumberOfComponents();
    return total;
}

vector<unique_ptr<MixtureModel>> clone_models(const vector<unique_ptr<MixtureModel>>& ms)
{
    vector<unique_ptr<MixtureModel>> submodels;
    for(auto& m: ms)
        submodels.push_back( unique_ptr<MixtureModel>(m->clone()) );
    return submodels;
}

ConcreteMixtureModel* ConcreteMixtureModel::clone() const
{
    return new ConcreteMixtureModel(*this);
}

std::pair<int,int> get_indices(const vector<unique_ptr<MixtureModel>>& ms, int i)
{
    int component = 0;
    for(auto& m: ms)
    {
        if (i < m->getNumberOfComponents())
            return {component, i};
        i -= m->getNumberOfComponents();
        component++;
    }
    std::abort();
}

TransitionProbabilityMatrix
ConcreteMixtureModel::calculateTransitionProbabilities(const Tree& tau,
                                                       int node_index,
                                                       int m,
                                                       double rate) const
{
    auto [component,index] = get_indices(sub_models, m);

    return sub_models[component]->calculateTransitionProbabilities(tau, node_index, index, rate);
}

bool ConcreteMixtureModel::simulateStochasticMapping(const Tree& tau, int node_index, int m, int rate, vector<size_t>& states, vector<double>& times)
{
    auto [component,index] = get_indices(sub_models, m);

    return sub_models[component]->simulateStochasticMapping(tau, node_index, index, rate, states, times);
}

vector<double> ConcreteMixtureModel::getRootFrequencies(int m) const
{
    auto [component, index] = get_indices(sub_models, m);

    return sub_models[component]->getRootFrequencies(index);
}

vector<double> ConcreteMixtureModel::componentProbs() const
{
    vector<double> ps;

    for(int i=0; i<fractions.size(); i++)
    {
        for(auto& sub_p: sub_models[i]->componentProbs())
            ps.push_back(fractions[i] * sub_p);
    }
    return ps;
}

void ConcreteMixtureModel::scale(double factor)
{
    for(auto& sub_model: sub_models)
        sub_model->scale(factor);
}

optional<double> ConcreteMixtureModel::rate() const
{
    double r = 0;
    for(int i=0;i<sub_models.size();i++)
    {
        if (auto sr = sub_models[i]->rate())
            r += fractions[i] * (*sr);
        else
            return {};
    }

    return r;
}

ConcreteMixtureModel::ConcreteMixtureModel(const vector<unique_ptr<MixtureModel>>& ms,
                                           const vector<double>& fs)
    :MixtureModel( sum_submodels(ms), ms[0]->getNumberOfStates() ),
     sub_models( clone_models(ms) ),
     fractions(fs)
{
    if (sub_models.size() != fractions.size())
        throw RbException()<<"ConcreteMixtureModel: got "<<sub_models.size()<<" components but "<<fractions.size()<<" weights.  These should be equal.";
}


ConcreteMixtureModel& ConcreteMixtureModel::operator=(const ConcreteMixtureModel& m)
{
    MixtureModel::operator=(m);
    sub_models = clone_models( m.sub_models );
    fractions = m.fractions;
    return *this;
}

ConcreteMixtureModel::ConcreteMixtureModel(const ConcreteMixtureModel& m)
    :MixtureModel(m),
     fractions(m.fractions)
{
    sub_models = clone_models( m.sub_models );
}

namespace RevBayesCore
{
ConcreteMixtureModel* scaled_mixture(const MixtureModel& sub_model, const std::vector<double>& fractions, const std::vector<double>& rates)
{
    vector<unique_ptr<MixtureModel>> sub_models;
    for(auto& rate: rates)
    {
        auto sub_model2 = std::unique_ptr<MixtureModel>(sub_model.clone());
        sub_model2->setRate(rate);
        sub_models.push_back( std::move(sub_model2) );
    }
    
    return new ConcreteMixtureModel(sub_models, fractions);    
}


}

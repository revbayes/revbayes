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
    int index = 0;
    for(auto& m: ms)
    {
        if (i < m->getNumberOfComponents())
            return {component, index};
        i -= m->getNumberOfComponents();
        component++;
    }
    std::abort();
}

void ConcreteMixtureModel::calculateTransitionProbabilities(const Tree& tau,
                                                        int node_index,
                                                        int m,
                                                        TransitionProbabilityMatrix& P) const
{
    const TopologyNode* node = tau.getNodes()[node_index];

    if ( node->isRoot() == true )
    {
        throw RbException("ConcreteMixtureModel called updateTransitionProbabilities for the root node\n");
    }

    double end_age = node->getAge();

    // if the tree is not a time tree, then the age will be not a number
    if ( RbMath::isFinite(end_age) == false )
    {
        // we assume by default that the end is at time 0
        end_age = 0.0;
    }
    double start_age = end_age + node->getBranchLength();
    
    auto p = get_indices(sub_models, m);
    int component = p.first;
    int index = p.second;

    sub_models[component]->calculateTransitionProbabilities( tau, node_index, index, P);
}

vector<double> ConcreteMixtureModel::getRootFrequencies(int m) const
{
    auto p = get_indices(sub_models, m);
    int component = p.first;
    int index = p.second;
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

void ConcreteMixtureModel::setRate(double r)
{
    scale( r / rate() );
}

double ConcreteMixtureModel::rate() const
{
    double r = 0;
    for(int i=0;i<getNumberOfComponents();i++)
        r += fractions[i] * sub_models[i]->rate();

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

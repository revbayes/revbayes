#include "SiteMixtureModel.h"

#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstddef>
#include <cstdint>

#include "MatrixReal.h"
#include "RbException.h"
#include "TypedDagNode.h"
#include "Tree.h"

using namespace RevBayesCore;

using std::vector;
using std::shared_ptr;

SiteMixtureModel* SiteMixtureModel::clone() const
{
    return new SiteMixtureModel(*this);
}

const SiteModel& SiteMixtureModel::getComponent(int m) const
{
    assert(m >= 0 and m < components.size());
    return *components[m];
}

const std::vector<double>& SiteMixtureModel::componentProbs() const
{
    return fractions;
}

int SiteMixtureModel::size() const
{
    assert(fractions.size() == components.size());
    return components.size();
}

int SiteMixtureModel::getNumberOfComponents() const
{
    return size();
}

int SiteMixtureModel::getNumberOfStates() const
{
    for(auto& component: components)
	assert(component->getNumberOfStates() == components[0]->getNumberOfStates());

    return getComponent(0).getNumberOfStates();
}

vector<TransitionProbabilityMatrix> SiteMixtureModel::calculateTransitionProbabilities(const Tree& t, int node, double rate) const
{
    vector<TransitionProbabilityMatrix> Ps;
    for(auto& component: components)
        Ps.push_back( component->calculateTransitionProbabilities(t, node, rate) );
    return Ps;
}

std::optional<double> SiteMixtureModel::rate() const
{
    double r = 0;
    for(int i=0;i<getNumberOfComponents();i++)
    {
        if (auto sr = components[i]->rate())
            r += fractions[i] * (*sr);
        else
            return {};
    }

    return r;
}

void SiteMixtureModel::scale(double factor)
{
    for(auto& component: components)
    {
	auto component2 = std::shared_ptr<SiteModel>(component->clone());
	component2->scale(factor);
	component = std::move(component2);
    }
}

void SiteMixtureModel::setRate(double r)
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

void SiteMixtureModel::executeMethod( const std::string &n, const std::vector<const DagNode*> &args, RbVector<RbVector<RbVector<double>>> &rv) const
{

    if (n == "getTransitionProbabilities")
    {
        // clear old values
        rv.clear();

        auto& tree = dynamic_cast<const TypedDagNode<Tree> *>( args[0] )->getValue();
        int node = dynamic_cast<const TypedDagNode<std::int64_t> *>( args[1] )->getValue() - 1; // Natural
        double rate = dynamic_cast<const TypedDagNode<double> *>( args[2] )->getValue();

        if (node < 0 or node >= tree.getNumberOfNodes())
            throw RbException()<<"SiteMixtureModel.getTransitionProbabilities(tree,node,rate): node "<<node<<" is out of range.  Tree only has "<<tree.getNumberOfNodes()<<" nodes.";

        auto Ps = calculateTransitionProbabilities( tree, node, rate );

        for(auto& P: Ps)
        {
            RbVector<RbVector<double>> matrix;
            for (size_t i = 0; i < P.getNumberOfStates(); i++)
            {
                RbVector<double> row;
                for (size_t j =0; j < P.getNumberOfStates(); j++)
                {
                    row.push_back(P[i][j]);
                }
                matrix.push_back(row);
            }
            rv.push_back(matrix);
        }
    }
}

void SiteMixtureModel::executeMethod( const std::string &n, const std::vector<const DagNode*> &args, double &rv) const
{

    if (n == "rate")
    {
	if (auto r = rate())
	    rv = *r;
	else
            throw RbException()<<"SiteMixtureModel.rate(): this model does not have a well-define rate.";
    }
}

SiteMixtureModel::SiteMixtureModel(const std::vector<std::shared_ptr<const SiteModel>>& c, const std::vector<double>& f)
    :components(c),
     fractions(f)
{
}

SiteMixtureModel::SiteMixtureModel(std::vector<std::shared_ptr<const SiteModel>>&& c, std::vector<double>&& f)
    :components(std::move(c)),
     fractions(std::move(f))
{
}

namespace RevBayesCore
{

shared_ptr<const SiteMixtureModel> scaled_mixture(const SiteMixtureModel& submodel, const std::vector<double>& fractions, const std::vector<double>& rates)
{
    const int n = fractions.size();
    assert(rates.size() == n);

    auto models = std::vector(n, std::shared_ptr<const SiteMixtureModel>(submodel.clone()));
    return mix_mixture( scale_models( models, rates), fractions);
}

shared_ptr<const SiteMixtureModel> mix_mixture(const vector<shared_ptr<const SiteMixtureModel>>& models, const std::vector<double>& fractions)
{
    const int n = fractions.size();

    vector<shared_ptr<const SiteModel>> components2;
    vector<double> fractions2;

    for(int i=0;i<n;i++)
    {
	for(int j=0;j<models[i]->size();j++)
	{
	    components2.push_back( std::shared_ptr<SiteModel>(models[i]->getComponent(j).clone()) );
	    fractions2.push_back( fractions[i] * models[i]->componentProbs()[j] );
	}
    }
    return std::make_shared<const SiteMixtureModel>( std::move(components2), std::move(fractions2) );
}

std::ostream& operator<<(std::ostream& o, const SiteMixtureModel& m)
{
    o<<"SiteMixtureModel with "<<m.getNumberOfComponents()<<" components and "<<m.getNumberOfStates()<<" states.";
    return o;
}

}

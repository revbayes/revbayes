#include "MixtureModel.h"

#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstddef>

#include "MatrixReal.h"
#include "MixtureModel.h"
#include "RbException.h"
#include "TypedDagNode.h"
#include "Tree.h"

using namespace RevBayesCore;

using std::vector;

SubstitutionMixtureModel::SubstitutionMixtureModel(int m, int n)
    :n_mixture_components(m),
     n_states(n)
{
    if (m <= 0)
        throw RbException()<<"Mixture model: number of components is "<<m<<", but must be at least 1";
}

int SubstitutionMixtureModel::getNumberOfComponents() const
{
    return n_mixture_components;
}

int SubstitutionMixtureModel::getNumberOfStates() const
{
    return n_states;
}

vector<TransitionProbabilityMatrix> SubstitutionMixtureModel::calculateTransitionProbabilities(const Tree& t, int node, double rate) const
{
    vector<TransitionProbabilityMatrix> Ps;
    for(int c=0; c < getNumberOfComponents(); c++)
        Ps.push_back( calculateTransitionProbabilities(t, node, c, rate) );
    return Ps;
}

void SubstitutionMixtureModel::setRate(double r)
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

void SubstitutionMixtureModel::executeMethod( const std::string &n, const std::vector<const DagNode*> &args, RbVector<RbVector<RbVector<double>>> &rv) const
{

    if (n == "getTransitionProbabilities")
    {
        // clear old values
        rv.clear();

        auto& tree = dynamic_cast<const TypedDagNode<Tree> *>( args[0] )->getValue();
        int node = dynamic_cast<const TypedDagNode<long> *>( args[1] )->getValue() - 1; // Natural
        double rate = dynamic_cast<const TypedDagNode<double> *>( args[2] )->getValue();

        if (node < 0 or node >= tree.getNumberOfNodes())
            throw RbException()<<"SubstitutionMixtureModel.getTransitionProbabilities(tree,node,rate): node "<<node<<" is out of range.  Tree only has "<<tree.getNumberOfNodes()<<" nodes.";

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

namespace RevBayesCore
{

std::ostream& operator<<(std::ostream& o, const SubstitutionMixtureModel& m)
{
    o<<"SubstitutionMixtureModel with "<<m.getNumberOfComponents()<<" components and "<<m.getNumberOfStates()<<" states.";
    return o;
}

}

#include "DagNodeTypeUtilities.h"
#include "SiteModel.h"

#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstddef>

#include "MatrixReal.h"
#include "RbException.h"
#include "TypedDagNode.h"
#include "Tree.h"

using namespace RevBayesCore;

using std::vector;

void SiteModel::setRate(double r)
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

void SiteModel::executeMethod( const std::string &n, const std::vector<const DagNode*> &args, RbVector<RbVector<double>> &rv) const
{
    if (n == "getTransitionProbabilities")
    {
        // clear old values
        rv.clear();

        auto& tree = RevBayesCore::getValueOfType<Tree>(args[0]);
        int node = RevBayesCore::getValueOfType<std::int64_t>(args[1]) - 1; // Natural
        double rate = RevBayesCore::getValueOfType<double>(args[2]);

        if (node < 0 or node >= tree.getNumberOfNodes())
            throw RbException()<<"SiteModel.getTransitionProbabilities(tree,node,rate): node "<<node<<" is out of range.  Tree only has "<<tree.getNumberOfNodes()<<" nodes.";

        auto P = calculateTransitionProbabilities( tree, node, rate );

	// Convert P to an RbVector<RbVector<double>>
	rv.clear();
	for (size_t i = 0; i < P.getNumberOfStates(); i++)
	{
	    RbVector<double> row;
	    for (size_t j =0; j < P.getNumberOfStates(); j++)
	    {
		row.push_back(P[i][j]);
	    }
	    rv.push_back(row);
	}
    }
}

namespace RevBayesCore
{

std::ostream& operator<<(std::ostream& o, const SiteModel& m)
{
    o<<"SiteModel with "<<m.getNumberOfStates()<<" states.";
    return o;
}

}

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

UnitMixtureModel* UnitMixtureModel::clone() const
{
    return new UnitMixtureModel(*this);
}

void UnitMixtureModel::calculateTransitionProbabilities(const Tree& tau,
                                                        int node_index,
                                                        int m,
                                                        TransitionProbabilityMatrix& P) const
{
    const TopologyNode* node = tau.getNodes()[node_index];

    if ( node->isRoot() == true )
    {
        throw RbException("UnitMixtureModel called updateTransitionProbabilities for the root node\n");
    }

    double end_age = node->getAge();

    // if the tree is not a time tree, then the age will be not a number
    if ( RbMath::isFinite(end_age) == false )
    {
        // we assume by default that the end is at time 0
        end_age = 0.0;
    }
    double start_age = end_age + node->getBranchLength();
    
    matrix->calculateTransitionProbabilities( start_age, end_age, 1.0, P);
}

vector<double> UnitMixtureModel::getRootFrequencies(int) const
{
    return matrix->getStationaryFrequencies();
}

vector<double> UnitMixtureModel::componentProbs() const { return {1.0};}

void UnitMixtureModel::scale(double factor)
{
    setRate( rate() * factor );
}

void UnitMixtureModel::setRate(double r)
{
    matrix->rescaleToAverageRate( r );
}

double UnitMixtureModel::rate() const
{
    return matrix->averageRate();
}
        
UnitMixtureModel::UnitMixtureModel(const RateMatrix& m)
    :MixtureModel(1, m.getNumberOfStates()),
     matrix( unique_ptr<RateMatrix>(m.clone()) )
{

}

UnitMixtureModel& UnitMixtureModel::operator=(const UnitMixtureModel& m)
{
    MixtureModel::operator=(m);
    matrix = unique_ptr<RateMatrix>(m.matrix->clone());
    return *this;
}

UnitMixtureModel::UnitMixtureModel(const UnitMixtureModel& m)
    :MixtureModel(m)
{
    matrix = unique_ptr<RateMatrix>(m.matrix->clone());
}




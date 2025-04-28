#include "DiscretizeGammaFromBetaQuantilesFunction.h"

#include "DistributionBeta.h"
#include "DistributionGamma.h"
#include "RbMathFunctions.h"

#include <cmath>

RevBayesCore::DiscretizeGammaFromBetaQuantilesFunction::DiscretizeGammaFromBetaQuantilesFunction(const TypedDagNode<double> *a, const TypedDagNode<double> *b, const TypedDagNode<double> *s, const TypedDagNode<double> *r, const TypedDagNode<std::int64_t> *nc, bool med) : TypedFunction< RbVector<double> >( new RbVector<double>(nc->getValue(), 1.0) ),
    alpha( a ),
    beta( b ),
    shape( s ),
    rate( r ),
    numCats(nc),
    median(med)
{

    addParameter( alpha );
    addParameter( beta );
    addParameter( shape );
    addParameter( rate );
    addParameter( numCats );

}



RevBayesCore::DiscretizeGammaFromBetaQuantilesFunction* RevBayesCore::DiscretizeGammaFromBetaQuantilesFunction::clone( void ) const
{

    return new DiscretizeGammaFromBetaQuantilesFunction(*this);
}


void RevBayesCore::DiscretizeGammaFromBetaQuantilesFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {

    if (oldP == alpha)
    {
        alpha = static_cast<const TypedDagNode<double>* >( newP );
    }

    if (oldP == beta)
    {
        beta = static_cast<const TypedDagNode<double>* >( newP );
    }

    if (oldP == shape)
    {
        shape = static_cast<const TypedDagNode<double>* >( newP );
    }

    if (oldP == rate)
    {
        rate = static_cast<const TypedDagNode<double>* >( newP );
    }

    if (oldP == numCats)
    {
        numCats = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }

}

void RevBayesCore::DiscretizeGammaFromBetaQuantilesFunction::update( void )
{
    double a = alpha->getValue();
    double b = beta->getValue();

    double s = shape->getValue();
    double r = rate->getValue();
// std::cout << "updating, a = " << a << ", b = " << b << ", s = " << s << ", r = " << r << std::endl;

    int nCats = (int)numCats->getValue();

    std::vector<double> beta_quants;
    beta_quants.push_back(0.0);
    for (size_t i=0; i<nCats-1; ++i)
    {
      double p = (i + 1.0)/nCats;
      beta_quants.push_back(RbStatistics::Beta::quantile(a, b, p));
    }
    beta_quants.push_back(1.0);

    if (median) {
        /* the median value for each category is used to represent all of the values
        in that category */
        for (size_t i=0; i<nCats; i++)
        {
            double p = beta_quants[i] + (beta_quants[i+1]-beta_quants[i])/2.0;
            (*value)[i] = RbStatistics::Gamma::quantile(s,r,p);
// std::cout << "evaluating gamma quantile at p = " << p << "; got " << (*value)[i] << std::endl;
        }

        double my_mean = 0.0;
        for (size_t i=0; i<nCats; i++)
        {
            my_mean += (*value)[i] * (beta_quants[i+1] - beta_quants[i]);
        }
// std::cout << "my mean is currently " << my_mean << std::endl;
        double factor = s / r / my_mean;
        for (size_t i=0; i<nCats; i++)
        {
            (*value)[i] *= factor;
        }
    }
    else
    {
        /* the mean value for each category is used to represent all of the values
        in that category */
        /* calculate the cumulative values */
        for (size_t i=0; i<nCats-1; i++)
        {
            (*value)[i] = RbStatistics::Gamma::quantile(s,r,beta_quants[i+1]);
      }
        for (size_t i=0; i<nCats-1; i++)
        {
            (*value)[i] = RbMath::incompleteGamma((*value)[i] * r, s + 1.0);
        }
        (*value)[nCats-1] = 1.0;

        // double factor = s / r * nCats;
        /* calculate the relative values and rescale */
        for (int i=nCats-1; i>0; i--)
        {
            (*value)[i] -= (*value)[i-1];
            (*value)[i] *= s / r / (beta_quants[i+1] - beta_quants[i]);
        }
        (*value)[0] *= s / r / beta_quants[1];
    }

}

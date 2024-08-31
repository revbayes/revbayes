#include "Trace.h"

#include <cstdlib>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "Cloner.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;


template <>
int Trace<double>::isCoveredInInterval(const std::string &v, double alpha, bool verbose)
{
    double sample = atof( v.c_str() );

    double smaller_values_count = 0;
    double equal_values_count   = 0;
    for (size_t j=0; j<values.size(); ++j)
    {
        if ( values[j] < sample )
        {
            ++smaller_values_count;
        }
        else if ( values[j] == sample )
        {
            ++equal_values_count;
        }
    }

    RandomNumberGenerator *rng = GLOBAL_RNG;

    double u = rng->uniform01();
    double quantile = (smaller_values_count + u*equal_values_count) / double(values.size());
    double lower = (1.0 - alpha) / 2.0;
    double upper = 1.0 - lower;
    bool covered = ( quantile >= lower && quantile <= upper );

    return (covered ? 0 : (quantile < lower ? -1 : 1) );
}


template <>
int Trace<long>::isCoveredInInterval(const std::string &v, double alpha, bool verbose)
{

    long sample = atof( v.c_str() );

    double smaller_values_count = 0;
    double equal_values_count   = 0;
    for (size_t j=0; j<values.size(); ++j)
    {

        if ( values[j] < sample )
        {
            ++smaller_values_count;
        }
        else if ( values[j] == sample )
        {
            ++equal_values_count;
        }

    }
    
    RandomNumberGenerator *rng = GLOBAL_RNG;

    double u = rng->uniform01();
    double quantile = (smaller_values_count + floor(u*(equal_values_count+1)) ) / double(values.size());
    double lower = (1.0 - alpha) / 2.0;
    double upper = 1.0 - lower;
    bool covered = ( quantile >= lower && quantile <= upper );

    return (covered ? 0 : (quantile < lower ? -1 : 1) );
}


template <>
int Trace<RbVector<double > >::isCoveredInInterval(const std::string &v, double i, bool verbose)
{

    RbVector<double> sample = RbVector<double>();
    sample.initFromString( v );

//    double alpha = 1.0 - std::pow(1.0-i,double(sample.size()));
    double alpha = i;

    std::vector<double> smaller_values_count = RbVector<double>(sample.size(), 0.0);
    std::vector<double> equal_values_count   = RbVector<double>(sample.size(), 0.0);
    for (size_t i=0; i<values.size(); ++i)
    {

        for (size_t j=0; j<sample.size(); ++j)
        {

            if (values[i][j] < sample[j] )
            {
                ++smaller_values_count[j];
            }
            else if ( values[i][j] == sample[j] )
            {
                ++equal_values_count[j];
            }

        }

    }


    bool covered = true;
    double num_covered = 0.0;
    for (size_t j=0; j<sample.size(); ++j)
    {
        double quantile = (smaller_values_count[j] + 0.5*equal_values_count[j]) / double(values.size());
        double lower = (1.0 - alpha) / 2.0;
        double upper = 1.0 - lower;
        if ( quantile >= lower && quantile <= upper )
        {
            ++num_covered;
        }
        covered &= ( quantile >= lower && quantile <= upper );
    }

    RandomNumberGenerator *rng = GLOBAL_RNG;
    double include_prob = num_covered / sample.size();
    covered = ( include_prob > rng->uniform01() );

    return (covered ? 0 : -1);
}


template <>
int Trace<Simplex>::isCoveredInInterval(const std::string &v, double i, bool verbose)
{

    Simplex sample = Simplex();
    sample.initFromString( v );

    //    double alpha = 1.0 - std::pow(1.0-i,double(sample.size()));
    double alpha = i;

    std::vector<double> smaller_values_count = std::vector<double>(sample.size(), 0.0);
    for (size_t i=0; i<values.size(); ++i)
    {

        for (size_t j=0; j<sample.size(); ++j)
        {

            if (values[i][j] < sample[j] )
            {
                ++smaller_values_count[j];
            }

        }

    }


    bool covered = true;
    double num_covered = 0.0;
    for (size_t j=0; j<sample.size(); ++j)
    {
        double quantile = smaller_values_count[j] / double(values.size());
        double lower = (1.0 - alpha) / 2.0;
        double upper = 1.0 - lower;
        if ( quantile >= lower && quantile <= upper )
        {
            ++num_covered;
        }
        covered &= ( quantile >= lower && quantile <= upper );
    }

    RandomNumberGenerator *rng = GLOBAL_RNG;
    double include_prob = num_covered / sample.size();
    covered = ( include_prob > rng->uniform01() );

    return (covered ? 0 : -1);
}

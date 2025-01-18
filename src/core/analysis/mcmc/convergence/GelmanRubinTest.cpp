#include "GelmanRubinTest.h"


#include "Cloner.h"
#include "TraceNumeric.h"

using namespace RevBayesCore;
using namespace std;

GelmanRubinTest::GelmanRubinTest(double r, size_t n) : ConvergenceDiagnosticContinuous(),
    R( r ),
    nBatches( n )
{
    
}

double GelmanRubinTest::assessConvergence(const TraceNumeric& trace)
{
    
    double withinBatchVariance  = 0;
    double betweenBatchVariance = 0;
    
    size_t totalSampleSize = trace.size(true);
    size_t batchSize       = size_t(totalSampleSize / (double)nBatches);
    size_t burnin          = trace.getBurnin();
    
    // get the mean between all traces
    double totalMean = trace.getMean();
    
    // get a mean and standard error for each block
    std::vector<double> batchMeans  = std::vector<double>(nBatches,0.0);
    for (size_t i=0; i<nBatches; i++)
    {
        batchMeans[i] = trace.getMean(i*batchSize+burnin, (i+1)*batchSize+burnin);
        
        // iterate over all samples from the chains
        for (size_t j=i*batchSize+burnin; j<(i+1)*batchSize+burnin; j++)
        {
            withinBatchVariance     += ( (trace.objectAt(j) - batchMeans[i])*(trace.objectAt(j) - batchMeans[i]) );
            betweenBatchVariance    += ( (trace.objectAt(j) - totalMean)*(trace.objectAt(i) - totalMean) );
        }
    }
    
    double psrf = ((totalSampleSize-nBatches) / (totalSampleSize-1.0)) * (betweenBatchVariance/withinBatchVariance);
    
    return psrf;
}

double GelmanRubinTest::assessConvergence(const std::vector<TraceNumeric>& traces)
{
    
    double within_chain_variance     = 0;
    double between_chain_variance    = 0;
    double total_mean                = 0;
    size_t total_sample_size         = 0;
    
    // get number of chains
    size_t nChains = traces.size();
    
    // get a mean and standard error for each block
    std::vector<double> chain_means  = std::vector< double >(nChains,0.0);
    for (size_t i=0; i<nChains; i++)
    {
        chain_means[i] = traces[i].getMean();

        total_mean += chain_means[i]*traces[i].size(true);
        total_sample_size += traces[i].size(true);
    }

    total_mean /= double(total_sample_size);

    // iterate over all chains
    for (size_t i=0; i<nChains; i++)
    {
        const TraceNumeric& chain    = traces[i];
        
        // iterate over all samples from the chains
        for (size_t j=chain.getBurnin(); j<chain.size(); j++)
        {
            within_chain_variance     += ( (chain.objectAt(j) - chain_means[i])*(chain.objectAt(j) - chain_means[i]) );
            between_chain_variance    += ( (chain.objectAt(j) - total_mean)*(chain.objectAt(j) - total_mean) );
        }
    }
    
    double psrf = ((total_sample_size-nChains) / (total_sample_size-1.0)) * (between_chain_variance/within_chain_variance);
    
    return psrf;
}

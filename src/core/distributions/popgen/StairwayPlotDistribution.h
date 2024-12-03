#ifndef StairwayPlotDistribution_H
#define StairwayPlotDistribution_H

#include "TypedDistribution.h"
#include "RbVector.h"

namespace RevBayesCore {
class DagNode;
class Simplex;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief StairwayPlot distribution class.
     *
     * The StairwayPlot distribution represents the StairwayPlot method by Xiaoming Liu & Yun-Xin Fu.
     * The StairwayPlot uses the SFS (site frequency spectrum) as observations and computes the
     * likelihood of a vector of theta (with theta = 4 * N_e * mu).
     * Given the theta vector and a mutation rate, we can construct the stairway plot as the
     * effective population sizes N_e and the expected coalescent times T, which are used the delimit
     * the length of the piecewise constant population size functions.
     *
     * References:
     * Liu & Fu: Exploring population size changes using SNP frequency spectra. Nature genetics, 2015, 47, 555-559
     * Liu & Fu: Stairway Plot 2: demographic history inference with folded SNP frequency spectra. Genome Biology, 2020, 21, 1-9
     *
     */
    class StairwayPlotDistribution : public TypedDistribution< RbVector<double> >, public MemberObject< RbVector<double> > {
        
    public:
        
        enum CODING { ALL, NO_MONOMORPHIC, NO_SINGLETONS };
        enum MONOMORPHIC_PROBABILITY { REST, TREE_LENGTH };
        
        StairwayPlotDistribution(const TypedDagNode< RbVector<double> > *th, long n, long n_ind, bool f=false, MONOMORPHIC_PROBABILITY mono=REST, CODING cod=ALL);
        virtual                                            ~StairwayPlotDistribution(void);                                                //!< Virtual destructor
        
        // public member functions
        StairwayPlotDistribution*                           clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;

    private:
        bool                                                calculateExpectedSFS(void) const;
        RbVector<double>                                    computeTimeBreakpoints(void) const;
        void                                                initialize(void);
        
        // members
        const TypedDagNode< RbVector<double> >*             theta;
        long                                                num_sites;
        bool                                                folded;
        long                                                num_individuals;
        mutable std::vector<double>                         expected_SFS;
        mutable std::vector< std::vector<double> >          prob_k;
        double                                              ln_factorial_num_sites_all;
        std::vector<double>                                 ln_factorial_num_sites;
        MONOMORPHIC_PROBABILITY                             monomorphic_probability;
        CODING                                              coding;
    };
    
}


#endif

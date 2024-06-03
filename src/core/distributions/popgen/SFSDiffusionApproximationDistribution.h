#ifndef SFSDiffusionApproximationDistribution_H
#define SFSDiffusionApproximationDistribution_H

#include "TypedDistribution.h"
#include "RbVector.h"

namespace RevBayesCore {
class DagNode;
class Simplex;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief SFSDiffusionApproximation distribution class.
     *
     * The SFSDiffusionApproximation distribution represents the SFSDiffusionApproximation method by Xiaoming Liu & Yun-Xin Fu.
     * The SFSDiffusionApproximation uses the SFS (site frequency spectrum) as observations and computes the
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
    class SFSDiffusionApproximationDistribution : public TypedDistribution< RbVector<double> >, public MemberObject< RbVector<double> > {
        
    public:
        
        enum CODING { ALL, NO_MONOMORPHIC, NO_SINGLETONS };
        
        SFSDiffusionApproximationDistribution(const TypedDagNode< RbVector<double> > *th, const TypedDagNode< RbVector<double> > *ls, long n, long n_ind, bool f=false, CODING cod=ALL);
        virtual                                            ~SFSDiffusionApproximationDistribution(void);                                                //!< Virtual destructor
        
        // public member functions
        SFSDiffusionApproximationDistribution*              clone(void) const;                                                          //!< Create an independent clone
        double                                              computeLnProbability(void);
        void                                                redrawValue(void);
        
    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;

    private:
        bool                                                calculateExpectedSFS(void) const;
        RbVector<double>                                    exponential_grid(long pts) const;
        RbVector<double>                                    extrap_SFS(void) const;
        RbVector<double>                                    from_phi(RbVector<double> phi, RbVector<double> grid) const;
        void                                                initialize(void);
        RbVector<double>                                    integration_one_pop(RbVector<double> phi, RbVector<double> grid, double epochlength, double nu, double theta0) const;
        RbVector<double>                                    n_epoch(long pts) const;
        RbVector<double>                                    phi_snm(double theta0, RbVector<double> grid) const;
        RbVector<double>                                    tridiag(RbVector<double> a, RbVector<double> b, RbVector<double> c, RbVector<double> r) const;

    
        // members
        const TypedDagNode< RbVector<double> >*             theta;
        const TypedDagNode< RbVector<double> >*             lengths;
        long                                                num_sites;
        bool                                                folded;
        long                                                num_individuals;
        mutable std::vector<double>                         expected_SFS;
        double                                              ln_factorial_num_sites_all;
        std::vector<double>                                 ln_factorial_num_sites;
        CODING                                              coding;
    };
    
}


#endif

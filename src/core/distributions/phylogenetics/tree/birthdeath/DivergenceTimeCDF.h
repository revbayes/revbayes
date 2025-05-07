#ifndef DivergenceTimeCDF_H
#define DivergenceTimeCDF_H

#include <cstdio>
#include <vector>

namespace RevBayesCore {

    template <class T>
    struct DivergenceTimeCDF_tolerance
    {
        DivergenceTimeCDF_tolerance(double d) { tol = d; }
        bool operator()(const T& a, const T& b) const { return (fabs(a-b) < tol);  }
        
        
        double tol;
    };
    
    class DivergenceTimeCDF {
    
        public:
    
            DivergenceTimeCDF(double o, double p, double u, double sf,
                              const std::vector<double> &r,
                              const std::vector<double> &b,
                              const std::vector<double> &d);
        
        
            double                          operator()(const double x);

        private:
        
            double                          computeProbabilitySurvival(double start, double end) const;
            size_t                          lower_index(double t) const;                                                                                  //!< Find the max index so that rateChangeTimes[index] < t < rateChangeTimes[index+1]
            size_t                          lower_index(double t, size_t min, size_t max) const;                                                                                  //!< Find the max index so that rateChangeTimes[index] < t < rateChangeTimes[index+1]
            double                          pSurvival(double start, double end, double r) const;                        //!< Compute the probability of survival of the process including uniform taxon sampling.
            double                          rateIntegral(double t_low, double t_high) const;

            double                          origin_age;
            double                          present_age;
            double                          u;

            std::vector<double>             rate_change_times;
            std::vector<double>             birth;
            std::vector<double>             death;
            double                          rho;


    };
}

#endif

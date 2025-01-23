#ifndef BDS_ODE_H
#define BDS_ODE_H

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief Birth-Death-Shift differential stepper function
     *
     */
    class BDS_ODE {
        
    public:
        
        BDS_ODE( 
                const std::vector<double> &l, 
                const std::vector<double> &m, 
                const size_t &n,
                const double &a,
                const double &b,
                const bool &f
                );
        
        void operator() ( const std::vector< double > &x, std::vector< double > &dxdt , const double t );
        
    private:
        const std::vector<double> mu;      //!< vector of extinction rate categories
        const std::vector<double> lambda;  //!< vector of speciation rate categories
        const size_t num_classes; //number of rate classes (not rate categories)
        const double alpha;
        const double beta;
        const bool forward;
    };
    
}

#endif

#include "SFSDiffusionApproximationDistribution.h"

#include <stddef.h>
#include <stdio.h>

#include "DistributionPoisson.h"
#include "RandomNumberFactory.h"
#include "RbMathFunctions.h"
#include "Cloneable.h"
#include "Simplex.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/*SFSDiffusionApproximation Distribution Constructor
 * @param p A simplex of the the probabilities for each category
 * @param n A long for the number of trials
 */

SFSDiffusionApproximationDistribution::SFSDiffusionApproximationDistribution(const TypedDagNode< RbVector<double> > *th, const TypedDagNode< RbVector<double> > *ls, long n, long n_ind, bool f, CODING c) : TypedDistribution< RbVector<double> >( new RbVector<double>( f ? (n_ind/2)+1 : n_ind, 1 ) ),
    theta( th ),
    lengths( ls ),
    num_sites( n ),
    folded( f ),
    num_individuals( n_ind ),
    coding( c )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( theta );
    addParameter( lengths );
    
    initialize();

    // calculateExpectedSFS();
    // long pts = 10; // 20, 30
    // std::cout << "pts: " << pts << "\n";  

    // RbVector<double> xx = exponential_grid(pts);
    // // for (size_t i = 0; i < xx.size(); i++)
    // // {
    // //     std::cout << xx[i] << "\t";
    // // }
    // // std::cout << "\n";
    // std::cout << "grid size: " << xx.size() << "\n";  

    // double theta0 = 1;
    // RbVector<double> phi = phi_snm(theta0, xx);
    // // for (size_t i = 0; i < phi.size(); i++)
    // // {
    // //     std::cout << phi[i] << "\t";
    // // }
    // // std::cout << "\n"; 
    // std::cout << "phi size: " << phi.size() << "\n";

    // // phi, grid, epochlength, nu, theta0
    // RbVector<double> phi_new = integration_one_pop(phi, xx, 0.7, 1.3, theta0);

    // // for (int i = 0; i < phi_new.size(); i++)
    // // {
    // //     std::cout << phi_new[i] << "\t";
    // // }
    // // std::cout << "\n";
    // std::cout << "phi_new size: " << phi_new.size() << "\n";

    // RbVector<double> model = from_phi(phi_new, xx);

    // RbVector<double> model = n_epoch(pts);
    // RbVector<double> model = extrap_SFS();
    // for (int i = 0; i < model.size(); i++)
    // {
    //     std::cout << model[i] << "\t";
    // }
    // std::cout << "\n";
    // std::cout << "model size: " << model.size() << "\n";

    // calculateExpectedSFS();
    // computeLnProbability();

}


SFSDiffusionApproximationDistribution::~SFSDiffusionApproximationDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too.
    // The model needs to do that!
}


bool SFSDiffusionApproximationDistribution::calculateExpectedSFS(void) const
{

    // special case one epoch?

    // What do I need?
    //
    // extrapolating function
    //
    // exponential grid -> need to set pts
    // phi -> integration
    // model_from_phi


    // get the thetas and epoch lengths for easier handling
    const RbVector<double>& th = theta->getValue();
    const RbVector<double>& ls = lengths->getValue();

    // long n_ind = num_individuals;

    // double theta_0 = th[0];
    // RbVector<long> pts = (n_ind + 10, n_ind + 20, n_ind + 30);

    // combine in n_epoch function - then extrapolate?
    // xx = exponential_grid(pts);
    // phi = phi_snm (theta0, xx);

    // for (size_t i = 1; i<th.size(); ++i)
    // {
    //     phi = integration_one_pop_(phi, xx, T = ls[i-1], nu = th[i]/theta_0, theta0 = theta_0);
    //     //std::cout << "epoch length " << ls[i-1] << " theta " << th[i] << "\n";
    // }

    // expected_SFS = from_phi(phi, n_ind, xx);

    RbVector<double> extrapolated_SFS = extrap_SFS();

    // get the data, i.e., the observed counts for the frequencies
    const RbVector<double>& obs_sfs_counts = *value;
    // double L_obs = std::accumulate(obs_sfs_counts.begin(), obs_sfs_counts.end(), 0.0);
    double L_obs = 0.0;

    for (int x = 0; x < obs_sfs_counts.size()-1; x++)
    {
        L_obs += obs_sfs_counts[x];
    }

    // std::cout << "L_obs: " << L_obs << "\n";
    // for (int i = 0; i < obs_sfs_counts.size(); i++)
    // {
    //     std::cout << obs_sfs_counts[i] << "\t";
    // }
    // std::cout << "\n";

    // for (int i = 0; i < extrapolated_SFS.size(); i++)
    // {
    //     std::cout << extrapolated_SFS[i] << "\t";
    // }
    // std::cout << "\n";
    // for (int i = 0; i < expected_SFS.size(); i++)
    // {
    //     std::cout << expected_SFS[i] << "\t";
    // }
    // std::cout << "\n";

    // "normalize" SFS by L:
    expected_SFS[0] = 1;
    for (int x = 1; x < expected_SFS.size(); x++)
    {
        expected_SFS[x] = extrapolated_SFS[x] / L_obs;
        if (x < expected_SFS.size()-1)
        {
            expected_SFS[0] -= expected_SFS[x];
        }
    }
    // expected_SFS[0] = 1 - std::accumulate(expected_SFS.begin()+1, expected_SFS.end(), 0.0);

    // for (int i = 0; i < expected_SFS.size(); i++)
    // {
    //     std::cout << expected_SFS[i] << "\t";
    // }
    // std::cout << "\n";


    // variable for the sum of all non-monomorphic frequencies
    double sum_expected_frequency = 0.0;
    
    // @Ronja:
    // here we need to compute a vector (with n+1 elements where n is the number of individuals)
    // of the expected SFS
    
    for (size_t i=1; i<num_individuals; ++i)
    {
        // @Ronja: here goes the model computation
        
//        double expected_frequency = 0.0;
//        for (size_t k=2; k<=(num_individuals-i+1); ++k)
//        {
//            double this_theta = th[k-2];
//            expected_frequency += this_theta * prob_k[i-1][k-2];
//        }
//        
       // add the current expected frequency to our sum
    //    sum_expected_frequency += expected_frequency;
       sum_expected_frequency += expected_SFS[i];
//        
//        // store the value of the expected frequency
//        expected_SFS[i] = expected_frequency;
    }
    
    // return false that the likelihood should be -Inf
    if ( sum_expected_frequency > 1 )
    {
        return false;
    }
    
    if ( coding != ALL )
    {
        double correction = 1.0;
        size_t min_allele_count = 1;
        size_t max_allele_count = num_individuals-1;

        if ( coding == NO_MONOMORPHIC )
        {
            correction = 1.0 - expected_SFS[0];
        }
        else if ( coding == NO_SINGLETONS )
        {
            correction = 1.0 - expected_SFS[0] - expected_SFS[1] - expected_SFS[num_individuals-1];
            min_allele_count = 2;
            max_allele_count = num_individuals-2;
        }
        
        // now normalize
        if ( coding == NO_MONOMORPHIC )
        {
            for (size_t i=min_allele_count; i<=max_allele_count; ++i)
            {
                // store the corrected frequency
                expected_SFS[i] = expected_SFS[i] / correction;
            }
        }
    }
    
    // return that our expected frequencies worked
    return true;
}



SFSDiffusionApproximationDistribution* SFSDiffusionApproximationDistribution::clone( void ) const
{
    return new SFSDiffusionApproximationDistribution( *this );
}


double SFSDiffusionApproximationDistribution::computeLnProbability( void )
{
    
    // initialize the probability
    double ln_prob = 0;
    
    // get the data, i.e., the observed counts for the frequencies
    const RbVector<double>& obs_sfs_counts = *value;
    // const RbVector<double> obs_sfs_counts (num_individuals+1,1.0);
    // const RbVector<double> obs_sfs_counts = extrap_SFS();

    // compute the expected SFS, i.e., the expected frequency of observing a site with frequency i
    bool success = calculateExpectedSFS();
    if ( success == false )
    {
        return RbConstants::Double::neginf;
    }
    
    size_t max_freq = num_individuals;
    if ( folded == true )
    {
        max_freq = floor( num_individuals / 2.0 ) + 1;
    }

    // check for the coding
    // only add the monomorphic probability of we use the coding "all"
    if ( coding == ALL )
    {
        ln_prob = (double)obs_sfs_counts[0] * log(expected_SFS[0]);
    }

    // shift the smallest allele count depending on coding
    size_t smallest_allele_count = 1;
    if ( coding == NO_SINGLETONS )
    {
        
        smallest_allele_count = 2;
        // also shift the max allele count if we don't allow for singletons
        if ( folded == false )
        {
            max_freq = num_individuals-1;
            ln_prob = (double)(obs_sfs_counts[0]+obs_sfs_counts[1]) * log(expected_SFS[0]+expected_SFS[1]);
        }
        else
        {
            ln_prob = (double)(obs_sfs_counts[0]+obs_sfs_counts[1]) * log(expected_SFS[0]+expected_SFS[1]+expected_SFS[num_individuals-1]);
        }
    }
    
    // compute the probability for all allele frequency counts
    for (size_t i=smallest_allele_count; i<max_freq; ++i)
    {
        
        // compute the multinomial probability for the SFS frequency
        if ( folded == false )
        {
            // ln_prob -= ln_factorial_num_sites[i];
            ln_prob += (double)obs_sfs_counts[i] * log(expected_SFS[i]);
        }
        else
        {
            // ln_prob -= ln_factorial_num_sites[i];
            if ( i == (num_individuals/2.0) )
            {
                ln_prob += (double)obs_sfs_counts[i] * log(expected_SFS[i]);
            }
            else
            {
                ln_prob += (double)obs_sfs_counts[i] * log(expected_SFS[i]+expected_SFS[num_individuals-i]);
            }
        }
    }
    
    // divide by the factorial of the total number of observation
    // ln_prob -= ln_factorial_num_sites_all;
    // std::cout << "ln_prob: " << ln_prob << "\n";
    
    return ln_prob;
}



void SFSDiffusionApproximationDistribution::executeMethod(const std::string &name, const std::vector<const DagNode *> &args, RbVector<double> &rv) const
{
   
    if ( name == "getExpectedAlleleFrequencies" )
    {
        rv = expected_SFS;
    }
    else
    {
        throw RbException("The SFSDiffusionApproximation does not have a member method called '" + name + "'.");
    }

}


void SFSDiffusionApproximationDistribution::initialize( void )
{
    
    // allocate/resize the expected SFS frequency vector
    expected_SFS = std::vector<double>( num_individuals+1, 0.0 );
        
    ln_factorial_num_sites.clear();
    ln_factorial_num_sites.resize( folded ? floor( num_individuals / 2.0 ) + 1 : num_individuals );

    
    // get the data, i.e., the observed counts for the frequencies
    const RbVector<double>& obs_sfs_counts = *value;
    
    for (size_t i=1; i<num_individuals; ++i)
    {
        if ( folded == false || i <= (num_individuals/2.0) )
        {
            ln_factorial_num_sites[i] = RbMath::lnGamma( obs_sfs_counts[i] + 1 );
        }
    }
    ln_factorial_num_sites[0] = RbMath::lnGamma( obs_sfs_counts[0] + 1 );
    
    ln_factorial_num_sites_all = RbMath::lnGamma( num_sites + 1 );
}


void SFSDiffusionApproximationDistribution::redrawValue( void )
{
    *value = RbVector<double>( num_individuals+1, 1 );
}


/** Swap a parameter of the distribution */
void SFSDiffusionApproximationDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == theta)
    {
        theta = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
}

RbVector<double> SFSDiffusionApproximationDistribution::extrap_SFS(void) const
{
    long n_ind = num_individuals;
    // std::cout << (n_ind + 10) << "\n";
    // RbVector<long> pts = ((n_ind + 10), (n_ind + 20), (n_ind + 30));
    RbVector<long> pts;

    for (int i = 10; i <= 30; i += 10)
    {
        pts.push_back(n_ind+i);
    }

    // std::cout << pts[0] << "\t" << pts[1] << "\t" << pts[2] << "\n";

    RbVector<double> y1 = n_epoch(pts[0]);
    RbVector<double> y2 = n_epoch(pts[1]);
    RbVector<double> y3 = n_epoch(pts[2]);

    RbVector<double> y1_log(y1.size());
    RbVector<double> y2_log(y2.size());
    RbVector<double> y3_log(y3.size());

    for (int x = 0; x < y1_log.size(); x++)
    {
        y1_log[x] = log(y1[x]);
        y2_log[x] = log(y2[x]);
        y3_log[x] = log(y3[x]);
        // std::cout << y1_log[x] << "\t" << y2_log[x] << "\t" << y3_log[x] << "\n";
    }

    double xl1 = exponential_grid(pts[0])[1];
    double xl2 = exponential_grid(pts[1])[1];
    double xl3 = exponential_grid(pts[2])[1];

    // std::cout << xl1 << "\t" << xl2 << "\t" << xl3 << "\n";

    RbVector<double> ex_result;

    for (int x = 0; x < y1_log.size(); x++)
    {
        double ex_result_value = xl2 * xl3 / ((xl1-xl2)*(xl1-xl3)) * y1_log[x] + xl1 * xl3 / ((xl2-xl1)*(xl2-xl3)) * y2_log[x] + xl1 * xl2 / ((xl3-xl1)*(xl3-xl2)) * y3_log[x];
        ex_result.push_back(ex_result_value);
        // std::cout << ex_result[x] << "\t";
    }
    // std::cout << "\n";

    RbVector<double> SFS_extrap;

    for (int x = 0; x < ex_result.size(); x++)
    {
        SFS_extrap.push_back(exp(ex_result[x]));
    }

    return (SFS_extrap);
}

RbVector<double> SFSDiffusionApproximationDistribution::n_epoch(long pts) const
{
    const RbVector<double>& th = theta->getValue();
    const RbVector<double>& ls = lengths->getValue();

    long n_ind = num_individuals;

    double theta_0 = th[0];

    RbVector<double> xx = exponential_grid(pts);
    RbVector<double> phi = phi_snm(theta_0, xx);

    for (size_t i = 1; i<th.size(); ++i)
    {
        phi = integration_one_pop(phi, xx, ls[i-1], th[i]/theta_0, theta_0);
    }

    RbVector<double> model = from_phi(phi, xx);

    return model;
}

RbVector<double> SFSDiffusionApproximationDistribution::exponential_grid(long pts) const
{
    double crwd = 8.0;

    double step = 2.0 / (pts-1);
    RbVector<double> grid;
    // RbVector<double> unif;
    double value_grid;
    double gridpoint;
    for (long i = 0; i < pts; i++)
    {
        value_grid = -1.0 + i * step;
        // std::cout << value_grid << "\t";
        //unif.push_back(value_grid);
        gridpoint = 1 / (1 + exp(-crwd * value_grid));
        // std::cout << gridpoint << "\t";
        grid.push_back(gridpoint);
    }
    // std::cout << "\n";

    // normalize
    RbVector<double> grid_norm;
    for (long i = 0; i < pts; i++)
    {
        gridpoint = (grid[i] - grid[0]) / (grid[pts-1] - grid[0]);
        // std::cout << gridpoint << "\t";
        grid_norm.push_back(gridpoint);
    }
    // std::cout << "\n";

    return grid_norm;
}

RbVector<double> SFSDiffusionApproximationDistribution::phi_snm(double theta0, RbVector<double> grid) const
{
    // double nu = 1.0;

    RbVector<double> phi;

    if (grid[0] == 0)
    {
        phi.push_back(theta0 / grid[1]);
        for (size_t x = 1; x < grid.size(); x++)
        {
            phi.push_back(theta0 / grid[x]);
        }
    }
    else
    {
        for (size_t x = 0; x < grid.size(); x++)
        {
            phi.push_back(theta0 / grid[x]);
        }
    }
    return phi;
}

RbVector<double> SFSDiffusionApproximationDistribution::integration_one_pop(RbVector<double> phi, RbVector<double> grid, double epochlength, double nu, double theta0) const
{
    RbVector<double> M(grid.size(), 0.0);
    // for (size_t x = 0; x < grid.size(); x++)
    // {
    //     M.push_back(0.0);
    // }
    // std::cout << "M size: " << M.size() << "\n";   
    RbVector<double> MInt = M;

    RbVector<double> V;
    for (size_t x = 0; x < grid.size(); x++)
    {
        double V_value = 1./nu * grid[x] * (1-grid[x]);
        // std::cout << V_value << "\t";
        V.push_back(V_value);
    }
    // std::cout << "\n";
    // std::cout << "V size: " << V.size() << "\n";

    RbVector<double> VInt;
    for (size_t x = 0; x < (grid.size()-1); x++)
    {
        double VInt_value = 1./nu * (grid[x] + grid[x+1])/2 * (1-(grid[x] + grid[x+1])/2);
        // std::cout << VInt_value << "\t";
        VInt.push_back(VInt_value);
    }
    // std::cout << "\n";
    // std::cout << "VInt size: " << VInt.size() << "\n";    
 
    // RbVector<double> diff = std::adjacent_difference(grid.begin(), grid.end(), grid.begin());
    RbVector<double> dx;
    for (size_t x = 0; x < (grid.size()-1); x++)
    {
        double dx_value = grid[x+1] - grid[x];
        dx.push_back(dx_value);
        // std::cout << dx[x] << "\t";
    }   
    // std::cout << "\n";
    // std::cout << "dx size: " << dx.size() << "\n";   

    RbVector<double> dfactor;
    dfactor.push_back(2/dx[0]);
    for (size_t x = 0; x < (dx.size()-1); x++)
    {
        double dfactor_value = 2 / (dx[x] + dx[x+1]);
        dfactor.push_back(dfactor_value);
    }
    dfactor.push_back(2/dx[dx.size()-1]);

    // for (size_t x = 0; x < dfactor.size(); x++)
    // {
    //     std::cout << dfactor[x] << "\t";
    // }    
    // std::cout << "\n";
    // std::cout << "dfactor size: " << dfactor.size() << "\n";   
    
    double delj = 0.5;

    RbVector<double> a;
    RbVector<double> c;
    RbVector<double> b;

    a.push_back(0.0);
    for (size_t x = 0; x < (phi.size()-1); x++)
    {
        // std::cout << "df: " << dfactor[x+1] << ", V: " << V[x] << ", dx: " << dx[x] << "\n";
        double a_value = dfactor[x+1] * (-MInt[x] * delj - V[x] / (2 * dx[x]));
        a.push_back(a_value);
        
        double c_value = - dfactor[x] * (-MInt[x] * (1-delj) + V[x+1] / (2 * dx[x]));
        c.push_back(c_value);

        double b_value = - dfactor[x] * (-MInt[x] * delj - V[x] / (2 * dx[x]));
        if (x > 0)
        {
            b_value += dfactor[x] * (-MInt[x-1] * (1-delj) + V[x] / (2 * dx[x-1]));
        }
        b.push_back(b_value);
    }
    c.push_back(0.0);
    double b_final = dfactor[phi.size()-1] * (-MInt[phi.size()-2] * delj + V[phi.size()-1] / (2 * dx[phi.size()-2]));
    b.push_back(b_final);

    if (M[0] <= 0)
    {
        b[0] = b[0] + (0.5/nu - M[0]) * 2/dx[0];
    }
    if (M[M.size()-1] <= 0)
    {
        b[b.size()-1] = b[b.size()-1] + (0.5/nu - M[M.size()-1]) * 2/dx[dx.size()-1];
    }   
    // for (size_t x = 0; x < phi.size(); x++)
    // {
    //     std::cout << "a: " << a[x] << ", c: " << c[x] << ", b: " << b[x] << "\n";
    // }
    // std::cout << "\n";
    // std::cout << "a size: " << a.size() << "\n";   
    // std::cout << "c size: " << c.size() << "\n";   
    // std::cout << "b size: " << c.size() << "\n";   

    double dt;
    double timescale_factor = 1e-3;

    // RbVector<double> VMvector = (0.25 / nu, 0.0, 0.0);
    // double maxVM = *std::max_element(VMvector.begin(), VMvector.end());
    double maxVM = 0.25 / nu;
    // std::cout << maxVM << "\n";

    if (maxVM > 0)
    {
        dt = timescale_factor / maxVM;
    }
    else
    {
        dt = std::numeric_limits<double>::infinity();
    }

    double current_t = 0.0;

    double this_dt;
    RbVector<double> phi_new = phi;
    RbVector<double> r(phi_new.size());
    RbVector<double> b_tridiag = b;
    while (current_t < epochlength)
    {
        this_dt = std::min(dt, epochlength - current_t);
        // std::cout << this_dt << "\n";
        phi_new[1] = phi_new[1] + this_dt / grid[1] * theta0/2.0 * 2.0/(grid[2] - grid[0]);
        // std::cout << phi_new[1];
        for (int x = 0; x < phi_new.size(); x++)
        {
            r[x] = phi_new[x]/this_dt;
            // std::cout << phi_new[x] << "\t";
            // std::cout << r[x] << "\t";
        }
        for (int x = 0; x < b.size(); x++)
        {
            b_tridiag[x] = b[x]+1/this_dt;
            // std::cout << b_tridiag[x] << "\t";
        }
        // std::cout << "\n";
        phi_new = tridiag(a, b_tridiag, c, r);
        current_t += this_dt;
    }

    return phi_new;
}


// copied and modified from the dadi bitbucket https://bitbucket.org/gutenkunstlab/dadi/src/master/dadi/
// https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
RbVector<double> SFSDiffusionApproximationDistribution::tridiag(RbVector<double> a, RbVector<double> b, RbVector<double> c, RbVector<double> r) const
{
    /*
    Based on Numerical Recipes in C tridiag function.
    */
    int n = a.size();
    std::vector<double> gam(n);
    std::vector<double> u(n);

    double bet = b[0];
    int j;

    u[0] = r[0]/bet;
    for(j=1; j <= n-1; j++){
        gam[j] = c[j-1]/bet;
        bet = b[j] - a[j]*gam[j];
        u[j] = (r[j]-a[j]*u[j-1])/bet;
    }
    
    for(j=(n-2); j >= 0; j--){
        u[j] -= gam[j+1]*u[j+1];
    }

    return(u);
}


RbVector<double> SFSDiffusionApproximationDistribution::from_phi(RbVector<double> phi, RbVector<double> grid) const
{
    double nsamp = num_individuals;

    RbVector<double> data (nsamp+1, 0.0);
    RbVector<double> xx = grid;
    for (int x = 0; x < xx.size(); x++)
    {
        xx[x] = std::min( std::max(xx[x],0.0) , 1.0);
        // std::cout << xx[x] << "\t";
    }
    RbVector<double> s;
    for (int x = 0; x < phi.size()-1; x++)
    {
        double s_value = (phi[x+1] - phi[x]) / (xx[x+1] - xx[x]);
        s.push_back(s_value);
        // std::cout << s[x] << "\t";
    }
    RbVector<double> c1;
    for (int x = 0; x < (phi.size()-1); x++)
    {
        double c1_value = (phi[x] - s[x] * xx[x]) / (nsamp + 1);
        c1.push_back(c1_value);
        // std::cout << c1[x] << "\t";        
    }
    RbVector<double> c2(s.size());
    RbVector<double> beta1(xx.size());
    RbVector<double> beta2(xx.size());
    RbVector<double> entries(xx.size());
    for (int d = 0; d < data.size(); d++)
    {
        for (int x = 0; x < s.size(); x ++)
        {
            double c2_value = s[x] * (d+1) / ((nsamp + 1) * (nsamp + 2));
            c2[x] = c2_value;
            // std::cout << c2[x] << "\t";
        }
        for (int x = 0; x < xx.size(); x ++)
        {
            double beta1_value = RbMath::incompleteBeta(d+1, nsamp-d+1, xx[x]);
            beta1[x] = beta1_value;
            double beta2_value = RbMath::incompleteBeta(d+2, nsamp-d+1, xx[x]);;
            beta2[x] = beta2_value;
            // std::cout << beta2[x] << "\t";
        }
        for (int x = 0; x < (xx.size()-1); x++){
            double entries_value = c1[x] * (beta1[x+1] - beta1[x]) + c2[x] * (beta2[x+1] - beta2[x]);
            entries[x] = entries_value;
            // std::cout << entries[x] << "\t";
        }
        // std::cout << "\n";
        data[d] = std::accumulate(entries.begin(), entries.end(),0.0);
    }
    // std::cout << "\n";

    return (data);
}
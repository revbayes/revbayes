#include "OccurrenceBirthDeathProcess.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include "RbMathMatrix.h"
#include "RbVector.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "TimeInterval.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;


/**
 * Constructor.
 * We delegate most parameters to the base class and initialize the members.
 *
 * \param[in]    t              Time of the origin/present/length of the process.
 * \param[in]    l              Speciation rate.
 * \param[in]    m              Extinction rate.
 * \param[in]    p              Extinction sampling rate.
 * \param[in]    o              Rate of occurrence observations.
 * \param[in]    rho            Sampling probability at present time.
 * \param[in]    r              Removal probability after sampling.
 * \param[in]    cdt            Condition of the process (none/survival/#Taxa).
 * \param[in]    tn             Taxa.
 * \param[in]    uo             If true ra is the origin time otherwise the root age of the process.
 */
OccurrenceBirthDeathProcess::OccurrenceBirthDeathProcess(                         const TypedDagNode<double> *t,
                                                                                  const TypedDagNode<double> *l,
                                                                                  const TypedDagNode<double> *m,
                                                                                  const TypedDagNode<double> *p,
                                                                                  const TypedDagNode<double> *o,
                                                                                  const TypedDagNode<double> *rho,
                                                                                  const TypedDagNode<double> *r,


                                                                                  const std::string& cdt,
                                                                                  const std::vector<Taxon> &tn,
                                                                                  const TypedDagNode< RbVector<double> > *tau,
                                                                                  bool uo,

                                                                                  TypedDagNode<Tree> *tr) : AbstractBirthDeathProcess( t, cdt, tn, uo ),
    tor( t ),
    lambda( l ),
    mu( m ),
    psi( p ),
    omega( o ),
    rho( rho ),
    removalPr( r ),
    dn_time_points ( tau )

{
    addParameter( tor );
    addParameter( lambda );
    addParameter( mu );
    addParameter( psi );
    addParameter( rho );


    if (tr != NULL)
    {
      delete value;
      value = &(tr->getValue());
    }
    else
    {
      simulateTree();
    }

}


OccurrenceBirthDeathProcess::~OccurrenceBirthDeathProcess( void ){
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
OccurrenceBirthDeathProcess* OccurrenceBirthDeathProcess::clone( void ) const
{

    return new OccurrenceBirthDeathProcess( *this );
}



/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double OccurrenceBirthDeathProcess::computeLnProbabilityDivergenceTimes( void ) const
{
  std::cout << "computeLnProbabilityDivergenceTimes in " << std::endl;

    // prepare the probability computation
    prepareProbComputation();
  std::cout << "compute ln prob prepared" << std::endl;

    // step 1: Create the list of all event times
    poolTimes();
    std::cout << "pool times" << std::endl;

    // variable declarations and initialization
    //double a1 = ComputeLt();
  	double lnProbTimes = ComputeMt();
    std::cout << "ComputeMt" << std::endl;

  //  double lnProbTimes = computeLnProbabilityTimes();

    return lnProbTimes;
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 * \return    The log-probability density.
 */


void OccurrenceBirthDeathProcess::poolTimes( void ) const
{
    // get node/time variables
    size_t num_nodes = value->getNumberOfNodes();
    std::cout << "number of nodes" << std::endl;
    std::cout << num_nodes << std::endl;

    int extant = 0;

    // classify nodes
    events.clear();
    std::cout << tor->getValue() << std::endl ;
    events.push_back(Event(tor->getValue(), "origin", 0));

    for (size_t i = 0; i < num_nodes; i++)
    {
        const TopologyNode& n = value->getNode( i );
//isFossil is an optional condition to obtain sampled ancestor node ages
/* p is a sampled ancestor
o sampled ancestor
r extant leaf
careful with internal nodes !!
q = "false" bifurcation
s = "true" bifurcation

 __|s__
|      |
o      |
       |
      q|___ p
       |
      r|

 1. Pick a fossil among those with brl > 0 (prob = 1/m)
 2. Set brl = 0
 */

        if ( n.isFossil() && n.isSampledAncestor() )
        {
            // node is sampled ancestor
            events.push_back(Event(n.getAge(),"sampled ancestor",0)) ;

        }
        else if ( n.isFossil() && !n.isSampledAncestor() )
        {
            // node is fossil leaf
            // named terminal non-removed in Lt
            // events.push_back(Event(n.getAge(),"fossil leaf", 0) ;
            events.push_back(Event(n.getAge(),"terminal non-removed",0)) ;
        }
        else if ( n.isTip() && !n.isFossil() )
        {
            // node is extant leaf
            // events.push_back(Event(n.getAge(),"extant leaf", 0) ;
            //events.push_back(Event(n.getAge(), "terminal non-removed", 0)) ;
            extant++;
        }
        else if ( n.isInternal() && !n.getChild(0).isSampledAncestor() && !n.getChild(1).isSampledAncestor() )
        {
            if (!n.isRoot())
            {
                // node is a "true" bifurcation event
                events.push_back(Event(n.getAge(),"branching time",0)) ;
            }
        }

    }

}

  //Pool all observations together into

    // // add the log probability for the serial sampling events
    // if (ps == 0.0)
    // {
    //     if ( serial_tip_ages.size() + sampled_ancestors_ages.size() > 0 )
    //     {
    //         return RbConstants::Double::neginf;
    //         //throw RbException("The serial sampling rate is zero, but the tree has serial sampled tips.");
    //     }
    // }
    // else
    // {
    //     lnProbTimes += (serial_tip_ages.size() + sampled_ancestors_ages.size() ) * log( ps );
    // }
    //
    // // add the log probability for sampling the extant taxa
    // if (num_extant_taxa > 0)
    // {
    //     lnProbTimes += num_extant_taxa * log( 4.0 * rh );
    // }
    //
    // // add the log probability of the initial sequences
    // lnProbTimes += -lnQ(process_time, c1, c2) * num_initial_lineages;
    //
    // // add the log probability for the internal node ages
    // lnProbTimes += internal_node_ages.size() * log( birth );
    // for (size_t i=0; i<internal_node_ages.size(); i++)
    // {
    //     lnProbTimes -= lnQ(internal_node_ages[i], c1, c2);
    // }
    //
    // // add the log probability for the serial tip ages
    // for (size_t i=0; i < serial_tip_ages.size(); i++)
    // {
    //     double t = serial_tip_ages[i];
    //     lnProbTimes += log(pZero(t, c1, c2)) + lnQ(t, c1, c2);
    // }
    //
    // // condition on survival
    // if ( condition == "survival")
    // {
    //     lnProbTimes -= num_initial_lineages * log(1.0 - pHatZero(process_time));
    // }
    // // condition on nTaxa
    // else if ( condition == "nTaxa" )
    // {
    //     lnProbTimes -= lnProbNumTaxa( value->getNumberOfTips(), 0, process_time, true );
    // }



double OccurrenceBirthDeathProcess::lnProbTreeShape(void) const
{
    // the birth death divergence times density is derived for a (ranked) unlabeled oriented tree
    // so we convert to a (ranked) labeled non-oriented tree probability by multiplying by 2^{n+m-1} / n!
    // where n is the number of extant tips, m is the number of sampled extinct tips

    int num_taxa = (int)value->getNumberOfTips();
    std::cout << "Number of tips" << std::endl;
    int num_extinct = (int)value->getNumberOfExtinctTips();
    std::cout << "Number of extinct tips" << std::endl;

    int num_sa = (int)value->getNumberOfSampledAncestors();
    std::cout << "Number of sampled ancestors" << std::endl;


    return (num_taxa - num_sa - 1) * RbConstants::LN2 - RbMath::lnFactorial(num_taxa - num_extinct);
}


/**
 * Compute the probability of survival if the process starts with one species at time start and ends at time end.
 *
 *
 * \param[in]    start      Start time of the process.
 * \param[in]    end        End/stopping time of the process.
 *
 * \return Speciation rate at time t.
 */
 double OccurrenceBirthDeathProcess::pSurvival(double start, double end) const
 {
        return 1.0 - pHatZero(end);
 }



/**
 * Simulate new speciation times.
 */
double OccurrenceBirthDeathProcess::simulateDivergenceTime(double origin, double present) const
{
    // incorrect placeholder for constant FBDP

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // get the parameters
    double age = origin - present;
    double b = lambda->getValue();
    double d = mu->getValue();
    double r = rho->getValue();

    // get a random draw
    double u = rng->uniform01();


    // compute the time for this draw
    // see Hartmann et al. 2010 and Stadler 2011

    double t = 0.0;
    if ( b > d )
    {
        t = ( log( ( (b-d) / (1 - (u)*(1-((b-d)*exp((d-b)*age))/(r*b+(b*(1-r)-d)*exp((d-b)*age) ) ) ) - (b*(1-r)-d) ) / (r * b) ) )  /  (b-d);
    }
    else
    {
        t = ( log( ( (b-d) / (1 - (u)*(1-(b-d)/(r*b*exp((b-d)*age)+(b*(1-r)-d) ) ) ) - (b*(1-r)-d) ) / (r * b) ) )  /  (b-d);
    }

    return present + t;
}


// double OccurrenceBirthDeathProcess::pZero(double t, double c1, double c2) const
// {
// //work in progress function u and p
//     // get helper variables
//     double a = birth - death - ps;
//     double c1 = std::fabs(sqrt(a * a + 4 * birth * ps));
//     double c2 = -(a - 2 * birth * rh) / c1;
//     // roots of the polynmial ODE
//     double Delta = pow(gamma,2) - 4*birth*death;
//     double x1 = (gamma - sqrt(Delta))/(2*birth);
//     double x2 = (gamma + sqrt(Delta))/(2*birth);
//
//     double b = lambda->getValue();
//     double d = mu->getValue();
//     double f = psi->getValue();
//     double v1 = exp(-c1*t) * (1.0-c2) - (1.0+c2);
//     double v2 = exp(-c1*t) * (1.0-c2) + (1.0+c2);
//     double v = (b + d + f + c1 * (v1 / v2)) / (2.0 * b);
//     return v;
// }



// double OccurrenceBirthDeathProcess::lnQ(double t, double c1, double c2) const
// {
//     //return log( 2*(1-c2*c2) + exp(-c1*t)*(1-c2)*(1-c2) + exp(c1*t)*(1+c2)*(1+c2) );
//
//     // numerically safe code
//     return c1*t + 2 * log( exp(-c1*t) * (1-c2) + (1+c2));
// }


double OccurrenceBirthDeathProcess::pHatZero(double t) const
 {
    double b = lambda->getValue();
    double d = mu->getValue();
    double r = rho->getValue();
    double val = r * (b-d) / (b*r + (b*(1-r)-d)*exp(-(b-d)*t));
    return 1.0 - val;
}




/**
 * Swap the parameters held by this distribution.
 *
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void OccurrenceBirthDeathProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == lambda)
    {
        lambda = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == mu)
    {
        mu = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == psi)
    {
        psi = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == rho)
    {
        rho = static_cast<const TypedDagNode<double>* >( newP );
    }
    else
    {
        // delegate the super-class
        AbstractBirthDeathProcess::swapParameterInternal(oldP, newP);
    }

}


double OccurrenceBirthDeathProcess::computeLnProbabilityTimes( void ) const
 {
 double lnProbTimes = ComputeMt();
 std::cout << "ComputeLnProbabilitytimes" << std::endl;
 std::cout << lnProbTimes << std::endl;
 return lnProbTimes;
 }

//from Warnock & Manceau computeLt function
double OccurrenceBirthDeathProcess::ComputeMt( void ) const
{
	// order times oldest to youngest -> revisit this
	std::sort( events.begin(), events.end(), AgeCompareReverse() );
  std::cout << "events" << std::endl;
  for ( int i = 0; i < events.size(); i++){
	std::cout << events[i].type	<< std::endl;

	}



	const double birth = lambda->getValue();
  std::cout << "lambda" << std::endl;
  std::cout << birth << std::endl;

	const double death = mu->getValue();
  std::cout << "mu" << std::endl;
  std::cout << death << std::endl;

	const double ps = psi->getValue();
  std::cout << "psi" << std::endl;
  std::cout << ps << std::endl;

	const double om = omega->getValue();
	const double rh = rho->getValue();
	const double rp = removalPr->getValue();

	const double gamma = birth + death + ps + om;

	size_t N = 20; // accuracy of the algorithm
	size_t S = dn_time_points->getValue().size();

	// step 2. initialize an empty matrix
	MatrixReal B(S, (N + 1), 0.0);

	size_t k = 1;

	// step 3. // this is the first entry of B
	RbVector<double> Mt;
	for ( int i = 0; i < (N + 1); i++){
		if (i == 0) Mt.push_back( 1.0 );
		else Mt.push_back( 0.0 );
	}

	// step 4.
	//-> calculate the state of the process just before the next time point

	MatrixReal A( (N+1), (N+1), 0.0 );

	double thPlusOne = tor->getValue();

	size_t indxJ = S-1;

	for(int h = 0; h < events.size(); h++){
		// deal with t > tor

		double th = events[h].time;

		if(th > tor->getValue()) continue;

		if( th != thPlusOne ){

			for(int i = 0; i < (N + 1); i++){
				for(int j = 0; j < (N + 1); j++){
					if(j == i) A[i][i] = gamma * (k + i);
					else if (j == i+1) A[i][i+1] = -death * (i + 1);
					else if (j == i-1) A[i][i-1] = -birth * (2*k + i - 1);
				}
			}

			//A = A * (th - thMinusOne);
			for(int i = 0; i < (N + 1); i++){
				for(int j = 0; j < (N + 1); j++){
					A[i][j] = A[i][j] * (th - thPlusOne);
				}
			}

		 RbMath::expMatrixPade(A, A, 4);

		 // A * Mt
		 RbVector<double> MtPrime;
		 // for each row in matrix A
		 for(int i = 0; i < Mt.size(); i++){
			 double sum = 0.0;
			 // for each entry in the vector Mt
			 for(int j = 0; j < Mt.size(); j++){
				 sum += A[i][j] * Mt[j];
			 }
			 MtPrime.push_back( sum );
		 }

		 //Mt = MtPrime;
		 for(int i = 0; i < Mt.size(); i++){
       std::cout << " MtPrime" <<std::endl ;
       std::cout << MtPrime <<std::endl ; 
				Mt[i] = MtPrime[i];
		 }

		}

		std::string type = events[h].type;

		// step 7-10. sample Mt
		if(type == "time slice"){
			for(int i = 0; i < Mt.size(); i++){
				B[indxJ][i] = Mt[i];
			}
			indxJ -= 1;
		}

    // step 13-14.
		if(type == "terminal removed"){
      std::cout << "terminal removed" << std::endl;
			for(int i = 0; i < Mt.size(); i++){
				Mt[i] = Mt[i] * ps * rp;
        std::cout << Mt[i] << std::endl;
			}
			k -= 1;
		}

		// step 15-16.
		if(type == "terminal non-removed"){
      std::cout << "terminal non-removed" << std::endl;

			for(int i = Mt.size()-1; i > 0; i--){
				Mt[i] = Mt[i-1] * ps * (1-rp);
      std::cout << Mt[i] << std::endl;
			}
			Mt[0] = 0;
			k -= 1;
		}

		// step 17-18.
		if(type == "sampled ancestor"){
      std::cout << "sampled ancestor" << std::endl;

			for(int i = 0; i < Mt.size(); i++){
				Mt[i] = Mt[i] * ps * (1-rp);
        std::cout << Mt[i] << std::endl;
			}
		}

		// step 19-20.
		if(type == "occurrence removed"){
      std::cout << "occurence removed" << std::endl;
			for(int i = 0; i < Mt.size()-1; i++){
				Mt[i] = Mt[i+1] * (i+1) * om * rp;

			}
		 }

		// step 21-22.
		if(type == "occurrence non-removed"){
      std::cout << "occurence non removed" << std::endl;

			for(int i = 0; i < Mt.size(); i++){
				Mt[i] = Mt[i] * (k+i) * om * (1-rp);
			}
		}

		// step 23-24.
		if(type == "branching time"){
      std::cout << "branching time" << std::endl;
			for(int i = 0; i < Mt.size(); i++){
				Mt[i] = Mt[i] * birth;
        std::cout << Mt[i] << std::endl;
			}
			k += 1;
		}

		thPlusOne = th;

	}

  std::cout<< "Mt[0]" << std::endl;
  std::cout<< Mt[0] << std::endl;
  std::cout<< "likelihoods" << std::endl;
  double test = 10e-50 ;
  std::cout << test << std::endl;
  double likelihood = Mt[0];
  for(int i = 1; i < Mt.size(); i++){
    std::cout << likelihood << std::endl;
    std::cout<< "Mt[i]s" << std::endl;
    std::cout << Mt[i] * pow(rh,k) * pow(1-rh,i) << std::endl;
    likelihood += Mt[i] * pow(rh,k) * pow(1-rh,i);




  }

  std::cout << "Compute mt output" << std::endl;
  std::cout << log(likelihood) << std::endl;
	return log(likelihood);

}





double OccurrenceBirthDeathProcess::ComputeLt( void ) const
{
// order times youngest to oldest
  std::cout << "Sorting events" << std::endl;
	std::sort( events.begin(), events.end(), AgeCompare() );

	const double birth = lambda->getValue();
	const double death = mu->getValue();
	const double ps = psi->getValue();
	const double om = omega->getValue();
	const double rh = rho->getValue();
	const double rp = removalPr->getValue();

	const double gamma = birth + death + ps + om;

	size_t N = 10; // accuracy of the algorithm
	size_t S = dn_time_points->getValue().size();
  std::cout << "S setting time points for comput" << std::endl;
  std::cout << S << std::endl;


	// step 2. initialize an empty matrix
	MatrixReal B(S, (N + 1), 0.0);

	size_t k = extant;

	// step 3. // this is the first entry of B
	RbVector<double> Lt;
	for ( int i = 0; i < (N + 1); i++){
		if (i == 0) Lt.push_back( pow(rh, k) );
		else Lt.push_back( pow(rh, k) * pow((1-rh), i) );
	}
  //std::cout << Lt << std::endl;

	// step 4.
	//-> calculate the state of the process just before the next time point

	MatrixReal A( (N+1), (N+1), 0.0 );

	double thMinusOne = 0.0;

	size_t indxJ = 0;

	for(int h = 0; h < events.size(); h++){

		double th = events[h].time;

		// if time is not the same recalculate A.
		// step 5-6.
		if( th != thMinusOne ){

			for(int i = 0; i < (N + 1); i++){
				for(int j = 0; j < (N + 1); j++){
					if(j == i) A[i][i] = -gamma * (k + i);
					else if (j == i+1) A[i][i+1] = birth * ( (2 * k) + i );
					else if (j == i-1) A[i][i-1] = death * i;
				}
			}

			//A = A * (th - thMinusOne);
			for(int i = 0; i < (N + 1); i++){
				for(int j = 0; j < (N + 1); j++){
					A[i][j] = A[i][j] * (th - thMinusOne);
				}
			}

		 RbMath::expMatrixPade(A, A, 4);

		 // A * Lt
		 RbVector<double> LtPrime;
		 // for each row in matrix A
		 for(int i = 0; i < Lt.size(); i++){
			 double sum = 0.0;
			 // for each entry in the vector Lt
			 for(int j = 0; j < Lt.size(); j++){
				 sum += A[i][j] * Lt[j];
			 }
			 LtPrime.push_back( sum );
		 }

		 //Lt = LtPrime;
		 for(int i = 0; i < Lt.size(); i++){
				Lt[i] = LtPrime[i];
		 }

		}

		std::string type = events[h].type;

		// step 7-10. sample Lt
		if(type == "time slice"){
			for(int i = 0; i < Lt.size(); i++){
				B[indxJ][i] = Lt[i];
			}
			indxJ += 1;
		}

		// step 11-12. end algorithm 1
		if(type == "origin"){
			continue;
		}

		// step 13-14.
		if(type == "terminal removed"){
			for(int i = 0; i < Lt.size(); i++){
				Lt[i] = Lt[i] * ps * rp;
        std::cout << Lt[i] << std::endl;
			}
			k += 1;
		}

		// step 15-16.
		if(type == "terminal non-removed"){
			for(int i = 0; i < Lt.size()-1; i++){
				Lt[i] = Lt[i+1] * ps * (1-rp);
        std::cout << Lt[i] << std::endl;

			}
			//Lt[N] = 0;
			k += 1;
		}

		// step 17-18.
		if(type == "sampled ancestor"){
			for(int i = 0; i < Lt.size(); i++){
				Lt[i] = Lt[i] * ps * (1-rp);
        std::cout << Lt[i] << std::endl;

			}
		}

		// step 19-20.
		if(type == "occurrence removed"){
		 	for(int i = Lt.size()-1; i > 0; i--){
		 		Lt[i] = Lt[i-1] * i * om * rp;
        std::cout << Lt[i] << std::endl;

		 	}
		 	Lt[0] = 0;
		 }

		// step 21-22.
		if(type == "occurrence non-removed"){
			for(int i = 0; i < Lt.size(); i++){
				Lt[i] = Lt[i] * (k+i) * om * (1-rp);
        std::cout << Lt[i] << std::endl;

			}
		}

		// step 23-24.
		if(type == "branching time"){
			for(int i = 0; i < Lt.size(); i++){
				Lt[i] = Lt[i] * birth;
        std::cout << Lt[i] << std::endl;

			}
			k -= 1;
		}

		thMinusOne = th;

	}
  std::cout << "Compute Lt output" << std::endl;
  std::cout << Lt[0] << std::endl;
	return log(Lt[0]);

}

//
//  ComputeLtFunction.cpp
//
//  Created by Rachel Warnock, Marc Manceau 30.01.2020.
//
#include "ComputeLtFunction.h"

#include <vector>
#include <iostream>
#include <cmath>

#include "RbMathMatrix.h"
#include "RbVector.h"
#include "TypedDagNode.h"
#include "Tree.h"
#include "TopologyNode.h"
// other includes

namespace RevBayesCore { class DagNode; } //tt

using namespace RevBayesCore;

/** ComputeLtFunction
 * this is temporary until we work out where this function should live in the long term
 * @param v the vector of branching times
 */
ComputeLtFunction::ComputeLtFunction(
	const TypedDagNode< RbVector<double> > *a, // branching times
	const TypedDagNode< RbVector<double> > *b, // SAs
	const TypedDagNode< RbVector<double> > *c, // terminal non-removed
	const TypedDagNode< RbVector<double> > *d, // terminal removed
	const TypedDagNode< RbVector<double> > *e, // occ non-removed
	const TypedDagNode< RbVector<double> > *f, // occ removed
	const TypedDagNode< double > *t,
	const TypedDagNode< double > *l,
	const TypedDagNode< double > *m,
	const TypedDagNode< double > *p,
	const TypedDagNode< double > *o,
	const TypedDagNode< double > *rho,
	const TypedDagNode< double > *r,
	const TypedDagNode< RbVector<double> > *g // times
 	) : TypedFunction<double>( new double(0.0) ),
	listA ( a ), listB ( b ), listC ( c ), listD ( d ), listE ( e ), listF ( f ),
	tor ( t ), lambda ( l ), mu ( m ), psi ( p ), omega ( o ), rho ( rho ), removalPr ( r ),
	listG ( g )
{

	this->addParameter( tor );
	this->addParameter( lambda );
	this->addParameter( mu );
	this->addParameter( psi );
	this->addParameter( omega );
	this->addParameter( rho );
	this->addParameter( removalPr ); // I'm not sure if we need to do this

	useTree = false;

	poolTimes(); // step 1.
	//update();

}

ComputeLtFunction::ComputeLtFunction(
	const TypedDagNode<Tree> *tr,
	const TypedDagNode< RbVector<double> > *b,
	const TypedDagNode< RbVector<double> > *c,
	const TypedDagNode< RbVector<double> > *d, // terminal removed
	const TypedDagNode< RbVector<double> > *e, // occ non-removed
	const TypedDagNode< RbVector<double> > *f, // occ removed
	const TypedDagNode< double > *t,
	const TypedDagNode< double > *l,
	const TypedDagNode< double > *m,
	const TypedDagNode< double > *p,
	const TypedDagNode< double > *o,
	const TypedDagNode< double > *rho,
	const TypedDagNode< double > *r,
	const TypedDagNode< RbVector<double> > *g // times
 	) : TypedFunction<double>( new double(0.0) ),
	tree ( tr ), listB ( b ), listC ( c ), listD ( d ), listE ( e ), listF ( f ),
	tor ( t ), lambda ( l ), mu ( m ), psi ( p ), omega ( o ), rho ( rho ), removalPr ( r ),
	listG ( g )
{

	this->addParameter( tor );
	this->addParameter( lambda );
	this->addParameter( mu );
	this->addParameter( psi );
	this->addParameter( omega );
	this->addParameter( rho );
	this->addParameter( removalPr );

	// see also Tree class -> addNodeParameter
	// getTipNodeWithName
	// renameNodeParameter
	// getNumberOfExtantTips
	// getNumberOfSampledAncestors

	useTree = true;

	size_t num_nodes = tree->getValue().getNumberOfNodes();

	for (size_t i = 0; i < num_nodes; i++){

      const TopologyNode& n = tree->getValue().getNode( i );

			double t = n.getAge();

			if ( n.isFossil() && n.isSampledAncestor() ){
				// 1. sampled ancestor
				SA_ages.push_back(t);
			} else if ( n.isFossil() && !n.isSampledAncestor() ){
				// 2. node is serial leaf
				tip_ages.push_back(t);
			} // else if ( n.isTip() && !n.isFossil() ){
				// 3. node is at present
				// tip_ages.push_back(t);
			// }
			else if ( n.isInternal() && !n.getChild(0).isSampledAncestor() && !n.getChild(1).isSampledAncestor() ){
				// 4. true branching event
				branching_times.push_back(t);
			}
	}

	for(int i = 0; i < branching_times.size(); i++){
		std::cout << "bt " << i << " : " << branching_times[i] << std::endl;
	}
	for(int i = 0; i < SA_ages.size(); i++){
		std::cout << "SA " << i << " : " << SA_ages[i] << std::endl;
	}
	for(int i = 0; i < tip_ages.size(); i++){
		std::cout << "tr " << i << " : " << tip_ages[i] << std::endl;
	}

	poolTimes(); // step 1.
	//update();

}

ComputeLtFunction::~ComputeLtFunction( void ){
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}

/** clone function */
ComputeLtFunction* ComputeLtFunction::clone( void ) const
{
    return new ComputeLtFunction( *this );
}

// required -> update function, where you do the work
void ComputeLtFunction::update( void ) {

	double a1 = ComputeLt();
	double a2 = ComputeMt();

	double out = a2;

	*this->value = out; // this will eventually be the lk of the tree + occurrences

}

double ComputeLtFunction::ComputeMt( void ) {

	// order times oldest to youngest -> revisit this
	std::sort( events.begin(), events.end(), AgeCompareReverse() );

	const double birth = lambda->getValue();
	const double death = mu->getValue();
	const double ps = psi->getValue();
	const double om = omega->getValue();
	const double rh = rho->getValue();
	const double rp = removalPr->getValue();

	const double gamma = birth + death + ps + om;

	size_t N = 4; // accuracy of the algorithm
	size_t S = listG->getValue().size();

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
			for(int i = 0; i < Mt.size(); i++){
				Mt[i] = Mt[i] * ps * rp;
			}
			k -= 1;
		}

		// step 15-16.
		if(type == "terminal non-removed"){
			for(int i = Mt.size()-1; i > 0; i--){
				Mt[i] = Mt[i-1] * ps * (1-rp);
			}
			Mt[0] = 0;
			k -= 1;
		}

		// step 17-18.
		if(type == "sampled ancestor"){
			for(int i = 0; i < Mt.size(); i++){
				Mt[i] = Mt[i] * ps * (1-rp);
			}
		}

		// step 19-20.
		if(type == "occurrence removed"){
			for(int i = 0; i < Mt.size()-1; i++){
				Mt[i] = Mt[i+1] * (i+1) * om * rp;
			}
		 }

		// step 21-22.
		if(type == "occurrence non-removed"){
			for(int i = 0; i < Mt.size(); i++){
				Mt[i] = Mt[i] * (k+i) * om * (1-rp);
			}
		}

		// step 23-24.
		if(type == "branching time"){
			for(int i = 0; i < Mt.size(); i++){
				Mt[i] = Mt[i] * birth;
			}
			k += 1;
		}

		thPlusOne = th;

	}

  double likelihood = 0.0;
  for(int i = 0; i < Mt.size(); i++){
    likelihood += Mt[i] * pow(rh,k) * pow(1-rh,i);
  }

	return likelihood;

}

double ComputeLtFunction::ComputeLt( void ) {

	// order times youngest to oldest
	std::sort( events.begin(), events.end(), AgeCompare() );

	const double birth = lambda->getValue();
	const double death = mu->getValue();
	const double ps = psi->getValue();
	const double om = omega->getValue();
	const double rh = rho->getValue();
	const double rp = removalPr->getValue();

	const double gamma = birth + death + ps + om;

	size_t N = 4; // accuracy of the algorithm
	size_t S = listG->getValue().size();

	// step 2. initialize an empty matrix
	MatrixReal B(S, (N + 1), 0.0);

	size_t k = extant;

	// step 3. // this is the first entry of B
	RbVector<double> Lt;
	for ( int i = 0; i < (N + 1); i++){
		if (i == 0) Lt.push_back( pow(rh, k) );
		else Lt.push_back( pow(rh, k) * pow((1-rh), i) );
	}

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
			}
			k += 1;
		}

		// step 15-16.
		if(type == "terminal non-removed"){
			for(int i = 0; i < Lt.size()-1; i++){
				Lt[i] = Lt[i+1] * ps * (1-rp);
			}
			//Lt[N] = 0;
			k += 1;
		}

		// step 17-18.
		if(type == "sampled ancestor"){
			for(int i = 0; i < Lt.size(); i++){
				Lt[i] = Lt[i] * ps * (1-rp);
			}
		}

		// step 19-20.
		if(type == "occurrence removed"){
		 	for(int i = Lt.size()-1; i > 0; i--){
		 		Lt[i] = Lt[i-1] * i * om * rp;
		 	}
		 	Lt[0] = 0;
		 }

		// step 21-22.
		if(type == "occurrence non-removed"){
			for(int i = 0; i < Lt.size(); i++){
				Lt[i] = Lt[i] * (k+i) * om * (1-rp);
			}
		}

		// step 23-24.
		if(type == "branching time"){
			for(int i = 0; i < Lt.size(); i++){
				Lt[i] = Lt[i] * birth;
			}
			k -= 1;
		}

		thMinusOne = th;

	}

	return Lt[0];

}

// required -> swap function. See Variance Function for inspiration
void ComputeLtFunction::swapParameterInternal( const DagNode *oldP, const DagNode *newP ) {

	if ( oldP == listB	)
	{
			listB = static_cast<const TypedDagNode< RbVector<double> >* >( newP ); //tt
	}

}

void ComputeLtFunction::poolTimes(void) {

	//const std::vector<double> &a = listA->getValue(); // branching times from a vector
	//const std::vector<double> &b = listB->getValue(); // sampled ancestors
	//const std::vector<double> &c = listC->getValue(); // terminal non-removed

	// we should be able to extract listD ages from the tree too
	const std::vector<double> &d = listD->getValue(); // terminal removed

	const std::vector<double> &e = listE->getValue(); // occurrences non-removed
	const std::vector<double> &f = listF->getValue(); // occurrences removed
	const std::vector<double> &g = listG->getValue(); // timeslices

	// this seems overly complex
	if(useTree){

		extant = branching_times.size() + 1 - tip_ages.size() - d.size();

		for(int i = 0; i < branching_times.size(); i++){
		    events.push_back( Event(branching_times[i], "branching time", 0) );
		}

		for(int i = 0; i < SA_ages.size(); i++){
		    events.push_back( Event(SA_ages[i], "branching time", 0) );
		}

		for(int i = 0; i < tip_ages.size(); i++){
		    events.push_back( Event(tip_ages[i], "terminal non-removed", 0) );
		}

	} else {

		extant = listA->getValue().size() + 1 - listC->getValue().size() - d.size();

		for(int i = 0; i < listA->getValue().size(); i++){
	    events.push_back( Event(listA->getValue()[i], "branching time", 0) );
		}

		for(int i = 0; i < listB->getValue().size(); i++){
				events.push_back( Event(listB->getValue()[i], "sampled ancestor", 0) );
		}

		for(int i = 0; i < listC->getValue().size(); i++){
		    events.push_back( Event(listC->getValue()[i], "terminal non-removed", 0) );
		}

	}

	for(int i = 0; i < d.size(); i++){
	    events.push_back( Event(d[i], "terminal removed", 0) );
	}

	for(int i = 0; i < e.size(); i++){
	    events.push_back( Event(e[i], "occurrence non-removed", 0) );
	}

	for(int i = 0; i < f.size(); i++){
	    events.push_back( Event(f[i], "occurrence removed", 0) );
	}

	for(int i = 0; i < g.size(); i++){
	    events.push_back( Event(g[i], "time slice", 0) );
	}

	events.push_back( Event(tor->getValue(), "origin", 0) );

}

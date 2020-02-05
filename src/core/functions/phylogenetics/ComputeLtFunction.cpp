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
// other includes

namespace RevBayesCore { class DagNode; } //tt

using namespace RevBayesCore;

/** ComputeLtFunction
 * this is temporary until we work out where this function should live in the long term
 * @param v the vector of branching times
 */
ComputeLtFunction::ComputeLtFunction(
	const TypedDagNode< RbVector<double> > *a,
	const TypedDagNode< RbVector<double> > *b,
	const TypedDagNode< RbVector<double> > *c,
	const TypedDagNode< RbVector<double> > *d,
	const TypedDagNode< RbVector<double> > *e,
	const TypedDagNode< RbVector<double> > *f,
	const TypedDagNode< double > *t,
	const TypedDagNode< double > *l,
	const TypedDagNode< double > *m,
	const TypedDagNode< double > *p,
	const TypedDagNode< double > *o,
	const TypedDagNode< double > *rho,
	const TypedDagNode< double > *r,
	const TypedDagNode< RbVector<double> > *g
 	) : TypedFunction<double>( new double(0.0) ),
	listA ( a ), listB ( b ), listC ( c ), listD ( d ), listE ( e ), listF ( f ),
	tor ( t ), lambda ( l ), mu ( m ), psi ( p ), omega ( o ), rho ( rho ), removalPr ( r ),
	listG ( g )
{

	//poolTimes(a, b, b, c, d, e, f, g);

	//add the parameter as parents
	//this->addParameter( listA );
	//this->addParameter( listB );
	//this->addParameter( listC );
	//this->addParameter( listD );
	//this->addParameter( listE );
	//this->addParameter( listF );

	this->addParameter( tor );
	this->addParameter( lambda );
	this->addParameter( mu );
	this->addParameter( psi );
	this->addParameter( omega );
	this->addParameter( rho );
	this->addParameter( removalPr ); // I'm not sure if we need to do this

	poolTimes(); // step 1.
	update();

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

	const double birth = lambda->getValue();
	const double death = mu->getValue();
	const double ps = psi->getValue();
	const double om = omega->getValue();
	const double rh = rho->getValue();
	const double rp = removalPr->getValue();

	const double gamma = birth + death + ps + om;

	// debugging
	//std::cout << "gamma " << gamma << std::endl;

	size_t N = 4; // accuracy of the algorithm
	size_t S = listG->getValue().size();

	// step 2. initialize an empty matrix
	MatrixReal B(S, (N + 1), 0.0);

	// debugging
	std::cout << "B 1" << B << std::endl;

	size_t k = extant;

	// debugging
	//std::cout << "k " << k << std::endl;

	// step 3. // I think this is supposed to be the first entry of B
	//std::vector<double> Lt;
	RbVector<double> Lt;
	for ( int i = 0; i < (N + 1); i++){
		if (i == 0) Lt.push_back( pow(rh, k) );
		else Lt.push_back( pow(rh, k) * pow((1-rh), i) );
	}

	// debugging
	std::cout << "Lt 1" << Lt << std::endl;

	// step 4.
	//-> calculate the state of the process just before the next time point

	MatrixReal A( (N+1), (N+1), 0.0 );

	double thMinusOne = 0.0;

	size_t indxJ = 0;

	//for(int h = 0; h < S; h++){
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

			//std::cout << "A 1" << A << std::endl;

			//A = A * (th - thMinusOne);
			for(int i = 0; i < (N + 1); i++){
				for(int j = 0; j < (N + 1); j++){
					A[i][j] = A[i][j] * (th - thMinusOne);
				}
			}

			//std::cout << std::setw(6) << "A 2" << A << std::endl;

			RbMath::expMatrixPade(A, A, 4);

			//std::cout << "A 3" << A << std::endl;

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
			//std::cout << "B 2" << B << std::endl;
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
			//debug
			for(int i = 0; i < Lt.size(); i++){
				std::cout << "step 13 " << Lt[i] << std::endl;
			}
		}

		// step 15-16.
		if(type == "terminal non-removed"){
			for(int i = 0; i < Lt.size()-1; i++){
				Lt[i] = Lt[i+1] * ps * (1-rp);
			}
			//Lt[N] = 0;
			k += 1;
			//debug
			for(int i = 0; i < Lt.size(); i++){
				std::cout << "step 15 " << Lt[i] << std::endl;
			}
		}

		// step 17-18.
		if(type == "sampled ancestor"){
			for(int i = 0; i < Lt.size(); i++){
				Lt[i] = Lt[i] * ps * (1-rp);
			}
			//debug
			for(int i = 0; i < Lt.size(); i++){
				std::cout << "step 17 " << Lt[i] << std::endl;
			}
		}

		// step 19-20.
		if(type == "occurrence removed"){
		 	for(int i = Lt.size()-1; i > 0; i--){
		 		Lt[i] = Lt[i-1] * i * om * rp;
		 	}
		 	Lt[0] = 0;
		 	//debug
		 	for(int i = 0; i < Lt.size(); i++){
		 		std::cout << "step 19 " << Lt[i] << std::endl;
		 	}
		 }

		// step 21-22.
		if(type == "occurrence non-removed"){
			for(int i = 0; i < Lt.size(); i++){
				Lt[i] = Lt[i] * (k+i) * om * (1-rp);
			}
			//debug
			for(int i = 0; i < Lt.size(); i++){
				std::cout << "step 21 " << Lt[i] << std::endl;
			}
		}

		// step 23-24.
		if(type == "branching time"){
			for(int i = 0; i < Lt.size(); i++){
				Lt[i] = Lt[i] * birth;
			}
			k -= 1;
			//debug
			for(int i = 0; i < Lt.size(); i++){
				std::cout << "step 23 " << Lt[i] << std::endl;
			}
		}

		thMinusOne = th;

	}

	//double out = log(Lt[0]);

	double out = Lt[0];

	*this->value = out; // this will eventually be the lk of the tree + occurrences

}

// required -> swap function. See Variance Function for inspiration
void ComputeLtFunction::swapParameterInternal( const DagNode *oldP, const DagNode *newP ) {

	if ( oldP == listA	)
	{
			listA = static_cast<const TypedDagNode< RbVector<double> >* >( newP ); //tt
	}

}

void ComputeLtFunction::poolTimes(void) {

	const std::vector<double> &a = listA->getValue(); // branching times
	const std::vector<double> &b = listB->getValue(); // sampled ancestors
	const std::vector<double> &c = listC->getValue(); // terminal non-removed
	const std::vector<double> &d = listD->getValue(); // terminal removed
	const std::vector<double> &e = listE->getValue(); // occurrences non-removed
	const std::vector<double> &f = listF->getValue(); // occurrences removed
	const std::vector<double> &g = listG->getValue(); // timeslices

	extant = a.size() + 1 - c.size() - d.size();

	for(int i = 0; i < a.size(); i++){
	    events.push_back( Event(a[i], "branching time", 0) );
	}

	for(int i = 0; i < b.size(); i++){
	    events.push_back( Event(b[i], "sampled ancestor", 0) );
	}

	for(int i = 0; i < c.size(); i++){
	    events.push_back( Event(c[i], "terminal non-removed", 0) );
	}

	for(int i = 0; i < d.size(); i++){
	    events.push_back( Event(d[i], "terminal removed", 0) );
	}

	for(int i = 0; i < e.size(); i++){
	    events.push_back( Event(e[i], "terminal non-removed", 0) );
	}

	for(int i = 0; i < f.size(); i++){
	    events.push_back( Event(f[i], "occurrence removed", 0) );
	}

	for(int i = 0; i < g.size(); i++){
	    events.push_back( Event(g[i], "time slice", 0) );
	}

	events.push_back( Event(tor->getValue(), "origin", 0) );

	// order times youngest to oldest
	std::sort( events.begin(), events.end(), AgeCompare() );

	// debugging code
	for(int i = 0; i < events.size(); i++){
	 	Event e = events[i];
	 	//std::cout << i << " " << e.time << std::endl;
	}

}

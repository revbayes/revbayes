#include "EigenSystem.h"
#include "MatrixReal.h"
#include "RateMatrix_MPQ.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RandomNumberGenerator.h"
#include "RandomNumberFactory.h"
#include "DistributionBeta.h"
#include "DistributionDirichlet.h"
#include "RbMathVector.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>

#define MIN_FREQ    10e-8

/* Smallest component of a Dirichlet draw we are willing to build a reverse
   concentration parameter from; see isDrawUsable below. */
#define MIN_DIRICHLET_DRAW  1.0e-290

/* Draws needed to land inside the polyhedron when initializing the non-reversible
   model. About one draw in five is admissible, so this is a very loose bound whose
   only purpose is to fail loudly instead of hanging. */
#define MAX_INITIALIZATION_ATTEMPTS  100000
#define A           0
#define C           1
#define G           2
#define T           3

// indices into the six backbone weights returned by calculateWeights()
#define W_AC        0
#define W_AG        1
#define W_AT        2
#define W_CG        3
#define W_CT        4
#define W_GT        5


/* ---------------------------------------------------------------------------
   Weight coordinates
   ------------------

   The state of this rate matrix is the pair (pi, w), where w_ij = pi_i q_ij is
   the average rate of change from nucleotide i to nucleotide j.  In these
   coordinates every constraint of the model is linear and, crucially, free of
   pi:

       sum_{i != j} w_ij  =  1                        (average rate is one)
       sum_{i != j} w_ij  =  sum_{k != j} w_jk        (pi is stationary)
       w_ij >= 0                                      (validity)

   so the non-reversible model is the polytope of non-negative unit circulations
   on the complete digraph K4 (dimension 8), and the time-reversible model is
   its "no net circulation" face w_ij = w_ji (dimension 5).  Neither polytope
   moves when pi moves.  The prior is therefore flat on each polytope, with
   density 1/V_R = 3840 and 1/V_N = 1935360 respectively, and the reversible
   jump Jacobian J_RN = 64 w_CG w_CT w_GT is taken in the same free coordinates.

   Every move below is written as a move on (pi, w).  This matters: a move
   written in terms of the exchangeability rates r, or of the raw q_ij, changes
   w through a non-linear map and needs the corresponding reparameterization
   Jacobian.  Omitting it does not show up in the acceptance rate, only as a
   quiet skew in the sampled weights.  Working directly in w removes the whole
   class of error, and removes every power of pi from the prior and from the
   jump Jacobian.
   --------------------------------------------------------------------------- */


using namespace RevBayesCore;

RateMatrix_MPQ::RateMatrix_MPQ(void) : RateMatrix(4) {

    q = new mpq_class[16];
    endBuffer = q + 16;
    
    the_rate_matrix = new MatrixReal(4);
    theEigenSystem  = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    pi.resize(4);
    r.resize(6);
    isReversible = false;

    /* needs_update was never initialized. It was also never set and never read,
       so the eigensystem was recomputed on every single call to
       calculateTransitionProbabilities: once per branch per site-rate category
       per likelihood evaluation, where once per accepted move would do. */
    needs_update       = true;
    eigen_update_count = 0;
}

RateMatrix_MPQ::RateMatrix_MPQ(const RateMatrix_MPQ& m) : RateMatrix(m) {
    
    this->isReversible = m.isReversible;
    
    q = new mpq_class[16];
    endBuffer = q + 16;
    
    the_rate_matrix = new MatrixReal(4);
    theEigenSystem       = new EigenSystem( *m.theEigenSystem );
    c_ijk                = m.c_ijk;
    cc_ijk               = m.cc_ijk;
    theEigenSystem->setRateMatrixPtr(the_rate_matrix);
    
    pi.resize(4);
    r.resize(6);
    
    mpq_class* r = q;
    for (mpq_class* p=m.q; p != m.endBuffer; p++)
        {
        *r = *p;
        r++;
        }
    for (int i=0; i<4; i++)
        this->pi[i] = m.pi[i];
    for (int i=0; i<6; i++)
        this->r[i] = m.r[i];
    
    moveToDouble();

    /* Take the source's staleness, not a blanket "stale". We copied its c_ijk, so
       if the source's cache was valid ours is too. This is what makes a rejected
       proposal free: undoProposal assigns the stored matrix back, and the stored
       matrix was clean, so no eigendecomposition is needed to resume. */
    needs_update       = m.needs_update;
    eigen_update_count = 0;
}

RateMatrix_MPQ::~RateMatrix_MPQ(void) {

    delete [] q;
    delete the_rate_matrix;
    delete theEigenSystem;
}

RateMatrix_MPQ& RateMatrix_MPQ::operator=(const RateMatrix_MPQ& rhs) {

    if (this != &rhs)
        {
        this->isReversible = rhs.isReversible;
            
        delete theEigenSystem;
            
        theEigenSystem       = new EigenSystem( *rhs.theEigenSystem );
        c_ijk                = rhs.c_ijk;
        cc_ijk               = rhs.cc_ijk;
            
        theEigenSystem->setRateMatrixPtr(the_rate_matrix);
            
        mpq_class* r = q;
        for (mpq_class* p=rhs.q; p != rhs.endBuffer; p++)
            {
            *r = *p;
            r++;
            }
        for (int i=0; i<4; i++)
            this->pi[i] = rhs.pi[i];
        for (int i=0; i<6; i++)
            this->r[i] = rhs.r[i];
            
        moveToDouble();

        // as in the copy constructor: c_ijk came from rhs, so rhs's staleness is ours
        needs_update = rhs.needs_update;
        }
    return *this;
}

void RateMatrix_MPQ::adjust(void) {

    if (isReversible == true)
        {
        std::vector<double> adjR(6);
        for (int i=0; i<6; i++)
            adjR[i] = r[i].get_d();
        for (int i=0; i<5; i++)
            r[i] = adjR[i];
        r[5] = 1 - (r[0] + r[1] + r[2] + r[3] + r[4]);
        std::vector<double> adjPi(4);
        for (int i=0; i<4; i++)
            adjPi[i] = pi[i].get_d();
        for (int i=0; i<3; i++)
            pi[i] = adjPi[i];
        pi[3] = 1 - (pi[0] + pi[1] + pi[2]);
            
        for (int i=0, k=0; i<4; i++)
            {
            for (int j=i+1; j<4; j++)
                {
                (*this)(i,j) = r[k] * pi[j];
                (*this)(j,i) = r[k] * pi[i];
                k++;
                }
            }
        mpq_class averageRate;
        for (int i=0; i<4; i++)
            {
            mpq_class sum;
            for (int j=0; j<4; j++)
                {
                if (i != j)
                    sum += (*this)(i,j);
                }
            (*this)(i,i) = -sum;
            averageRate += pi[i] * sum;
            }
        mpq_class factor = 1 / averageRate;
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                (*this)(i,j) *= factor;
                
        }
    else
        {
        mpq_class adjQ[4][4];
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                adjQ[i][j] = (*this)(i,j);
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                (*this)(i,j) = adjQ[i][j];
                
        calculateStationaryFrequencies(this->pi);
                
        mpq_class averageRate;
        for (int i=0; i<4; i++)
            {
            mpq_class sum;
            for (int j=0; j<4; j++)
                {
                if (i != j)
                    sum += (*this)(i,j);
                }
            (*this)(i,i) = -sum;
            averageRate += pi[i] * sum;
            }
        mpq_class factor = 1 / averageRate;
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                (*this)(i,j) *= factor;
        }
}

double RateMatrix_MPQ::averageRate(void) const {

    mpq_class tmp;
    calculateAverageRate( tmp );
    return tmp.get_d();
}                                                                //!< Calculate the average rate

void RateMatrix_MPQ::calculateAverageRate(mpq_class& ave) const {

    ave = 0;
    for (int i=0; i<4; i++)
        ave += -(pi[i] * (*this)(i,i));
}

/** Do precalculations on eigenvectors */
void RateMatrix_MPQ::calculateCijk(void) {
    
    if ( theEigenSystem->isComplex() == false )
        {
        // real case
        const MatrixReal& ev  = theEigenSystem->getEigenvectors();
        const MatrixReal& iev = theEigenSystem->getInverseEigenvectors();
        double* pc = &c_ijk[0];
        for (size_t i=0; i<num_states; i++)
            {
            for (size_t j=0; j<num_states; j++)
                {
                for (size_t k=0; k<num_states; k++)
                    {
                    *(pc++) = ev[i][k] * iev[k][j];
                    }
                }
            }
        }
    else
        {
        // complex case
        const MatrixComplex& cev  = theEigenSystem->getComplexEigenvectors();
        const MatrixComplex& ciev = theEigenSystem->getComplexInverseEigenvectors();
        std::complex<double>* pc = &cc_ijk[0];
        for (size_t i=0; i<num_states; i++)
            {
            for (size_t j=0; j<num_states; j++)
                {
                for (size_t k=0; k<num_states; k++)
                    {
                    *(pc++) = cev[i][k] * ciev[k][j];
                    }
                }
            }
        }
}

void RateMatrix_MPQ::calculateStationaryFrequencies(std::vector<mpq_class>& f) {

    // transpose the rate matrix (qMatrix) and put into QT
    RateMatrix_MPQ QT;
    transposeMatrix(*this, QT);

    // compute the LU decomposition of the transposed rate matrix
    RateMatrix_MPQ L;
    RateMatrix_MPQ U;
    computeLandU(QT, L, U);
    
    // back substitute into z = 0 to find un-normalized stationary frequencies
    // start with x_n = 1.0
    f[3] = 1;
    for (int i=4-2; i>=0; i--)
        {
        mpq_class dotProduct;
        for (int j=i+1; j<4; j++)
            dotProduct += U(i,j) * f[j];
        f[i] = (0 - dotProduct) / U(i,i);
        }
        
    // normalize the solution vector
    mpq_class sum;
    for (int i=0; i<4; i++)
        sum += f[i];
    for (int i=0; i<4; i++)
        f[i] /= sum;
  
    // make certain to initialize the instance variable
    for (int i=0; i<4; i++)
        this->pi[i] = f[i];
}

void RateMatrix_MPQ::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const {
    
    // recompute the eigensystem only if the rate matrix has changed since the
    // last time we did
    updateIfNeeded();
    
    double t = rate * (startAge - endAge);
    if ( theEigenSystem->isComplex() == false )
        {
        tiProbsEigens(t, P);
        }
    else
        {
        tiProbsComplexEigens(t, P);
        }
}

/* All twelve weights w_ij = pi_i q_ij, laid out as a 4 X 4 array with zeros on
   the diagonal. */
void RateMatrix_MPQ::calculateAllWeights(std::vector<mpq_class>& w) const {

    if (w.size() != 16)
        w.resize(16);
    for (int i=0; i<4; i++)
        {
        for (int j=0; j<4; j++)
            {
            if (i == j)
                w[i*4+j] = 0;
            else
                w[i*4+j] = pi[i] * (*this)(i,j);
            }
        }
}

/* Install a full set of twelve weights: q_ij = w_ij / pi_i, with the diagonal
   set so that each row sums to zero.  Returns false if any off-diagonal weight
   is negative, which puts the point outside the polytope and so outside the
   support of the prior. */
bool RateMatrix_MPQ::setRatesFromAllWeights(const std::vector<mpq_class>& w) {

    for (int i=0; i<4; i++)
        {
        for (int j=0; j<4; j++)
            {
            if (i != j && w[i*4+j] < 0)
                return false;
            }
        if (pi[i] <= 0)
            return false;
        }

    for (int i=0; i<4; i++)
        {
        mpq_class sum = 0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                {
                (*this)(i,j) = w[i*4+j] / pi[i];
                sum += (*this)(i,j);
                }
            }
        (*this)(i,i) = -sum;
        }
    return true;
}

/* Install a time-reversible rate matrix from the six backbone weights, which
   are assumed to sum to one half. */
void RateMatrix_MPQ::setReversibleRatesFromBackbone(const std::vector<mpq_class>& wR) {

    for (int i=0, k=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            {
            (*this)(i,j) = wR[k] / pi[i];
            (*this)(j,i) = wR[k] / pi[j];
            k++;
            }
        }

    for (int i=0; i<4; i++)
        {
        mpq_class sum = 0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += (*this)(i,j);
            }
        (*this)(i,i) = -sum;
        }

    setExchangeabilityRates();
}

/* The map h_RN from the time-reversible backbone and (u1,u2,u3) to the twelve
   non-reversible weights.  Writing a_k = 2 u_k - 1,

       w_AC = w^R_AC + w^R_CG a1 + w^R_CT a2      w_CA = w^R_AC - w^R_CG a1 - w^R_CT a2
       w_AG = w^R_AG - w^R_CG a1 + w^R_GT a3      w_GA = w^R_AG + w^R_CG a1 - w^R_GT a3
       w_AT = w^R_AT - w^R_CT a2 - w^R_GT a3      w_TA = w^R_AT + w^R_CT a2 + w^R_GT a3
       w_CG = 2 w^R_CG u1                         w_GC = 2 w^R_CG (1 - u1)
       w_CT = 2 w^R_CT u2                         w_TC = 2 w^R_CT (1 - u2)
       w_GT = 2 w^R_GT u3                         w_TG = 2 w^R_GT (1 - u3)

   Each pair sums to twice the corresponding backbone weight, so the total flow
   is preserved at one and every node stays balanced.  Returns false when the
   point lies outside the polyhedron of admissible (u1,u2,u3), that is, when
   some weight would go negative. */
bool RateMatrix_MPQ::buildNonReversibleFromBackbone(const std::vector<mpq_class>& wR, const mpq_class& u1, const mpq_class& u2, const mpq_class& u3) {

    mpq_class a1 = 2 * u1 - 1;
    mpq_class a2 = 2 * u2 - 1;
    mpq_class a3 = 2 * u3 - 1;

    std::vector<mpq_class> w(16);
    w[A*4+C] = wR[W_AC] + wR[W_CG] * a1 + wR[W_CT] * a2;
    w[C*4+A] = wR[W_AC] - wR[W_CG] * a1 - wR[W_CT] * a2;
    w[A*4+G] = wR[W_AG] - wR[W_CG] * a1 + wR[W_GT] * a3;
    w[G*4+A] = wR[W_AG] + wR[W_CG] * a1 - wR[W_GT] * a3;
    w[A*4+T] = wR[W_AT] - wR[W_CT] * a2 - wR[W_GT] * a3;
    w[T*4+A] = wR[W_AT] + wR[W_CT] * a2 + wR[W_GT] * a3;
    w[C*4+G] = 2 * wR[W_CG] * u1;
    w[G*4+C] = 2 * wR[W_CG] * (1 - u1);
    w[C*4+T] = 2 * wR[W_CT] * u2;
    w[T*4+C] = 2 * wR[W_CT] * (1 - u2);
    w[G*4+T] = 2 * wR[W_GT] * u3;
    w[T*4+G] = 2 * wR[W_GT] * (1 - u3);

    return setRatesFromAllWeights(w);
}

/* Recover the point (u1,u2,u3) in the polyhedron that the current non-reversible
   matrix corresponds to.  Since w_CG = 2 w^R_CG u1 and w_GC = 2 w^R_CG (1-u1),
   u1 is just the share of the C <-> G flow that runs from C to G, and likewise
   for u2 and u3.  Each is automatically in (0,1). */
bool RateMatrix_MPQ::recoverU(mpq_class& u1, mpq_class& u2, mpq_class& u3) const {

    mpq_class wCG = pi[C] * (*this)(C,G);
    mpq_class wGC = pi[G] * (*this)(G,C);
    mpq_class wCT = pi[C] * (*this)(C,T);
    mpq_class wTC = pi[T] * (*this)(T,C);
    mpq_class wGT = pi[G] * (*this)(G,T);
    mpq_class wTG = pi[T] * (*this)(T,G);

    mpq_class sCG = wCG + wGC;
    mpq_class sCT = wCT + wTC;
    mpq_class sGT = wGT + wTG;
    if (sCG <= 0 || sCT <= 0 || sGT <= 0)
        return false;

    u1 = wCG / sCG;
    u2 = wCT / sCT;
    u3 = wGT / sGT;
    return true;
}

/* Turn a vector of doubles into exact rationals summing to target.  A double is
   itself a binary rational, so the conversion is exact and the rescaling by a
   single common factor keeps the state exactly proportional to what was drawn.
   The alternative of building the last component as a residual, target minus the
   sum of the others, lets that component drift away from the drawn value and,
   when the concentration is small enough for a draw to underflow, go negative
   outright.  Returns false on a degenerate draw so the caller can reject. */
/* Is a Dirichlet draw safe to build a reverse concentration parameter from?

   A Dirichlet variate is generated as a ratio of gamma variates, and a gamma
   variate with a small shape parameter underflows: at shape 0.005, a quarter of
   a percent of draws come back as denormals and two and a half percent as exact
   zeros.  A denormal passes a plain "> 0" test and then reappears, multiplied by
   the concentration, as the argument to lnGamma, which refuses anything below
   about 1e-300.  Keeping every concentration parameter at or above the proposal
   offset is what actually prevents this, but the check is cheap and it stops a
   future retuning that sets the offset to zero from reintroducing the failure.

   Rejecting here removes proposals of probability on the order of 1e-1000 once
   the offset is in place, which is far below anything the sampler can resolve. */
bool RateMatrix_MPQ::isDrawUsable(const std::vector<double>& x) {

    for (size_t i=0; i<x.size(); i++)
        {
        if (std::isfinite(x[i]) == false || x[i] < MIN_DIRICHLET_DRAW)
            return false;
        }
    return true;
}

bool RateMatrix_MPQ::exactlyNormalize(const std::vector<double>& x, std::vector<mpq_class>& out, const mpq_class& target) {

    size_t n = x.size();
    out.resize(n);
    mpq_class sum = 0;
    for (size_t i=0; i<n; i++)
        {
        if (x[i] <= 0.0 || std::isfinite(x[i]) == false)
            return false;
        out[i] = x[i];
        sum += out[i];
        }
    if (sum <= 0)
        return false;

    mpq_class factor = target / sum;
    for (size_t i=0; i<n; i++)
        out[i] *= factor;
    return true;
}

void RateMatrix_MPQ::calculateWeights(std::vector<mpq_class>& wts) const {

    if (wts.size() != 6)
        throw(RbException("Weights array must have 6 elements"));
        
    // wts[0] = (pi[0] * (*this)(0,1) + pi[1] * (*this)(1,0)) / 2; // (pi[A] * (*this)(A,C) + pi[C] * (*this)(C,A)) / 2
    // wts[1] = (pi[0] * (*this)(0,2) + pi[2] * (*this)(2,0)) / 2; // (pi[A] * (*this)(A,G) + pi[G] * (*this)(G,A)) / 2
    // wts[2] = (pi[0] * (*this)(0,3) + pi[3] * (*this)(3,0)) / 2; // (pi[A] * (*this)(A,T) + pi[T] * (*this)(T,A)) / 2
    // wts[3] = (pi[1] * (*this)(1,2) + pi[2] * (*this)(2,1)) / 2; // (pi[C] * (*this)(C,G) + pi[G] * (*this)(G,C)) / 2
    // wts[4] = (pi[1] * (*this)(1,3) + pi[3] * (*this)(3,1)) / 2; // (pi[C] * (*this)(C,T) + pi[T] * (*this)(T,C)) / 2
    // wts[5] = (pi[2] * (*this)(2,3) + pi[3] * (*this)(3,2)) / 2; // (pi[G] * (*this)(G,T) + pi[T] * (*this)(T,G)) / 2
    for (int i=0, k=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            wts[k++] = (pi[i] * (*this)(i,j) + pi[j] * (*this)(j,i)) / 2;
        }
}

bool RateMatrix_MPQ::check(void) {

    mpq_class averageRate;
    for (int i=0; i<4; i++)
        {
        mpq_class sum;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += (*this)(i,j);
            }
        if ((*this)(i,i) != -sum)
            return false;
        averageRate += pi[i] * sum;
        }
    if (averageRate != 1)
        return false;
    return true;
}

RateMatrix_MPQ* RateMatrix_MPQ::clone( void ) const {

    return new RateMatrix_MPQ( *this );
}

void RateMatrix_MPQ::computeLandU(RateMatrix_MPQ& aMat, RateMatrix_MPQ& lMat, RateMatrix_MPQ& uMat) {

    for (int j=0; j<4; j++)
        {
        for (int k=0; k<j; k++)
            for (int i=k+1; i<j; i++)
                aMat(i,j) = aMat(i,j) - aMat(i,k) * aMat(k,j);

        for (int k=0; k<j; k++)
            for (int i=j; i<4; i++)
                aMat(i,j) = aMat(i,j) - aMat(i,k) * aMat(k,j);

        for (int m=j+1; m<4; m++)
            aMat(m,j) /= aMat(j,j);
        }

    for (int row=0; row<4; row++)
        {
        for (int col=0; col<4; col++)
            {
            if ( row <= col )
                {
                uMat(row,col) = aMat(row,col);
                lMat(row,col) = (row == col ? 1 : 0);
                }
            else
                {
                lMat(row,col) = aMat(row,col);
                uMat(row,col) = 0;
                }
            }
        }
}

double RateMatrix_MPQ::getRate(size_t from, size_t to, double age, double rate) const {

    return (*this)(from,to).get_d() * rate;
}

double RateMatrix_MPQ::getRate(size_t from, size_t to, double rate) const {

    return (*this)(from,to).get_d() * rate;
}

std::vector<double> RateMatrix_MPQ::getRates(void) const {

    std::vector<double> tmp(12);
    size_t k=0;
    for (int i=0;i<4;i++)
        {
        for (int j=0;j<4;j++)
            {
            if ( i != j )
                tmp[k++] = (*this)(i,j).get_d();
            }
        }
    return tmp;
}

std::vector<double> RateMatrix_MPQ::getStationaryFrequencies(void) const {

    std::vector<double> tmp(4);
    for (int i=0;i<4;i++)
        tmp[i] = pi[i].get_d();
    return tmp;
}

/* Draw a time-reversible rate matrix from the prior.

   The prior is Dirichlet(alpha) on the stationary frequencies, independently of a
   uniform draw on the polytope of the model. Those really are independent: pi and
   w are separate coordinates, because reversibility (w_ij = w_ji) and stationarity
   (the circulation condition) are both properties of the weights alone and say
   nothing about pi.

   In weight coordinates the time-reversible polytope is

       { w_ij = w_ji >= 0,  sum of the six backbone weights = 1/2 }

   which is a 5-simplex, so a uniform draw on it is a flat Dirichlet rescaled to
   sum to one half. No rejection is needed. Its volume is (1/2)^5/5! = 1/3840,
   which is the constant computeLnProbability uses, so the two agree by
   construction rather than by coincidence.

   The previous version drew the exchangeability rates from a flat Dirichlet, which
   is NOT uniform on the backbone: r and w are related through a non-linear map
   with a non-constant Jacobian. That only affected where the chain started, so it
   could not bias the posterior, but it did mean this class could not be used with
   the Validation analysis, which starts each replicate from a draw taken by this
   function and checks that the posterior covers the value it started from. */
void RateMatrix_MPQ::initializeTimeReversibleModel(const std::vector<double>& alpha, RandomNumberGenerator* rng) {

    isReversible = true;

    mpq_class one = 1;
    mpq_class oneHalf(1, 2);

    // the stationary frequencies
    std::vector<double> bf = RbStatistics::Dirichlet::rv(alpha, *rng);
    while (isDrawUsable(bf) == false || exactlyNormalize(bf, this->pi, one) == false)
        bf = RbStatistics::Dirichlet::rv(alpha, *rng);

    // the six backbone weights, uniform on the half-simplex
    std::vector<double> flat6(6, 1.0);
    std::vector<mpq_class> wR;
    std::vector<double> bw = RbStatistics::Dirichlet::rv(flat6, *rng);
    while (isDrawUsable(bw) == false || exactlyNormalize(bw, wR, oneHalf) == false)
        bw = RbStatistics::Dirichlet::rv(flat6, *rng);

    setReversibleRatesFromBackbone(wR);

    mpq_class averageRate;
    calculateAverageRate(averageRate);
    if (averageRate != 1)
        throw(RbException("Average rate should be one after initializing the time reversible model"));
}

/* Draw a non-reversible rate matrix from the prior.

   Again pi is independent and drawn from its Dirichlet. The rest has to be uniform
   on the polytope of non-negative unit circulations on K4, which is eight
   dimensional and is not a simplex, so it cannot be hit with a single Dirichlet.

   Carrying it as (backbone, u1, u2, u3) instead, the uniform density picks up the
   Jacobian of that reparameterization, 64 w_CG w_CT w_GT, on the region where the
   point is admissible. The Jacobian can be absorbed exactly rather than handled by
   rejection: multiplying a flat Dirichlet by w_CG w_CT w_GT gives a Dirichlet with
   those three concentrations raised by one, so drawing the backbone from
   Dirichlet(1,1,1,2,2,2) and (u1,u2,u3) uniformly on the cube reproduces the
   target up to the admissibility indicator alone.

   What is left to reject is therefore only the points where some weight would go
   negative, which is about four draws in five. The alternative, rejecting on the
   Jacobian as well, would have thrown away about ninety-nine in a hundred.

   The previous version drew the twelve rates from a flat Dirichlet and derived pi
   from them. That is neither uniform on the polytope nor independent of pi. */
void RateMatrix_MPQ::initializeNonReversibleModel(const std::vector<double>& alpha, RandomNumberGenerator* rng) {

    isReversible = false;

    mpq_class one = 1;
    mpq_class oneHalf(1, 2);

    // the stationary frequencies
    std::vector<double> bf = RbStatistics::Dirichlet::rv(alpha, *rng);
    while (isDrawUsable(bf) == false || exactlyNormalize(bf, this->pi, one) == false)
        bf = RbStatistics::Dirichlet::rv(alpha, *rng);

    /* Concentrations (1,1,1,2,2,2) against the backbone order
       (AC, AG, AT, CG, CT, GT): the three raised entries are exactly the three
       weights that appear in the Jacobian. */
    std::vector<double> tilted6(6, 1.0);
    tilted6[W_CG] = 2.0;
    tilted6[W_CT] = 2.0;
    tilted6[W_GT] = 2.0;

    std::vector<mpq_class> wR;
    for (int attempt = 0; attempt < MAX_INITIALIZATION_ATTEMPTS; attempt++)
        {
        std::vector<double> bw = RbStatistics::Dirichlet::rv(tilted6, *rng);
        if (isDrawUsable(bw) == false || exactlyNormalize(bw, wR, oneHalf) == false)
            continue;

        mpq_class u1 = rng->uniform01();
        mpq_class u2 = rng->uniform01();
        mpq_class u3 = rng->uniform01();

        /* buildNonReversibleFromBackbone checks every weight before it writes any
           of them, so a point outside the polyhedron leaves the matrix untouched
           and the loop can simply try again. */
        if (buildNonReversibleFromBackbone(wR, u1, u2, u3) == true)
            {
            mpq_class averageRate;
            calculateAverageRate(averageRate);
            if (averageRate != 1)
                throw(RbException("Average rate should be one after initializing the non-reversible model"));
            return;
            }
        }

    throw(RbException("Failed to draw a non-reversible rate matrix from the prior."));
}

void RateMatrix_MPQ::moveToDouble( void ) const {
    
    for (int i=0; i<4;i++)
        for (int j=0;j<4;j++)
            (*the_rate_matrix)[i][j] = (*this)(i,j).get_d();
}

void RateMatrix_MPQ::nonreversibilize(mpq_class& u1, mpq_class& u2, mpq_class& u3) {

    if (isReversible == false)
        throw(RbException("Cannot make a non-reversible model non-reversible (again)"));

    // the backbone of a time-reversible matrix is just its six weights
    std::vector<mpq_class> wR(6);
    calculateWeights(wR);

    // we weren't non-reversible before, but we are now (or will be in a microsecond)
    isReversible = false;

    if (buildNonReversibleFromBackbone(wR, u1, u2, u3) == false)
        throw(RbException("The point (u1,u2,u3) lies outside the polyhedron of valid non-reversible models"));

    // the total flow, and hence the average rate, should still be one
    mpq_class averageRate;
    calculateAverageRate(averageRate);
    if (averageRate != 1)
        throw(RbException("Average rate should be one when moving to nonreversible model"));
}

void RateMatrix_MPQ::print(void) {

    std::cout << std::fixed << std::setprecision(8);
    for (int i=0; i<4; i++)
        {
        for (int j=0; j<4; j++)
            {
            if ((*this)(i,j) >= 0)
                std::cout << " ";
            std::cout << (*this)(i,j).get_d() << " ";
            }
        std::cout << std::endl;
        }
}

void RateMatrix_MPQ::reversibilize(void) {

    if (isReversible == true)
        throw(RbException("Cannot reversibilize a time-reversible rate matrix"));
        
    // we weren't reversible before, but we are now (or will be in a microsecond)
    isReversible = true;

    // set off diagonals
    mpq_class w;
    for (int i=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            {
            w = this->pi[i] * (*this)(i,j) + this->pi[j] * (*this)(j,i);
            (*this)(i,j) = w / (2 * this->pi[i]);
            (*this)(j,i) = w / (2 * this->pi[j]);
            }
        }
        
    // set diagonals
    mpq_class sum;
    mpq_class averageRate;
    for (int i=0; i<4; i++)
        {
        sum = 0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += (*this)(i,j);
            }
        (*this)(i,i) = -sum;
        averageRate += pi[i] * sum;
        }
        
    setExchangeabilityRates();
        
    // make certain average rate is one
    if (averageRate != 1)
        throw(RbException("Average rate should be one when moving to a reversible model"));
}

void RateMatrix_MPQ::setExchangeabilityRates(void) {

    if (isReversible == false)
        throw(RbException("Cannot set exchangeability rates for a non-reversible rate matrix"));
        
    // this->r[0] = (*this)(0,1) / pi[1]; // r_AC = Q(A,C) / pi[C]
    // this->r[1] = (*this)(0,2) / pi[2]; // r_AG = Q(A,G) / pi[G]
    // this->r[2] = (*this)(0,3) / pi[3]; // r_AT = Q(A,T) / pi[T]
    // this->r[3] = (*this)(1,2) / pi[2]; // r_CG = Q(C,G) / pi[G]
    // this->r[4] = (*this)(1,3) / pi[3]; // r_CT = Q(C,T) / pi[T]
    // this->r[5] = (*this)(2,3) / pi[3]; // r_GT = Q(G,T) / pi[T]
        
    mpq_class sum;
    for (int i=0, k=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            {
            this->r[k] = (*this)(i,j) / pi[j];
            sum += this->r[k];
            k++;
            }
        }
    for (int i=0; i<6; i++)
        this->r[i] /= sum;
}

void RateMatrix_MPQ::setPi(std::vector<mpq_class>& f) {

    for (int i=0; i<4; i++)
        this->pi[i] = f[i];
}

/** Calculate the transition probabilities for the real case */
void RateMatrix_MPQ::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const {
    
    // get a reference to the eigenvalues
    const std::vector<double>& eigenValue = theEigenSystem->getRealEigenvalues();
    
    // precalculate the product of the eigenvalue and the branch length
    std::vector<double> eigValExp(num_states);
    for (size_t s=0; s<num_states; s++)
        {
        eigValExp[s] = exp(eigenValue[s] * t);
        }
    
    // calculate the transition probabilities
    const double* ptr = &c_ijk[0];
    double*         p = P.theMatrix;
    for (size_t i=0; i<num_states; i++)
        {
        double rowsum = 0.0;
        for (size_t j=0; j<num_states; j++, ++p)
            {
            double sum = 0.0;
            for (size_t s=0; s<num_states; s++)
                {
                sum += (*ptr++) * eigValExp[s];
                }
            
            sum = (sum < 0.0) ? 0.0 : sum;
            rowsum += sum;
            (*p) = sum;
            }

        // Normalize transition probabilities for row to sum to 1.0
        double* p2 = p - num_states;
        for (size_t j=0; j<num_states; j++, ++p2)
            *p2 /= rowsum;
        }
}

/** Calculate the transition probabilities for the complex case */
void RateMatrix_MPQ::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const {
    
    // get a reference to the eigenvalues
    const std::vector<double>& eigenValueReal = theEigenSystem->getRealEigenvalues();
    const std::vector<double>& eigenValueComp = theEigenSystem->getImagEigenvalues();
    
    // precalculate the product of the eigenvalue and the branch length
    std::vector<std::complex<double> > ceigValExp(num_states);
    for (size_t s=0; s<num_states; s++)
        {
        std::complex<double> ev = std::complex<double>(eigenValueReal[s], eigenValueComp[s]);
        ceigValExp[s] = exp(ev * t);
        }
    
    // calculate the transition probabilities
    const std::complex<double>* ptr = &cc_ijk[0];
    for (size_t i=0; i<num_states; i++)
        {
        double rowsum = 0.0;
        for (size_t j=0; j<num_states; j++)
            {
            std::complex<double> sum = std::complex<double>(0.0, 0.0);
            for (size_t s=0; s<num_states; s++)
                sum += (*ptr++) * ceigValExp[s];

            double real_sum = (sum.real() < 0.0) ? 0.0 : sum.real();
            P[i][j] = real_sum;
            rowsum += real_sum;
            }
            
        // normalize transition probabilities for row to sum to 1.0
        for (size_t j=0; j<num_states; j++)
            P[i][j] /= rowsum;
        }
}

void RateMatrix_MPQ::transposeMatrix(const RateMatrix_MPQ& a, RateMatrix_MPQ& t) {
    
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            t(j,i) = a(i,j);
}

/** Update the eigen system */
void RateMatrix_MPQ::updateEigenSystem(void) {
    
    theEigenSystem->update();
    calculateCijk();
}

void RateMatrix_MPQ::update(void) {
    
    updateIfNeeded();
}


/* Bring the double-precision copy of the matrix and its eigensystem into line
   with the rational matrix, if anything has touched the latter since we last did.

   Every write to a rate goes through the non-const operator(), which sets
   needs_update, so "nothing has touched it" is a claim this class can actually
   make rather than one it has to trust its callers for.

   Compile with -DMPQ_VERIFY_EIGEN_CACHE to check that claim on every call: the
   eigensystem is then recomputed even when believed current, and compared against
   what was cached. A stale cache is otherwise silent, giving wrong transition
   probabilities and a wrong likelihood with no error, so the first run of anything
   new is worth doing with this switched on. */
void RateMatrix_MPQ::updateIfNeeded(void) const {

    RateMatrix_MPQ* self = const_cast<RateMatrix_MPQ*>(this);

    if ( needs_update == true )
        {
        moveToDouble();
        self->updateEigenSystem();
        eigen_update_count++;
        needs_update = false;
        return;
        }

#   ifdef MPQ_VERIFY_EIGEN_CACHE
        {
        std::vector<double> cached_c   = c_ijk;
        std::vector<std::complex<double> > cached_cc = cc_ijk;
        bool cached_is_complex = theEigenSystem->isComplex();

        moveToDouble();
        self->updateEigenSystem();
        eigen_update_count++;

        if ( theEigenSystem->isComplex() != cached_is_complex )
            throw(RbException("Stale eigensystem cache: the eigenvalues changed from real to complex or back while the matrix was believed unchanged."));

        for (size_t i=0; i<cached_c.size(); i++)
            {
            if ( fabs(cached_c[i] - c_ijk[i]) > 1e-12 )
                throw(RbException("Stale eigensystem cache: a real c_ijk entry changed while the matrix was believed unchanged. Some code is modifying the rate matrix without going through operator()."));
            }
        for (size_t i=0; i<cached_cc.size(); i++)
            {
            if ( std::abs(cached_cc[i] - cc_ijk[i]) > 1e-12 )
                throw(RbException("Stale eigensystem cache: a complex c_ijk entry changed while the matrix was believed unchanged. Some code is modifying the rate matrix without going through operator()."));
            }
        }
#   endif
}

/* ---------------------------------------------------------------------------
   Moves.

   The state is (pi, w).  For the time-reversible model the free coordinates are
   pi and the six backbone weights, which sum to one half.  For the non-reversible
   model they are pi and the eight free circulation weights, which we carry in the
   equivalent (backbone, u1, u2, u3) parameterization: the backbone is recovered
   as w^R_ij = (w_ij + w_ji)/2 and the u's as the share of each pair's flow that
   runs in the forward direction.

   The prior is flat in w.  In the (backbone, u) parameterization it is therefore
   proportional to the Jacobian of that reparameterization, 64 w^R_CG w^R_CT w^R_GT,
   which is constant in u but not in the backbone.  That is the only non-trivial
   density term any of these moves has to carry.
   --------------------------------------------------------------------------- */

/* Propose new stationary frequencies holding all twelve weights fixed.

   This works for either model.  Reversibility is the statement w_ij = w_ji and
   stationarity is the circulation condition; neither mentions pi, so holding w
   fixed keeps the matrix in whichever model it was already in, and keeps the
   average rate sum_ij w_ij at one.  The state coordinate w does not move, so
   there is no Jacobian and the Hastings ratio is just the ratio of the two
   Dirichlet proposal densities.  (The old version of this function held the
   exchangeability rates fixed instead, which drags w along a non-linear path
   and needs a Jacobian that was not there.) */
double RateMatrix_MPQ::updateStationaryFrequencies(RandomNumberGenerator* rng, double alpha0, double offset) {

    mpq_class one = 1;

    std::vector<mpq_class> w;
    calculateAllWeights(w);

    std::vector<double> oldFreqs(4);
    std::vector<double> alphaForward(4);
    std::vector<double> alphaReverse(4);
    for (int i=0; i<4; i++)
        oldFreqs[i] = pi[i].get_d();
    for (int i=0; i<4; i++)
        alphaForward[i] = oldFreqs[i] * alpha0 + offset;

    std::vector<double> newFreqs = RbStatistics::Dirichlet::rv(alphaForward, *rng);
    if (isDrawUsable(newFreqs) == false)
        return RbConstants::Double::neginf;
    std::vector<mpq_class> newPi;
    if (exactlyNormalize(newFreqs, newPi, one) == false)
        return RbConstants::Double::neginf;
    for (int i=0; i<4; i++)
        alphaReverse[i] = newFreqs[i] * alpha0 + offset;

    pi = newPi;
    if (setRatesFromAllWeights(w) == false)
        return RbConstants::Double::neginf;
    if (isReversible == true)
        setExchangeabilityRates();

    return RbStatistics::Dirichlet::lnPdf(alphaReverse, oldFreqs) -
           RbStatistics::Dirichlet::lnPdf(alphaForward, newFreqs);
}

/* The single-element version of the move above: pick one nucleotide, redraw its
   frequency, and scale the other three by (1 - pi'_k) / (1 - pi_k) so that they
   keep their relative proportions and the four still sum to one.  The Jacobian
   of that scaling on the remaining two free coordinates is factor^(4-2).

   NB the previous version scaled the others by pi_k / pi'_k, which does not put
   pi'_k at the value that was drawn, so the Beta density was being evaluated at
   a point the chain never visited. */
double RateMatrix_MPQ::updateStationaryFrequenciesSingle(RandomNumberGenerator* rng, double alpha0, double offset) {

    mpq_class one = 1;

    std::vector<mpq_class> w;
    calculateAllWeights(w);

    size_t index = size_t(rng->uniform01() * 4);
    if (index > 3)
        index = 3;

    std::vector<double> oldVals(2);
    std::vector<double> alphaForward(2);
    std::vector<double> alphaReverse(2);
    oldVals[0] = pi[index].get_d();
    oldVals[1] = 1.0 - oldVals[0];
    alphaForward[0] = oldVals[0] * alpha0 + offset;
    alphaForward[1] = oldVals[1] * alpha0 + offset;

    std::vector<double> newVals = RbStatistics::Dirichlet::rv(alphaForward, *rng);
    if (isDrawUsable(newVals) == false)
        return RbConstants::Double::neginf;
    alphaReverse[0] = newVals[0] * alpha0 + offset;
    alphaReverse[1] = newVals[1] * alpha0 + offset;

    double factor = newVals[1] / oldVals[1];
    std::vector<double> proposed(4);
    for (int i=0; i<4; i++)
        {
        if (i == (int)index)
            proposed[i] = newVals[0];
        else
            proposed[i] = pi[i].get_d() * factor;
        }

    std::vector<mpq_class> newPi;
    if (exactlyNormalize(proposed, newPi, one) == false)
        return RbConstants::Double::neginf;

    pi = newPi;
    if (setRatesFromAllWeights(w) == false)
        return RbConstants::Double::neginf;
    if (isReversible == true)
        setExchangeabilityRates();

    double lnProposalProb = RbStatistics::Dirichlet::lnPdf(alphaReverse, oldVals) -
                            RbStatistics::Dirichlet::lnPdf(alphaForward, newVals);
    lnProposalProb += (4 - 2) * log(factor);
    return lnProposalProb;
}

/* Propose all six backbone weights of a time-reversible matrix, holding pi fixed.

   The backbone is the state coordinate and the prior is flat on it, so a Dirichlet
   proposal on the simplex needs nothing else: the fixed rescaling between the unit
   simplex the Dirichlet lives on and the one-half simplex the weights live on
   contributes the same constant in both directions and cancels.  (The old version
   proposed exchangeability rates, for which the target is not flat.) */
double RateMatrix_MPQ::updateBackboneWeights(RandomNumberGenerator* rng, double alpha0, double offset) {

    if (isReversible == false)
        throw(RbException("Can only update the backbone weights directly for time reversible models"));

    mpq_class oneHalf(1, 2);

    std::vector<mpq_class> wR(6);
    calculateWeights(wR);

    std::vector<double> oldW(6);
    std::vector<double> alphaForward(6);
    std::vector<double> alphaReverse(6);
    for (int i=0; i<6; i++)
        oldW[i] = 2.0 * wR[i].get_d();
    for (int i=0; i<6; i++)
        alphaForward[i] = oldW[i] * alpha0 + offset;

    std::vector<double> newW = RbStatistics::Dirichlet::rv(alphaForward, *rng);
    if (isDrawUsable(newW) == false)
        return RbConstants::Double::neginf;
    std::vector<mpq_class> newWR;
    if (exactlyNormalize(newW, newWR, oneHalf) == false)
        return RbConstants::Double::neginf;
    for (int i=0; i<6; i++)
        alphaReverse[i] = newW[i] * alpha0 + offset;

    setReversibleRatesFromBackbone(newWR);

    return RbStatistics::Dirichlet::lnPdf(alphaReverse, oldW) -
           RbStatistics::Dirichlet::lnPdf(alphaForward, newW);
}

/* The single-element version: redraw one backbone weight and scale the other five
   to keep the sum at one half.  Jacobian factor^(6-2). */
double RateMatrix_MPQ::updateBackboneWeightsSingle(RandomNumberGenerator* rng, double alpha0, double offset) {

    if (isReversible == false)
        throw(RbException("Can only update the backbone weights directly for time reversible models"));

    mpq_class oneHalf(1, 2);

    std::vector<mpq_class> wR(6);
    calculateWeights(wR);

    size_t index = size_t(rng->uniform01() * 6);
    if (index > 5)
        index = 5;

    std::vector<double> oldVals(2);
    std::vector<double> alphaForward(2);
    std::vector<double> alphaReverse(2);
    oldVals[0] = 2.0 * wR[index].get_d();
    oldVals[1] = 1.0 - oldVals[0];
    alphaForward[0] = oldVals[0] * alpha0 + offset;
    alphaForward[1] = oldVals[1] * alpha0 + offset;

    std::vector<double> newVals = RbStatistics::Dirichlet::rv(alphaForward, *rng);
    if (isDrawUsable(newVals) == false)
        return RbConstants::Double::neginf;
    alphaReverse[0] = newVals[0] * alpha0 + offset;
    alphaReverse[1] = newVals[1] * alpha0 + offset;

    double factor = newVals[1] / oldVals[1];
    std::vector<double> proposed(6);
    for (int i=0; i<6; i++)
        {
        if (i == (int)index)
            proposed[i] = newVals[0];
        else
            proposed[i] = 2.0 * wR[i].get_d() * factor;
        }

    std::vector<mpq_class> newWR;
    if (exactlyNormalize(proposed, newWR, oneHalf) == false)
        return RbConstants::Double::neginf;

    setReversibleRatesFromBackbone(newWR);

    double lnProposalProb = RbStatistics::Dirichlet::lnPdf(alphaReverse, oldVals) -
                            RbStatistics::Dirichlet::lnPdf(alphaForward, newVals);
    lnProposalProb += (6 - 2) * log(factor);
    return lnProposalProb;
}

/* Propose a new backbone for a non-reversible matrix, holding (u1,u2,u3) and pi
   fixed.

   The state is the eight free non-reversible weights, but the move is expressed
   in (backbone, u) coordinates, in which the target is not flat: the Jacobian of
   that reparameterization is 64 w^R_CG w^R_CT w^R_GT.  Proposing a new backbone
   at fixed u therefore needs the ratio of that quantity at the new and old
   backbones.  Without it the chain under-weights CG, CT and GT relative to AC,
   AG and AT by several percent.

   The proposed backbone may put the retained u outside the polyhedron, in which
   case some weight goes negative, the point is outside the support of the prior,
   and the move is rejected. */
double RateMatrix_MPQ::updateNonReversibleBackbone(RandomNumberGenerator* rng, double alpha0, double offset) {

    if (isReversible == true)
        throw(RbException("Can only update the non-reversible backbone for non-reversible models"));

    mpq_class oneHalf(1, 2);

    std::vector<mpq_class> wR(6);
    calculateWeights(wR);

    mpq_class u1, u2, u3;
    if (recoverU(u1, u2, u3) == false)
        return RbConstants::Double::neginf;

    std::vector<double> oldW(6);
    std::vector<double> alphaForward(6);
    std::vector<double> alphaReverse(6);
    for (int i=0; i<6; i++)
        oldW[i] = 2.0 * wR[i].get_d();
    for (int i=0; i<6; i++)
        alphaForward[i] = oldW[i] * alpha0 + offset;

    std::vector<double> newW = RbStatistics::Dirichlet::rv(alphaForward, *rng);
    if (isDrawUsable(newW) == false)
        return RbConstants::Double::neginf;
    std::vector<mpq_class> newWR;
    if (exactlyNormalize(newW, newWR, oneHalf) == false)
        return RbConstants::Double::neginf;
    for (int i=0; i<6; i++)
        alphaReverse[i] = newW[i] * alpha0 + offset;

    if (buildNonReversibleFromBackbone(newWR, u1, u2, u3) == false)
        return RbConstants::Double::neginf;

    double lnProposalProb = RbStatistics::Dirichlet::lnPdf(alphaReverse, oldW) -
                            RbStatistics::Dirichlet::lnPdf(alphaForward, newW);

    // the reparameterization Jacobian, 64 w_CG w_CT w_GT, at the new backbone
    // relative to the old one
    lnProposalProb += log(newWR[W_CG].get_d()) + log(newWR[W_CT].get_d()) + log(newWR[W_GT].get_d());
    lnProposalProb -= log(wR[W_CG].get_d())    + log(wR[W_CT].get_d())    + log(wR[W_GT].get_d());

    return lnProposalProb;
}

/* Propose a new point (u1,u2,u3) in the polyhedron, holding the backbone and pi
   fixed.  In (backbone, u) coordinates the target is proportional to
   64 w^R_CG w^R_CT w^R_GT, which does not involve u at all, so a symmetric
   random walk on u has a Hastings ratio of one.  A step that leaves the unit
   cube, or that leaves the polyhedron and so drives a weight negative, falls
   outside the support and is rejected. */
double RateMatrix_MPQ::updateNonReversibleU(RandomNumberGenerator* rng, double delta) {

    if (isReversible == true)
        throw(RbException("Can only update (u1,u2,u3) for non-reversible models"));

    std::vector<mpq_class> wR(6);
    calculateWeights(wR);

    mpq_class u[3];
    if (recoverU(u[0], u[1], u[2]) == false)
        return RbConstants::Double::neginf;

    for (int i=0; i<3; i++)
        {
        double proposed = u[i].get_d() + delta * (rng->uniform01() - 0.5);
        if (proposed <= 0.0 || proposed >= 1.0)
            return RbConstants::Double::neginf;
        u[i] = proposed;
        }

    if (buildNonReversibleFromBackbone(wR, u[0], u[1], u[2]) == false)
        return RbConstants::Double::neginf;

    return 0.0;
}

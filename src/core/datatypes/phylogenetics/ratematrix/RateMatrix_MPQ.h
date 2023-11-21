#ifndef RateMatrix_MPQ_h
#define RateMatrix_MPQ_h

#include <gmpxx.h>
#include <iomanip>
#include <iostream>
#include <vector>


#include "RateMatrix.h"


namespace RevBayesCore {


class MatrixReal;
class RandomNumberGenerator;
class TransitionProbabilityMatrix;

class RateMatrix_MPQ : public RateMatrix {
    
    /**
     * This class is used to hold a 4 X 4 rate matrix for a continuous-time Markov model. The
     * implementation of this class relies on GMP rationals to hold rates and to calculate
     * stationary frequencies or exchangeability parameters. This class also has functionality
     * to propose a reversible model from a non-reversible one and also to do the reverse.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (John Huelsenbeck)
     * @since 2014-11-18, version 1.0
     */
    
public:
    RateMatrix_MPQ(void);
    RateMatrix_MPQ(const RateMatrix_MPQ& m);
    ~RateMatrix_MPQ(void);
    mpq_class&              operator()(size_t r, size_t c) { return this->q[r * 4 + c]; }
    const mpq_class&        operator()(size_t r, size_t c) const { return this->q[r * 4 + c]; }
    RateMatrix_MPQ&         operator=(const RateMatrix_MPQ& rhs);
    void                    adjust(void);
    void                    calculateAverageRate(mpq_class& ave) const;
    void                    calculateStationaryFrequencies(std::vector<mpq_class>& f);
    void                    calculateWeights(std::vector<mpq_class>& wts);
    bool                    check(void);
    std::vector<mpq_class>& getExchangeabilityRates(void) { return r; }
    bool                    getIsReversible(void) { return isReversible; }
    std::vector<mpq_class>& getPi(void) { return pi; }
    std::vector<double>     getRates(void) const ;
    void                    initializeTimeReversibleModel(const std::vector<double>& alpha, RandomNumberGenerator* rng);
    void                    initializeNonReversibleModel(const std::vector<double>& alpha, RandomNumberGenerator* rng);
    void                    nonreversibilize(mpq_class& u1, mpq_class& u2, mpq_class& u3);
    void                    print(void);
    void                    reversibilize(void);
    void                    setExchangeabilityRates(void);
    void                    setIsReversible(bool tf) { isReversible = tf; }
    void                    setPi(std::vector<mpq_class>& f);
    double                  updateNonReversibleRates(RandomNumberGenerator* rng, double alpha0);
    double                  updateNonReversibleRatesSingle(RandomNumberGenerator* rng, double alpha0);
    double                  updateExchangeabilityRates(RandomNumberGenerator* rng, double alpha0);
    double                  updateExchangeabilityRatesSingle(RandomNumberGenerator* rng, double alpha0);
    double                  updateStationaryFrequencies(RandomNumberGenerator* rng, double alpha0);
    double                  updateStationaryFrequenciesSingle(RandomNumberGenerator* rng, double alpha0);

    // virtual methods from RateMatrix
    double                  averageRate(void) const;                                                                //!< Calculate the average rate
    void                    calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrixrate matrix
    RateMatrix_MPQ*         clone(void) const;
    double                  getRate(size_t from, size_t to, double age, double rate) const;                         //!< Calculate the rate from state i to state j over the given time interval scaled by a rate
    double                  getRate(size_t from, size_t to, double rate=1.0) const;
    std::vector<double>     getStationaryFrequencies(void) const;                                                   //!< Return the stationary frequencies
    void                    rescaleToAverageRate(double r) { throw RbException("We do not support rescaling of the non-reversible rate matrix"); }                                                         //!< Rescale the rate matrix such that the average rate is "r"
    void                    setDiagonal(void) { throw RbException("We do not support to set the diagonal of a non-reversible rate matrix."); }                                                                      //!< Set the diagonal such that each row sums to zero
    void                    update(void);                                                                           //!< Update the rate entries of the matrix (is needed if stationarity freqs or similar have changed)

    
private:
    void                    computeLandU(RateMatrix_MPQ& aMat, RateMatrix_MPQ& lMat, RateMatrix_MPQ& uMat);
    void                    transposeMatrix(const RateMatrix_MPQ& a, RateMatrix_MPQ& t);
    mpq_class*              q;                     // elements of the rate matrix
    mpq_class*              endBuffer;             // memory one past the end of the rate matrix array
    bool                    isReversible;          // flag indicating whether or not this rate matrix is time reversible
    std::vector<mpq_class>  pi;                    // the stationary frequencies of the rate matrix
    std::vector<mpq_class>  r;                     // the exchangeability parameters, if time reversible
    
    void                    moveToDouble(void) const;
    
    void                                calculateCijk(void);                                                                //!< Do precalculations on eigenvectors and their inverse
    void                                tiProbsEigens(double t, TransitionProbabilityMatrix& P) const;                      //!< Calculate transition probabilities for real case
    void                                tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const;               //!< Calculate transition probabilities for complex case
    void                                updateEigenSystem(void);                                                            //!< Update the system of eigenvalues and eigenvectors
    
    MatrixReal*                         the_rate_matrix;                                                                    //!< Holds the rate matrix
    EigenSystem*                        theEigenSystem;                                                                     //!< Holds the eigen system
    std::vector<double>                 c_ijk;                                                                              //!< Vector of precalculated product of eigenvectors and their inverse
    std::vector<std::complex<double> >  cc_ijk;                                                                             //!< Vector of precalculated product of eigenvectors and thier inverse for complex case
    bool                                needs_update;

    
    friend std::ostream& operator<<(std::ostream& os, RateMatrix_MPQ& m);
};

}

inline std::ostream& operator<<(std::ostream& os, RevBayesCore::RateMatrix_MPQ& m) {

    std::vector<mpq_class> bf(4);
    m.calculateStationaryFrequencies(bf);
    os << std::fixed << std::setprecision(10);
    mpq_class averageRate;
    for (int i=0; i<4; i++)
        {
        mpq_class sum;
        for (int j=0; j<4; j++)
            {
            if (m(i,j) > 0)
                os << " ";
            os << m(i,j) << " ";
            sum += m(i,j);
            if (i != j)
                averageRate += bf[i] * m(i,j);
            }
        os << "(" << sum << ")";
        os << std::endl;
        }
    os << "Average Rate = " << averageRate << std::endl;
    return os;
}

#endif

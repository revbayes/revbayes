/**
 * @brief
 * This file contains the declaration of the Multivariate Normally distributed random variable class.
 * This class is derived from the stochastic node and each instance will represent a random variable
 * from a normal distribution in the model graph.
 * The Multivariate Normal distribution represents a family of distributions
 * defined on the real numbers. The Multivariate Normal distribution has a varying number of parameters depending on the dimensionality of the distribution:
 *
 *@param Mu A location parameter for each dimension
 *@param Sigma a variance-covariance matrix
 */



#ifndef BirthDeathJeffreysDistribution_H
#define BirthDeathJeffreysDistribution_H

#include <cstddef>
#include <vector>
#include <iosfwd>

#include "RbVector.h"
#include "TypedDistribution.h"
#include "MatrixReal.h"
#include "RevPtr.h"

namespace RevLanguage { class RevVariable; }

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;
    
    class BirthDeathJeffreysDistribution : public TypedDistribution< RbVector<double> > {
        
    public:
        BirthDeathJeffreysDistribution(const TypedDagNode<double> *ra,
                					   const std::string &cdt,
                                       bool uo,
									   const TypedDagNode<double> *r,
									   double l);
        virtual                                                         ~BirthDeathJeffreysDistribution(void);                                                //!< Virtual destructor
        
        
        // public member functions
        BirthDeathJeffreysDistribution*                                 clone(void) const;
        double                                                          computeLnProbability(void);
        void                                                            redrawValue(void);

    protected:

        // Parameter management functions
        void                                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Swap a parameter
        
    private:
        
        // members
        const TypedDagNode<double>*       process_age;                                                                                        //!< Time since the start of the process.
        bool                              use_origin;
        std::string                       condition;                                          //!< The condition of the process (none/survival/#taxa).
        const TypedDagNode<double>*       rho;                                                                                        //!< Time since the start of the process.
        double                            limit;

    };
    
}

#endif

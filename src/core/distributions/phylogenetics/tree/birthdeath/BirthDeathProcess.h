#ifndef BirthDeathProcess_H
#define BirthDeathProcess_H

#include "Taxon.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "AbstractBirthDeathProcess.h"

namespace RevBayesCore {

    class Clade;

    /**
     * @file
     * This file contains the declaration of the random variable class for constant rate birth-death process.
     *
     * @brief Declaration of the abstract Birth-Death process class.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-17, version 1.0
     *
     */
    class BirthDeathProcess : public AbstractBirthDeathProcess {

    public:
        BirthDeathProcess(const TypedDagNode<double> *ro,
                          const TypedDagNode<double> *rh,
                          const TypedDagNode<double> *mp,
                          const std::string& ss,
                          const std::vector<Clade> &ic,
                          const std::string &cdt,
                          const std::vector<Taxon> &tn,
                          Tree *t);

        // pure virtual member functions
        virtual BirthDeathProcess*                          clone(void) const = 0;                                                      //!< Create an independent clone


    protected:
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        virtual void                                        restoreSpecialization(const DagNode *restorer);
        virtual void                                        touchSpecialization(const DagNode *toucher, bool touchAll);

        // pure virtual helper functions
        virtual double                                      lnSpeciationRate(double t) const = 0;                                       //!< Get the log-transformed speciation rate at time t.
        virtual double                                      rateIntegral(double t_low, double t_high) const = 0;                        //!< Compute the rate integral.
        virtual double                                      computeProbabilitySurvival(double start, double end) const = 0;             //!< Compute the probability of survival of the process (without incomplete taxon sampling).

        virtual void                                        prepareRateIntegral(double end) const;                                      //!< Compute the rate integral.
        virtual void                                        prepareSurvivalProbability(double end, double r) const;                     //!< Compute the rate integral.
        
        
        // helper functions
        virtual double                                      computeLnProbabilityTimes(void) const;                                      //!< Compute the log-transformed probability of the current value.
        size_t                                              getNumberOfTaxaAtPresent(void) const;
        double                                              lnP1(double T, double r) const;
        double                                              lnP1(double t, double T, double r) const;
        double                                              lnProbNumTaxa(size_t n, double start, double end, bool MRCA) const;         //!< Compute the log-transformed probability of the number of taxa.
        double                                              lnProbSurvival(double start, double end) const;                             //!< Compute the probability of survival of the process (without incomplete taxon sampling).
        double                                              lnProbSurvival(double start, double end, double r) const;                   //!< Compute the probability of survival of the process including uniform taxon sampling.
        double                                              pSurvival(double start, double end) const;                                  //!< Compute the probability of survival of the process (without incomplete taxon sampling).
        double                                              pSurvival(double start, double end, double r) const;                        //!< Compute the probability of survival of the process including uniform taxon sampling.

        // members
        const TypedDagNode<double>*                         rho;                                                                        //!< Sampling probability of each species.
        const TypedDagNode<double>*                         sampling_mixture_proportion;                                                //!< The proportion of sampling that is diversified in sampling mixture model.
        std::string                                         sampling_strategy;                                                          //!< The incomplete taxon sampling strategy (uniform/diversified).
        std::vector<int>                                    missing_species;
        std::vector<Clade>                                  incomplete_clades;                                                          //!< Topological constrains.
        std::vector<double>                                 incomplete_clade_ages;                                                      //!< Topological constrains.
        
        mutable std::vector<double>                         log_p_survival;                                                             //!< Topological constrains.
        mutable std::vector<double>                         rate_integral;                                                              //!< Topological constrains.
        
        
    };

}

#endif

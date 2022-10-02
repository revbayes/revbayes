#ifndef UNIT_MIXTURE_MODEL_H
#define UNIT_MIXTURE_MODEL_H

#include <vector>

#include "MixtureModel.h"
#include "RateMatrix.h"

namespace RevBayesCore {

    class Tree;

    struct UnitMixtureModel: public MixtureModel
    {
        std::unique_ptr<RateMatrix> matrix;

        UnitMixtureModel* clone() const;
        void calculateTransitionProbabilities(const Tree& tau, int node_index, int m, TransitionProbabilityMatrix& P) const;

        std::vector<double> getRootFrequencies(int) const;

        std::vector<double> componentProbs() const;

        double rate() const;
        void scale(double factor);
        void setRate(double r);
        
        UnitMixtureModel& operator=(const UnitMixtureModel&);
        
        UnitMixtureModel(const RateMatrix& m);
        UnitMixtureModel(const UnitMixtureModel& m);
    };

}

#endif

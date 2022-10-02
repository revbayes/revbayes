#ifndef RATE_MIXTURE_MODEL_H
#define RATE_MIXTURE_MODEL_H

#include <vector>

#include "MixtureModel.h"

namespace RevBayesCore {

    class Tree;

    struct ConcreteMixtureModel: public MixtureModel
    {
        std::vector<std::unique_ptr<MixtureModel>> sub_models;
        std::vector<double> fractions;

        ConcreteMixtureModel* clone() const;
        void calculateTransitionProbabilities(const Tree& tau, int node_index, int m, TransitionProbabilityMatrix& P) const;

        std::vector<double> getRootFrequencies(int) const;

        std::vector<double> componentProbs() const;

        void scale(double factor);

        void setRate(double r);

        double rate() const;
        
        ConcreteMixtureModel& operator=(const ConcreteMixtureModel&);
        ConcreteMixtureModel(const ConcreteMixtureModel&);

        ConcreteMixtureModel(const std::vector<std::unique_ptr<MixtureModel>>& m, const std::vector<double>& fs);
    };


    ConcreteMixtureModel* scaled_mixture(const MixtureModel& sub_model, const std::vector<double>& fs, const std::vector<double>& rs);
}

#endif

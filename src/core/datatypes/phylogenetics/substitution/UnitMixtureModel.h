#ifndef UNIT_MIXTURE_MODEL_H
#define UNIT_MIXTURE_MODEL_H

#include "MixtureModel.h"
#include "RateGenerator.h"
#include "RateMatrix.h"

namespace RevBayesCore {

    class Tree;

    class UnitMixtureModel: public MixtureModel
    {
        std::unique_ptr<RateGenerator> generator;
        std::vector<double> frequencies;
        double _scale;

    public:
        UnitMixtureModel* clone() const;
        TransitionProbabilityMatrix calculateTransitionProbabilities(const Tree& tau, int node_index, int m, double rate) const;

        bool simulateStochasticMapping(const Tree&, int node, int mixture_component, int rate, std::vector<size_t>&, std::vector<double>&);

        std::vector<double> getRootFrequencies(int) const;

        std::vector<double> componentProbs() const;

        std::optional<double> rate() const;
        void scale(double factor);
        
        UnitMixtureModel& operator=(const UnitMixtureModel&);
        UnitMixtureModel(const UnitMixtureModel&);

        UnitMixtureModel& operator=(UnitMixtureModel&&);
        UnitMixtureModel(UnitMixtureModel&&);
        
        UnitMixtureModel(const RateMatrix& m, double s = 1);
        UnitMixtureModel(const RateMatrix& m, const std::vector<double>& freqs, double s = 1);
        UnitMixtureModel(const RateGenerator& g, const std::vector<double>& freqs, double s = 1);
    };

}

#endif

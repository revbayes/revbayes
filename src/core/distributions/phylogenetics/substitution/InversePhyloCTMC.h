#ifndef InversePhyloCTMC_H
#define InversePhyloCTMC_H

#include "Dist_PhyloCTMC.h"

class InversePhyloCTMC {
public:
    static std::shared_ptr<Dist_phyloCTMC> create(const std::shared_ptr<Dist_phyloCTMC>& base_dist) {
        // Create a new phyloCTMC instance
        auto new_dist = std::make_shared<>(*base_dist); // Assuming copy constructor

        // Override calcLnProbability
        new_dist->calcLnProbability = [base_dist](void) {
            return -base_dist->calcLnProbability();
        };

        return new_dist;
    }
};

#endif

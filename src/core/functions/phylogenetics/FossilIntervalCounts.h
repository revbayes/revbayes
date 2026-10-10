/**
 * @file
 * This file contains the per-interval taxon counts that the classical paleontological rate
 * estimators are built from, and the presence matrix they are derived from.
 *
 * @brief Per-interval taxon counts for the classical fossil rate estimators.
 *
 * @author June Walker
 * @since 2026-10-08, version 1.3
 *
 */

#ifndef FossilIntervalCounts_H
#define FossilIntervalCounts_H

#include <cstddef>
#include <vector>

#include "Taxon.h"

namespace RevBayesCore {

    //!< Where an occurrence bin that straddles an interval boundary is counted.
    enum FossilAgeAmbiguity { FOSSIL_AGE_MIDPOINT, FOSSIL_AGE_OVERLAP, FOSSIL_AGE_EXCLUDE };

    /**
     * @brief Counts of the taxon categories in each interval, Foote's four and then Alroy's.
     */
    struct FossilIntervalCounts {
        std::vector<double> singleton;                                                  //!< N_FL, first and last appearance both inside
        std::vector<double> bottom_last;                                                //!< N_bL, crosses the older edge, last appearance inside
        std::vector<double> first_top;                                                  //!< N_Ft, first appearance inside, crosses the younger edge
        std::vector<double> through;                                                    //!< N_bt, crosses both edges
        std::vector<double> two_timer_older;                                            //!< N_2t,i, sampled here and in the older neighbour
        std::vector<double> two_timer_younger;                                          //!< N_2t,i+1, sampled here and in the younger neighbour
        std::vector<double> three_timer;                                                //!< N_3t, sampled here and in both neighbours
        std::vector<double> part_timer;                                                 //!< N_pt, sampled in both neighbours but not here
    };

    /*
     * Boundaries are the interval start ages, youngest first, as the fossilized birth-death
     * processes take them: the present, then each rate shift time. Interval j is
     * [boundaries[j], boundaries[j+1]), half-open at its older edge, and the oldest interval
     * runs to top, which is infinite unless the caller bounds it.
     */
    std::vector<std::vector<bool> > fossilIntervalPresence(const std::vector<Taxon> &taxa,
                                                           const std::vector<double> &boundaries,
                                                           double top,
                                                           FossilAgeAmbiguity policy);  //!< Which intervals each taxon was sampled in.
    FossilIntervalCounts            countFossilIntervals(const std::vector<std::vector<bool> > &presence);           //!< The category counts, from that presence matrix.
    size_t                          fossilIntervalOf(double age, const std::vector<double> &boundaries, double top); //!< The interval holding an age, or the interval count if it falls outside.
    void                            checkFossilIntervals(const std::vector<double> &boundaries);                     //!< Throws unless the boundaries are finite and ascending.
    bool                            fossilIntervalsEqualWidth(const std::vector<double> &boundaries, double top);    //!< The three-timer method assumes they are.

}

#endif

/**
 * @file
 * This file contains the declaration of Signals, which is
 * used to hold global signals.
 *
 * @brief Declaration of Signals
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @since 2010-09-19
 *
 * $Id$
 */

#ifndef Signals_H
#define Signals_H

#include <cstdint>

class Signals {

    public:
        enum                    RbFlags {   BREAK       = 1L << 0,
                                            CONTINUE    = 1L << 1,
                                            RETURN      = 1L << 2 };                            //!< RevBayes flags

        void                    set(const std::uint64_t flag) { flags |= flag; }                         //!< Set a flag
        void                    clearFlags(void) { flags = 0L; }                                //!< Clear flags
        bool                    isGood(void) const { return flags == 0; }                       //!< Test all flags
        bool                    isSet(const std::uint64_t flag) const { return (flags & flag) != 0; }    //!< Test a flag

        /** Get unique signals instance */
        static Signals& getSignals() {
                static Signals theSignals = Signals();
                return theSignals;
        }

    private:
                                Signals(void) : flags(0L) {}                                    //!< Prevent construction
                                Signals(const Signals& w) {}                                    //!< Prevent copy

        Signals&                operator=(const Signals& x);                                    //!< Prevent assignment

        std::uint64_t                    flags;                                                          //!< The flags
};

#endif


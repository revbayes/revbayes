#include "ScreenMonitor.h"

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>

#include "DagNode.h"
#include "Model.h"
#include "Cloneable.h"
#include "StringUtilities.h"

using namespace RevBayesCore;

/* Constructor */
ScreenMonitor::ScreenMonitor(DagNode *n, std::uint64_t g, bool pp, bool l, bool pr) : Monitor(g,n),
    posterior( pp ),
    prior( pr ),
    likelihood( l ),
    printWaitingTime( true ),
    printElapsedTime( true ),
    prefixSeparator("   "),
    suffixSeparator("   |"),
    headerPrintingInterval( 20 ),
    startTime( 0 ),
    numCycles( 0 ),
    currentGen( 0 ),
    startGen( 0 ),
    replicateIndex( 0 )
{}


ScreenMonitor::ScreenMonitor(const std::vector<DagNode *> &n, std::uint64_t g, bool pp, bool l, bool pr) : Monitor(g,n),
    posterior( pp ),
    prior( pr ),
    likelihood( l ),
    printWaitingTime( true ),
    printElapsedTime( true ),
    prefixSeparator("   "),
    suffixSeparator("   |"),
    headerPrintingInterval( 20 ),
    startTime( 0 ),
    numCycles( 0 ),
    currentGen( 0 ),
    startGen( 0 ),
    replicateIndex( 0 )
{}


/* Clone the object */
ScreenMonitor* ScreenMonitor::clone(void) const
{
    
    return new ScreenMonitor(*this);
}

/** Is it a screen monitor ? Yes. **/
bool ScreenMonitor::isScreenMonitor( void ) const
{
    return true;
}


void ScreenMonitor::monitor(std::uint64_t gen)
{

    // only the enabled monitor is allowed to print
    if ( enabled == true )
    {
        
        // get the printing frequency
        std::uint64_t samplingFrequency = printgen;
        
        if ( gen > 0 && gen % (headerPrintingInterval*samplingFrequency) == 0 )
        {
            std::cout << std::endl;
            printHeader();
        }
        
        
        if (gen % samplingFrequency == 0)
        {
            // handy string and stringstream
            std::string s;
            std::stringstream ss;
            
            // set column width
            int columnWidth = 12;
            
            // set cycle column width
            int cycleWidth = floor( log10( numCycles ) ) + 1;
            cycleWidth = 9 > cycleWidth ? 9 : cycleWidth;

            // print the cycle number
            ss << gen;
            s = ss.str();
            StringUtilities::fillWithSpaces(s, cycleWidth, true);
            std::cout << s << suffixSeparator;
            ss.str("");

            if ( posterior )
            {
                const std::vector<DagNode*> &n = model->getDagNodes();
                double pp = 0.0;
                for (std::vector<DagNode*>::const_iterator it = n.begin(); it != n.end(); ++it)
                {
                    pp += (*it)->getLnProbability();
                }
                
                ss << pp;
                s = ss.str();
                StringUtilities::fillWithSpaces( s, columnWidth, false );
                std::cout << prefixSeparator << s << suffixSeparator;
                ss.str("");
            }
            
            if ( likelihood )
            {
                const std::vector<DagNode*> &n = model->getDagNodes();
                double pp = 0.0;
                for (std::vector<DagNode*>::const_iterator it = n.begin(); it != n.end(); ++it)
                {
                    if ( (*it)->isClamped() )
                    {
                        pp += (*it)->getLnProbability();
                    }
                }

                ss << pp;
                s = ss.str();
                StringUtilities::fillWithSpaces( s, columnWidth, false );
                std::cout << prefixSeparator << s << suffixSeparator;
                ss.str("");
            }
            
            if ( prior )
            {
                const std::vector<DagNode*> &n = model->getDagNodes();
                double pp = 0.0;
                for (std::vector<DagNode*>::const_iterator it = n.begin(); it != n.end(); ++it)
                {
                    if ( !(*it)->isClamped() )
                    {
                        pp += (*it)->getLnProbability();
                    }
                }

                ss << pp;
                s = ss.str();
                StringUtilities::fillWithSpaces( s, columnWidth, false );
                std::cout << prefixSeparator << s << suffixSeparator;
                ss.str("");
            }
            
            for (std::vector<DagNode*>::iterator i = nodes.begin(); i != nodes.end(); ++i)
            {
                // get the node
                DagNode *node = *i;
                
                // print the value
                node->printValue(ss, prefixSeparator + suffixSeparator, int( columnWidth ), false, false, true);
                std::cout << prefixSeparator << ss.str() << suffixSeparator;
                ss.str("");
            }
            
            if ( printElapsedTime )
            {
                
                size_t timeUsed = time(NULL) - startTime;
                
                size_t hours   = timeUsed / 3600;
                size_t minutes = timeUsed / 60 - hours * 60;
                size_t seconds = timeUsed - minutes * 60 - hours * 3600;
                
                ss << std::setw( 2 ) << std::setfill( '0' ) << hours << ":";
                ss << std::setw( 2 ) << std::setfill( '0' ) << minutes << ":";
                ss << std::setw( 2 ) << std::setfill( '0' ) << seconds;
                
                std::cout << prefixSeparator << ss.str() << suffixSeparator;
                ss.str("");
                
            }
            
            if ( printWaitingTime )
            {

                if ( (gen-startGen) <= samplingFrequency )
                {
                    std::cout << prefixSeparator << "--:--:--" << suffixSeparator;
                }
                else
                {
                    double done = gen - startGen;
                    double total = numCycles - startGen;
                    double progress = done / total;
                    size_t timeUsed = time(NULL) - startTime;
                
                    size_t waitTime = double( timeUsed ) / progress - double( timeUsed );
                    
                    size_t hours   = waitTime / 3600;
                    size_t minutes = waitTime / 60 - hours * 60;
                    size_t seconds = waitTime - minutes * 60 - hours * 3600;
                
                    ss << std::setw( 2 ) << std::setfill( '0' ) << hours << ":";
                    ss << std::setw( 2 ) << std::setfill( '0' ) << minutes << ":";
                    ss << std::setw( 2 ) << std::setfill( '0' ) << seconds;
                
                    std::cout << prefixSeparator << ss.str() << suffixSeparator;
                }
            
            }
        
            std::cout << std::endl;
            std::cout.flush();
        
        }
        
    } // end if enabled
    
    currentGen = gen;
}


void ScreenMonitor::printHeader( void )
{
    if ( enabled == true )
    {
        // print empty line first
        std::cout << std::endl;
    
        // print everything to a string stream
        std::stringstream ss;

        // print one column for the iteration number
        std::string header = "Iter";

        int width = 9;
    
        int numWidth = int( log10( numCycles ) ) + 1;
        width = width > numWidth ? width : numWidth;
    
        int columnWidth = 12;

        StringUtilities::fillWithSpaces( header, width, true );
        ss << header << suffixSeparator;
    
        if ( posterior )
        {
            header = "Posterior";
            StringUtilities::fillWithSpaces( header, columnWidth, false );
            ss << prefixSeparator << header << suffixSeparator;
        }
    
        if ( likelihood )
        {
            header = "Likelihood";
            StringUtilities::fillWithSpaces( header, columnWidth, false );
            ss << prefixSeparator << header << suffixSeparator;
        }
    
        if ( prior )
        {
            header = "Prior";
            StringUtilities::fillWithSpaces( header, columnWidth, false );
            ss << prefixSeparator << header << suffixSeparator;
        }
    
        for (std::vector<DagNode *>::const_iterator it=nodes.begin(); it!=nodes.end(); it++)
        {
            const DagNode* the_node = *it;

            // print the header
            if ( the_node->getName() != "" )
            {
                ss << prefixSeparator;
                the_node->printName( ss, prefixSeparator + suffixSeparator, int( columnWidth ), false );
                ss << suffixSeparator;
            }
            else
            {
                header = "Unnamed";
                StringUtilities::fillWithSpaces( header, columnWidth, false );
                ss << prefixSeparator << header << suffixSeparator;
            }
        }
        
        if ( printElapsedTime )
        {
            // We know it takes 8 characters to print the waiting time, so hard-set column
            // width to this number
            header = "elapsed";
            StringUtilities::fillWithSpaces( header, 8, false );
            ss << prefixSeparator << header << suffixSeparator;
        }
    
        if ( printWaitingTime )
        {
            // We know it takes 8 characters to print the waiting time, so hard-set column
            // width to this number
            header = "ETA";
            StringUtilities::fillWithSpaces( header, 8, false );
            ss << prefixSeparator << header << suffixSeparator;
        }

        std::cout << ss.str() << std::endl;

        for (size_t i=0; i<ss.str().size(); ++i)
        {
            std::cout << "-";
        }
    
        std::cout << std::endl;
        
    } // end if enabled
    
}


void ScreenMonitor::reset( size_t n )
{
    startGen = currentGen;
    numCycles = n;
    startTime = time( NULL );
    printWaitingTime = true;
}



/**
 * Set the replicate index.
 * If the index is not 1, then the monitor will be disabled.
 */
void ScreenMonitor::setReplicateIndex(size_t idx)
{
    
    // store the index for possible later uses
    replicateIndex = idx;
    
    enabled = replicateIndex == 1;
    
}



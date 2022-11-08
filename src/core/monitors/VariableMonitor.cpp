#include "VariableMonitor.h"

#include <fstream>
#include <string>

#include "DagNode.h"
#include "Model.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RbSettings.h"
#include "RbVersion.h"
#include "Cloneable.h"
#include "StringUtilities.h"

using namespace RevBayesCore;

/* Constructor */
VariableMonitor::VariableMonitor(DagNode *n, unsigned long g, const path &fname,
                                 const std::string &del, bool pp, bool l, bool pr, bool ap, bool wv) :
    AbstractFileMonitor(n,g,fname,ap,wv),
    posterior( pp ),
    prior( pr ),
    likelihood( l ),
    separator( del )
{
    
}


/* Constructor */
VariableMonitor::VariableMonitor(const std::vector<DagNode *> &n, unsigned long g, const path &fname, const std::string &del,
                                 bool pp, bool l, bool pr, bool ap, bool wv) :
    AbstractFileMonitor(n,g,fname,ap,wv),
    posterior( pp ),
    prior( pr ),
    likelihood( l ),
    separator( del )
{

}


/* Clone the object */
VariableMonitor* VariableMonitor::clone(void) const
{

    return new VariableMonitor(*this);
}

/**
 * Print header for monitored values
 */
void VariableMonitor::printHeader( void )
{

    if ( enabled == true )
    {

//    out_stream.open( working_file_name.c_str(), std::fstream::out | std::fstream::app);
        out_stream.seekg(0, std::ios::end);

        if ( write_version == true )
        {
            RbVersion version;
            out_stream << "#RevBayes version (" + version.getVersion() + ")\n";
            out_stream << "#Build from " + version.getGitBranch() + " (" + version.getGitCommit() + ") on " + version.getDate() + "\n";
        }

        // print one column for the iteration number
        out_stream << "Iteration";

        if ( posterior == true )
        {
            // add a separator before every new element
            out_stream << separator;
            out_stream << "Posterior";
        }

        if ( likelihood == true )
        {
            // add a separator before every new element
            out_stream << separator;
            out_stream << "Likelihood";
        }

        if ( prior == true )
        {
            // add a separator before every new element
            out_stream << separator;
            out_stream << "Prior";
        }

        // print the headers for the variables
        printFileHeader();

        out_stream << std::endl;

        out_stream.flush();

        //    out_stream.close();
    }

}

/**
 * Monitor
 */
void VariableMonitor::monitor(unsigned long gen)
{

    // get the printing frequency
    unsigned long samplingFrequency = printgen;

    if ( enabled == true && gen % samplingFrequency == 0 )
    {
//        out_stream.open( working_file_name.c_str(), std::fstream::out | std::fstream::app);
        out_stream.seekg(0, std::ios::end);

        // print the iteration number first
        out_stream << gen;
        
        std::streamsize previousPrecision = out_stream.precision();
        std::ios_base::fmtflags previousFlags = out_stream.flags();
        out_stream.precision(RbSettings::userSettings().getOutputPrecision());

        if ( posterior == true )
        {
            // add a separator before every new element
            out_stream << separator;

            const std::vector<DagNode*> &n = model->getDagNodes();
            double pp = 0.0;
            for (std::vector<DagNode*>::const_iterator it = n.begin(); it != n.end(); ++it)
            {
                pp += (*it)->getLnProbability();
            }
            out_stream << pp;
        }

        if ( likelihood == true )
        {
            // add a separator before every new element
            out_stream << separator;

            const std::vector<DagNode*> &n = model->getDagNodes();
            double pp = 0.0;
            for (std::vector<DagNode*>::const_iterator it = n.begin(); it != n.end(); ++it)
            {
                if ( (*it)->isClamped() )
                {
                    pp += (*it)->getLnProbability();
                }
            }
            out_stream << pp;
        }

        if ( prior == true )
        {
            // add a separator before every new element
            out_stream << separator;

            const std::vector<DagNode*> &n = model->getDagNodes();
            double pp = 0.0;
            for (std::vector<DagNode*>::const_iterator it = n.begin(); it != n.end(); ++it)
            {
                if ( !(*it)->isClamped() )
                {
                    pp += (*it)->getLnProbability();
                }
            }
            out_stream << pp;
        }
        
        out_stream.setf(previousFlags);
        out_stream.precision(previousPrecision);

        monitorVariables( gen );

        out_stream << std::endl;

        out_stream.flush();

    }

}

/**
 * Print additional header for monitored values
 */
void VariableMonitor::printFileHeader( void )
{

    for (std::vector<DagNode *>::const_iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
        // add a separator before every new element
        out_stream << separator;

        const DagNode* the_node = *it;

        // print the header
        if (the_node->getName() != "")
        {
            the_node->printName(out_stream,separator, -1, true, flatten);
        }
        else
        {
            out_stream << "Unnamed";
        }

    }

}

/**
 * Monitor value at generation gen
 */
void VariableMonitor::monitorVariables(unsigned long gen)
{

    for (std::vector<DagNode*>::iterator i = nodes.begin(); i != nodes.end(); ++i)
    {
        // add a separator before every new element
        out_stream << separator;

        // get the node
        DagNode *node = *i;

        // print the value
        node->printValue(out_stream, separator, -1, false, false, true, flatten);
    }

}

/**
 * Combine output for the monitor.
 * Overwrite this method for specialized behavior.
 */
void VariableMonitor::combineReplicates( size_t n_reps, MonteCarloAnalysisOptions::TraceCombinationTypes tc )
{

    if ( enabled == true )
    {

        std::fstream combined_output_stream;

        int sample_number = 0;

        // open the stream to the file
        combined_output_stream.open( filename.string(), std::fstream::out);
        combined_output_stream.close();
        combined_output_stream.open( filename.string(), std::fstream::in | std::fstream::out);

        if ( tc == MonteCarloAnalysisOptions::SEQUENTIAL )
        {
            for (size_t i=0; i<n_reps; ++i)
            {
                std::stringstream ss;
                ss << "_run_" << (i+1);
                std::string s = ss.str();
                path current_file_name = appendToStem(filename, s);

                std::ifstream current_input_stream( current_file_name.string() );

                if ( not current_input_stream )
                {
                    throw RbException()<<"Could not open file "<<current_file_name<<".";
                }

                std::string read_line = "";
                size_t lines_skipped = 0;
                size_t lines_to_skip = ( write_version == true ? 3 : 1 );
                while (std::getline(current_input_stream,read_line))
                {
                    ++lines_skipped;
                    if ( lines_skipped < lines_to_skip)
                    {
                        if ( i == 0 )
                        {
                            // write output
                            combined_output_stream << read_line;

                            // add a new line
                            combined_output_stream << std::endl;
                        }
                        continue;
                    }
                    else if ( lines_skipped == lines_to_skip )
                    {

                        if ( i == 0 )
                        {
                            std::vector<std::string> fields;
                            StringUtilities::stringSplit(read_line, separator, fields);

                            // write output
                            combined_output_stream << fields[0];

                            // add a separator before every new element
                            combined_output_stream << separator;
                            combined_output_stream << "Replicate_ID";

                            for (size_t j=1; j<fields.size(); ++j)
                            {
                                // add a separator before every new element
                                combined_output_stream << separator;

                                // write output
                                combined_output_stream << fields[j];
                            }

                            // add a new line
                            combined_output_stream << std::endl;
                        }
                        continue;

                    }

                    std::vector<std::string> fields;
                    StringUtilities::stringSplit(read_line, separator, fields);

                    // add the current sample number
                    combined_output_stream << sample_number;

                    // add a separator before every new element
                    combined_output_stream << separator;
                    combined_output_stream << i;

                    ++sample_number;
                    for (size_t j=1; j<fields.size(); ++j)
                    {
                        // add a separator before every new element
                        combined_output_stream << separator;

                        // write output
                        combined_output_stream << fields[j];
                    }
                    // add a new line
                    combined_output_stream << std::endl;

                }

                current_input_stream.close();

            }

        }
        else if ( tc == MonteCarloAnalysisOptions::MIXED )
        {
            std::vector< std::ifstream* > input_streams;

            for (size_t i=0; i<n_reps; ++i)
            {
                std::stringstream ss;
                ss << "_run_" << (i+1);
                std::string s = ss.str();

                path current_file_name = appendToStem(filename, s);

                std::ifstream * current_input_stream = new std::ifstream( current_file_name.string() );

                if ( not *current_input_stream )
                {
                    throw RbException()<<"Could not open file "<<current_file_name<<".";
                }

                input_streams.push_back( current_input_stream );
            }

            std::vector<std::string> read_lines(n_reps,"");
            size_t lines_skipped = 0;
            size_t lines_to_skip = ( write_version == true ? 3 : 1 );
            while ( std::getline(*input_streams[0],read_lines[0]) )
            {
                for (size_t i=1; i<n_reps; ++i)
                {
                    if ( !(std::getline(*input_streams[i],read_lines[i])) )
                    {
                        throw RbException("Cannot merge output trace files with unequal number of lines.");
                    }
                }

                ++lines_skipped;
                if ( lines_skipped < lines_to_skip)
                {
                    // write output
                    combined_output_stream << read_lines[0];

                    // add a new line
                    combined_output_stream << std::endl;
                    continue;
                }
                else if ( lines_skipped == lines_to_skip )
                {
                    std::vector<std::string> fields;
                    StringUtilities::stringSplit(read_lines[0], separator, fields);

                    // write output
                    combined_output_stream << fields[0];

                    // add a separator before every new element
                    combined_output_stream << separator;
                    combined_output_stream << "Replicate_ID";

                    for (size_t j=1; j<fields.size(); ++j)
                    {
                        // add a separator before every new element
                        combined_output_stream << separator;

                        // write output
                        combined_output_stream << fields[j];
                    }

                    // add a new line
                    combined_output_stream << std::endl;
                    continue;

                }


                for (size_t i=0; i<n_reps; ++i)
                {
                    std::vector<std::string> fields;
                    StringUtilities::stringSplit(read_lines[i], separator, fields);

                    // add the current sample number
                    combined_output_stream << sample_number;

                    // add a separator before every new element
                    combined_output_stream << separator;
                    combined_output_stream << i;

                    ++sample_number;
                    for (size_t j=1; j<fields.size(); ++j)
                    {
                        // add a separator before every new element
                        combined_output_stream << separator;

                        // write output
                        combined_output_stream << fields[j];
                    }

                    // add a new line
                    combined_output_stream << std::endl;

                }

            }

            for (size_t i=0; i<n_reps; ++i)
            {
                input_streams[i]->close();
                delete input_streams[i];
            }

        }

        combined_output_stream.close();

    }

}

/**
 * Set flag about whether to print the likelihood.
 *
 * \param[in]   tf   Flag if the likelihood should be printed.
 */
void VariableMonitor::setPrintLikelihood(bool tf)
{

    likelihood = tf;

}


/**
 * Set flag about whether to print the posterior probability.
 *
 * \param[in]   tf   Flag if the posterior probability should be printed.
 */
void VariableMonitor::setPrintPosterior(bool tf)
{

    posterior = tf;

}


/**
 * Set flag about whether to print the prior probability.
 *
 * \param[in]   tf   Flag if the prior probability should be printed.
 */
void VariableMonitor::setPrintPrior(bool tf)
{

    prior = tf;

}

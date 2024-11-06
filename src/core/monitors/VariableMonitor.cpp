#include "VariableMonitor.h"

#include <cstdint>
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
VariableMonitor::VariableMonitor(DagNode *n, std::uint64_t g, const path &fname,
                                 const SampleFormat& f, bool pp, bool l, bool pr, bool ap, bool wv) :
    AbstractFileMonitor(n,g,fname,ap,wv),
    posterior( pp ),
    prior( pr ),
    likelihood( l ),
    format( f )
{
    
}


/* Constructor */
VariableMonitor::VariableMonitor(const std::vector<DagNode *> &n, std::uint64_t g, const path &fname, const SampleFormat &f,
                                 bool pp, bool l, bool pr, bool ap, bool wv) :
    AbstractFileMonitor(n,g,fname,ap,wv),
    posterior( pp ),
    prior( pr ),
    likelihood( l ),
    format( f )
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
    if (not enabled) return;
    
    if (to<JSONFormat>(format))
    {
	std::vector<json> fields;
	// print one column for the iteration number
	fields.push_back("Iteration");

	if ( posterior ) fields.push_back("Posterior");

	if ( likelihood ) fields.push_back("Likelihood");

	if ( prior ) fields.push_back("Prior");

	json header;
	header["fields"] = fields;
	header["format"] = "MCON";
	header["version"] = "0.1";
	header["nested"] = true;
	header["atomic"] = false;
    
	if ( write_version == true )
	{
	    RbVersion version;

	    json rb;
	    rb["version"] = version.getVersion();
	    rb["branch"] = version.getGitBranch();
	    rb["commit"] = version.getGitCommit();
	    rb["builddate"] = version.getDate();
	    // build date?
	    header["RevBayes"] = rb;
	}

	out_stream.seekg(0, std::ios::end);
	out_stream<<header<<std::endl;
    }
    else if (auto f = to<SeparatorFormat>(format))
    {
	auto& separator = f->separator;
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
	    out_stream << separator << "Posterior";
	}

	if ( likelihood == true )
	{
	    // add a separator before every new element
	    out_stream << separator << "Likelihood";
	}

	if ( prior == true )
	{
	    // add a separator before every new element
	    out_stream << separator << "Prior";
	}

	// print the headers for the variables
	printFileHeader();

	out_stream << std::endl;
    }

    out_stream.flush();
}

/**
 * Monitor
 */
void VariableMonitor::monitor(std::uint64_t gen)
{
    if ( not enabled or gen % printgen != 0 ) return;

    out_stream.seekg(0, std::ios::end);

    if (to<SeparatorFormat>(format))
    {
        // print the iteration number "first", before we change the precision?
        out_stream << gen;
    }

    std::streamsize previousPrecision = out_stream.precision();
    std::ios_base::fmtflags previousFlags = out_stream.flags();
    out_stream.precision(RbSettings::userSettings().getOutputPrecision());

    double Posterior = 0;
    double Likelihood = 0;
    double Prior = 0;
    if (posterior or likelihood or prior)
    {
	for (auto& node: model->getDagNodes())
	{
	    double Pr = node->getLnProbability();
	    Posterior += Pr;
	    if (node->isClamped())
		Likelihood += Pr;
	    else
		Prior += Pr;
	}
    }
        

    if (to<JSONFormat>(format))
    {
	json line;

	line["Iteration"] = gen;
        
	if (Posterior) line["Posterior"] = Posterior;
	if (Likelihood) line["Likelihood"] = Likelihood;
	if (Prior) line["Prior"] = Prior;

	for (auto& node: nodes)
	{
	    auto name = node->getName();
	    if (name.empty())
		name = std::to_string((uintptr_t)node);
	    line[name] = node->getValueAsJSON();
	}

	out_stream << line << "\n";
    }
    else
    {
	auto& separator = to<SeparatorFormat>(format)->separator;

        if ( posterior == true )
        {
            // add a separator before every new element
            out_stream << separator << Posterior;
        }

        if ( likelihood == true )
        {
            // add a separator before every new element
            out_stream << separator << Likelihood;
        }

        if ( prior == true )
        {
            // add a separator before every new element
            out_stream << separator << Prior;
        }
        
        monitorVariables( gen );

        out_stream << std::endl;
    }

    out_stream.setf(previousFlags);
    out_stream.precision(previousPrecision);
    out_stream.flush();
}

/**
 * Print additional header for monitored values
 */
void VariableMonitor::printFileHeader( void )
{
    auto& separator = to<SeparatorFormat>(format)->separator;

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
void VariableMonitor::monitorVariables(std::uint64_t gen)
{
    auto& separator = to<SeparatorFormat>(format)->separator;

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
    if ( enabled == true and to<SeparatorFormat>(format))
    {

	auto & separator = to<SeparatorFormat>(format)->separator;

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

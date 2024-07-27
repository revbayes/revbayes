#include <cstdlib>
#include <boost/program_options.hpp> // IWYU pragma: keep
#include <string>
#include <vector>
#include <iostream>

namespace po = boost::program_options;
using po::variables_map;

#include "MpiUtilities.h"
#include "Parser.h"
#include "RbVersion.h"
#include "RbException.h"
#include "RbSettings.h"
#include "RevClient.h"
#include "RevLanguageMain.h"
#include "RlCommandLineOutputStream.h"
#include "RlUserInterface.h"
#include "StringUtilities.h"
#include "TerminalFormatter.h"
#include "Workspace.h"

#ifdef RB_MPI
#include <mpi.h>
#include "RandomNumberFactory.h" // IWYU pragma: keep
#include "RandomNumberGenerator.h" // IWYU pragma: keep
#endif

std::string usage()
{
    return "Usage: rb [OPTIONS]\n       rb [OPTIONS] <file1> [<file2> ...]";
    // Other usages not mentioned
}


std::string short_description()
{
    return "Bayesian phylogenetic inference using probabilistic graphical models and an interpreted language";
}

//
variables_map parse_cmd_line(int argc, char* argv[])
{
    
    using namespace po;

    // Put all options in one group for now.
    options_description general("Options");
    general.add_options()
	("help,h","Show information on flags.")
	("version,v","Show version and exit.")

	// implicit_value(1) means that -V => -V1
//      RevBayes doesn't use a global verbose_logging flag.
//	("verbose,V",value<int>()->implicit_value(1),"Log extra information for debugging.")

	("batch,b","Run in batch mode.")
        ("jupyter,j","Run in jupyter mode.")
        // multitoken means that `--args a1 a2 a3` works the same as `--args a1 --args a2 --args a3`
        ("args",value<std::vector<std::string> >()->multitoken(),"Command line arguments to initialize RevBayes variables.")
        // multitoken means that `--args a1 a2 a3` works the same as `--args a1 --args a2 --args a3`
        ("cmd",value<std::vector<std::string> >()->multitoken(),"Script and command line arguments to initialize RevBayes variables.")
	// composing means that --file can occur multiple times
        ("file",value<std::vector<std::string> >()->composing(),"File(s) to source.")
        ("setOption",value<std::vector<std::string> >()->composing(),"Set an option key=value.")
	;

    // Treat all positional options as "file" options.
    positional_options_description p;
    p.add("file", -1);

    // Parse the command line into variables_map 'args'
    variables_map args;

    try {
        store(command_line_parser(argc, argv).options(general).positional(p).run(), args);
    }
    catch(po::error& e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        std::cout << std::endl;
        
        int rank = 0;
#ifdef RB_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
        if (rank == 0)
        {
            std::cout << usage() << std::endl;
            std::cout << short_description() << std::endl;
            std::cout << std::endl;
            std::cout << general << std::endl;
            std::cout << "See http://revbayes.github.io for more information." << std::endl;
        }
#ifdef RB_MPI
        MPI_Finalize();
#endif
        std::exit(0);
        
    }
    notify(args);

    // Print flags and usage info in this function since we know the flags here.
    if ( args.count("help") > 0 )
    {
	// Do we want to avoid displaying --file here, since its a positional option also?

        int rank = 0;
#ifdef RB_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
        if (rank == 0)
        {
            std::cout<< usage() << std::endl;
            std::cout<< short_description() << std::endl;
            std::cout << std::endl;
            std::cout<< general << std::endl;
            std::cout<< "See http://revbayes.github.io for more information." << std::endl;
        }
#ifdef RB_MPI
        MPI_Finalize();
#endif
        std::exit(0);
    }

    return args;
}

#ifndef MAC_GUI // don't include main for the Mac GUI written in Swift

int main(int argc, char* argv[]) {

    using namespace po;

    // If we aren't using MPI, this will be zero.
    // If we are using MPI, it will be zero for the first process.
    int process_id = 0;

#   ifdef RB_MPI
    int num_processes = 0;
    try
    {
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
        MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

        RevBayesCore::MpiUtilities::synchronizeRNG( MPI_COMM_WORLD );

    }
    catch (char* str)
    {
        return -1;
    }
#else
    RevBayesCore::MpiUtilities::synchronizeRNG(  );
#endif

    /* Parse argv to get the command line arguments */
    variables_map args = parse_cmd_line(argc, argv);

    if ( args.count("version") > 0 )
    {
        std::cout << RbVersion().getVersion() << std::endl;
        exit(0);
    }

    if ( args.count("verbose") > 0 )
    {
        int verbosity = args["verbose"].as<int>();
    }

    if ( args.count("setOption") > 0 )
    {
        std::vector<std::string> options = args["setOption"].as<std::vector<std::string> >();
        for(int i=0;i<options.size();i++)
        {
            std::vector<std::string> tokens;
            StringUtilities::stringSplit(options[i], "=", tokens);
            if (tokens.size() != 2)
            {
                throw RbException("Option '"+options[i]+"' must have the form key=value");
            }
            else
            {
                RbSettings::userSettings().setOption( tokens[0], tokens[1], false );
            }
        }
    }
    
    /*default to interactive mode*/
    bool batch_mode = (args.count("batch") > 0);
    // FIXME -- the batch_mode variable appears to have no effect if true.

    /* seek out files from command line */
    std::vector<std::string> source_files;
    if ( args.count("file") > 0 )
    {
        source_files = args["file"].as<std::vector<std::string> >();
    }
    
    if ( args.count("args") && args.count("cmd"))
    {
        throw RbException("command line: received both --args and --cmd");
    }
    
    std::vector<std::string> rb_args;
    if ( args.count("args") > 0 )
    {
        rb_args = args["args"].as<std::vector<std::string> >();
    }
    else if ( args.count("cmd") > 0)
    {
        rb_args = args["cmd"].as<std::vector<std::string> >();
        source_files.push_back(rb_args[0]);
        rb_args.erase(rb_args.begin());

        // Let's make batch mode default to true for scripts.
        if ( args.count("batch") == 0 )
        {
            batch_mode = true;
        }
    }

    /* initialize environment */
    RevLanguageMain rl = RevLanguageMain(batch_mode);

    CommandLineOutputStream *rev_output = new CommandLineOutputStream();
    RevLanguage::UserInterface::userInterface().setOutputStream( rev_output );
    rl.startRevLanguageEnvironment(rb_args, source_files);

    if (args.count("jupyter"))
    {
        RevClient::startJupyterInterpreter();
    }
    else
    {
        enableTermAnsi();

        RevClient::startInterpreter();
    }

#   ifdef RB_MPI
    MPI_Finalize();
#   endif

    return 0;
}

#endif

#include <stdio.h>
#include <string>
#include <vector>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using po::variables_map;

#include "RbVersion.h"
#include "RbException.h"
#include "RbSettings.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RevClient.h"
#include "RevLanguageMain.h"
#include "RlCommandLineOutputStream.h"
#include "RlUserInterfaceOutputStream.h"
#include "RlUserInterface.h"
#include "Parser.h"
#include "Workspace.h"

#ifdef RB_MPI
#include <mpi.h>
#endif


using std::string;
using std::vector;

string usage()
{
    return "Usage: rb [OPTIONS]\n       rb [OPTIONS] <file1> [<file2> ...]";
    // Other usages not mentioned
}

string short_description()
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
	// composing means that --file can occur multiple times
	("file",value<vector<string> >()->composing(),"File(s) to source.")
	("setOption",value<vector<string> >()->composing(),"Set an option key=value.")
	;

    // Treat all positional options as "file" options.
    positional_options_description p;
    p.add("file", -1);

    // Parse the command line into variables_map 'args'
    variables_map args;

    store(command_line_parser(argc, argv).options(general).positional(p).run(), args);
    notify(args);

    // Print flags and usage info in this function since we know the flags here.
    if (args.count("help"))
    {
	// Do we want to avoid displaying --file here, since its a positional option also?

        int rank = 0;
#ifdef RB_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
        if (rank == 0)
        {
            std::cout<<usage()<<"\n";
            std::cout<<short_description()<<"\n\n";
            std::cout<<general<<"\n";
            std::cout<<"See http://revbayes.github.io for more information\n";
        }
#ifdef RB_MPI
        MPI_Finalize();
#endif
	std::exit(0);
    }

    return args;
}

#ifndef MAC_GUI

int main(int argc, char* argv[]) {

      using namespace po;
#   ifdef RB_MPI
    int process_id = 0;
    int num_processes = 0;
    try
    {
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
        MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

        unsigned int seed = 0;

        // sync the random number generators
        if ( process_id == 0 )
        {
            seed = RevBayesCore::GLOBAL_RNG->getSeed();

        }

        MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);

        RevBayesCore::GLOBAL_RNG->setSeed( seed );

    }
    catch (char* str)
    {
        return -1;
    }
#endif

    /* Parse argv to get the command line arguments */
    variables_map args = parse_cmd_line(argc, argv);

    if (args.count("version"))
    {
        std::cout<<RbVersion().getVersion()<<std::endl;
        exit(0);
    }

    if (args.count("verbose"))
    {
        int verbosity = args["verbose"].as<int>();
    }

    if (args.count("setOption"))
    {
        vector<string> options = args["setOption"].as<vector<string> >();
        for(int i=0;i<options.size();i++)
        {
            vector<string> tokens;
            StringUtilities::stringSplit(options[i], "=", tokens);
            if (tokens.size() != 2)
                    throw RbException("Option '"+options[i]+"' must have the form key=value");
            else
                RbSettings::userSettings().setOption( tokens[0], tokens[1], false );
        }
    }

    /*default to interactive mode*/
    bool batch_mode = args.count("batch");
    // FIXME -- the batch_mode variable appears to have no effect if true.

    /* seek out files from command line */
    vector<string> sourceFiles;
    if (args.count("file"))
	sourceFiles = args["file"].as<vector<string> >();

    /* initialize environment */
    RevLanguageMain rl = RevLanguageMain(batch_mode);

    CommandLineOutputStream *rev_output = new CommandLineOutputStream();
    RevLanguage::UserInterface::userInterface().setOutputStream( rev_output );
    rl.startRevLanguageEnvironment(sourceFiles);

#   ifdef RB_XCODE

#   ifndef RB_MPI
    int process_id = 0;

#   endif
    /* Declare things we need */
    int result = 0;
    std::string commandLine = "";
    std::string line = "";

    for (;;)
    {
        /* Print prompt based on state after previous iteration */
        if ( process_id == 0 )
        {
            if (result == 0 || result == 2)
            {
                std::cout << "> ";
            }
            else
            {
                std::cout << "+ ";
            }

            /* Get the line */
            std::istream& retStream = std::getline(std::cin, line);
            if (!retStream)
            {
#               ifdef RB_MPI
                MPI::Finalize();
#               endif
                exit(0);
            }

            if (result == 0 || result == 2)
            {
                commandLine = line;
            }
            else if (result == 1)
            {
                commandLine += ";" + line;
            }
        }

        size_t bsz = commandLine.size();
#       ifdef RB_MPI
        MPI_Bcast(&bsz, 1, MPI_INT, 0, MPI_COMM_WORLD);
#       endif

        char * buffer = new char[bsz+1];
        buffer[bsz] = 0;
        for (int i = 0; i < bsz; i++)
        {
            buffer[i] = commandLine[i];
        }
#       ifdef RB_MPI
        MPI_Bcast(buffer, (int)bsz, MPI_CHAR, 0, MPI_COMM_WORLD);
#       endif

        std::string tmp = std::string( buffer );

        result = RevLanguage::Parser::getParser().processCommand(tmp, &RevLanguage::Workspace::userWorkspace());
    }

#   else

    RevClient c;
    c.startInterpretor();

#   endif

#   ifdef RB_MPI
    MPI_Finalize();
#   endif

    return 0;
}

#else

    // do initialization here for Mac GUI
#endif

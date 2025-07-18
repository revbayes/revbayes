#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>

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
#include "CLI11.hpp" // CLI11 command-line parsing library.

#ifdef RB_MPI
#include <mpi.h>
#include "RandomNumberFactory.h" // IWYU pragma: keep
#include "RandomNumberGenerator.h" // IWYU pragma: keep
#endif

std::string usage()
{
    return "Usage: rb [OPTIONS]\n       rb [OPTIONS] file [args]\n       rb [options] -e expr [-e expr2 ...] [args]\n";
    // Other usages not mentioned
}


std::string short_description()
{
    return "Bayesian phylogenetic inference using probabilistic graphical models and an interpreted language";
}

struct ParsedOptions
{
    bool help = false;
    bool version = false;
    bool batch = false;
    bool no_header = false;
    bool jupyter = false;
    std::vector<std::string> options;
    std::vector<std::string> expressions;
    std::optional<std::string> filename;
    std::vector<std::string> args;
};

void show_help(const CLI::App& app)
{
    int rank = 0;

#ifdef RB_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if (rank == 0)
    {
        std::cerr << short_description() << std::endl;
        std::cerr << usage() << std::endl;
//        std::cerr << std::endl;
        std::cerr << app.help() << std::endl;
        std::cerr << "See http://revbayes.github.io for more information." << std::endl;
    }
}

//
ParsedOptions parse_cmd_line(int argc, char* argv[])
{
    CLI::App stage1(short_description());

    ParsedOptions options;

//    stage1.add_option("-h,--help",    options.help,      "Show information on flags");
    
    stage1.add_flag("-v,--version",  options.version,    "Show version and exit");
    stage1.add_flag("-b,--batch",    options.batch,      "Run in batch mode");
    stage1.add_flag("--no-header",   options.no_header,  "Suppress header");
    stage1.add_flag("-j,--jupyter",  options.jupyter,    "Run in jupyter mode");

    stage1.add_option("--setOption",  options.options,   "Set an option key=value.  See ?setOption for the list of available keys and their associated values.");

    try {
        stage1.parse(argc, argv);
    }
    catch (const CLI::ParseError &e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << std::endl;
        
        show_help(stage1);

#ifdef RB_MPI
        MPI_Finalize();
#endif

        std::exit(e.get_exit_code());
    }

    // Print flags and usage info in this function since we know the flags here.
    if ( options.help )
    {
	// Do we want to avoid displaying --file here, since its a positional option also?

        show_help(stage1);

#ifdef RB_MPI
        MPI_Finalize();
#endif

        std::exit(0);
    }

    
    
    return options;
}

#ifndef MAC_GUI // don't include main for the Mac GUI written in Swift

int main(int argc, char* argv[])
{
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
    auto options = parse_cmd_line(argc, argv);

    if ( options.version )
    {
        std::cout << RbVersion().getVersion() << std::endl;
        exit(0);
    }

    for(auto& option: options.options)
    {
        std::vector<std::string> tokens;
        StringUtilities::stringSplit(option, "=", tokens);
        if (tokens.size() != 2)
        {
            throw RbException() << "Option '" << option << "' must have the form key=value"; 
        }
        else
        {
            RbSettings::userSettings().setOption( tokens[0], tokens[1], false );
        }
    }

    /*default to interactive mode*/
    bool batch_mode = options.batch;
    // FIXME -- the batch_mode variable appears to have no effect if true.

    /* seek out files from command line */
    std::vector<std::string> source_files;
    if ( options.filename )
    {
        source_files = { *options.filename };
    }

    std::vector<std::string> rb_args = options.args;;

    /* initialize environment */
    bool show_header = not options.no_header;
    RevLanguageMain rl = RevLanguageMain(batch_mode, show_header);

    CommandLineOutputStream *rev_output = new CommandLineOutputStream();
    RevLanguage::UserInterface::userInterface().setOutputStream( rev_output );
    rl.startRevLanguageEnvironment(rb_args, source_files);

    if ( options.jupyter )
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

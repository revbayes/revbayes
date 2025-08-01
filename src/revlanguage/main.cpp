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
#include "RandomNumberGenerator.h" // for setting the seed.
#include "RandomNumberFactory.h"

#ifdef RB_MPI
#include <mpi.h>
#include "RandomNumberFactory.h" // IWYU pragma: keep
#include "RandomNumberGenerator.h" // IWYU pragma: keep
#endif

std::string usage()
{
    return "Usage: rb [options]\n       rb [options] file [args]\n       rb [options] -e expr [-e expr2 ...] [args]\n";
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
    bool error_exit = false;
    bool echo = true;
    bool interactive = false;
    bool quiet = false;
    bool jupyter = false;
    std::vector<std::string> options;

    std::optional<std::uint64_t> seed;

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
        std::cerr << app.help() << std::endl;
    }
}

class RBFormatter : public CLI::Formatter {
  public:
    std::string make_usage(const CLI::App *app, std::string name) const override
        {
            return usage();
        }
    std::string make_footer(const CLI::App *app) const override
        {
            return "\nExpressions (one or more '-e <expr>') may be used *instead* of 'file'.\nSee http://revbayes.github.io for more information.";
        }
};

ParsedOptions parse_cmd_line(int argc, char* argv[])
{
    ParsedOptions options;

    // Stage 1: Parse options until we see something we don't recognize.
    CLI::App stage1(short_description());
    stage1.formatter(std::make_shared<RBFormatter>());
    stage1.get_formatter()->column_width(35);
    stage1.get_formatter()->right_column_width(45);

    stage1.add_flag("-v,--version",          options.version,   "Show version and exit");
    std::optional<bool> interactive;
    stage1.add_flag("-i,--interactive",      interactive,       "Force interactive with expressions or file");
    stage1.add_flag("-q,--quiet",            options.quiet,     "Don't print startup message");
    stage1.add_flag("-b,--batch",                               "Deprecated");
    stage1.add_flag("-j,--jupyter",          options.jupyter,   "Run in jupyter mode");

    std::optional<bool> error_exit;
    stage1.add_option("-x,--error-exit",     error_exit,        "Exit on the first error");
    std::optional<bool> echo;
    stage1.add_option("-p,--echo",           echo,              "Echo commands to the screen");
    stage1.add_option("-s,--seed",           options.seed,      "Random seed");

    stage1.add_option("-o,--setOption",  options.options,   "Set an option key=value  (See ?setOption for the list of available keys and their associated values)")->allow_extra_args(false);

    try {
        stage1.parse(argc, argv);
    }
    catch (const CLI::ExtrasError &e)
    {
        // Extra arguments are fine.
    }
    catch (const CLI::CallForHelp &e)
    {
        show_help(stage1);

#ifdef RB_MPI
        MPI_Finalize();
#endif

        std::exit(0);
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

    {
        auto remaining = stage1.remaining();

        // Stage 2: Parse a sequence of "-e expr" until we find something that isn't "-e".
        int i = 0;
        for(; i<remaining.size() and remaining[i] == "-e"; i+=2)
        {
            if (i+1 < remaining.size())
                options.expressions.push_back(remaining[i+1]);
            else
            {
                std::cerr<<"-e not followed by expression\n";
                std::exit(1);
            }
        }

        // Stage 3: Put everything else into the args vector.
        for(; i<remaining.size(); i++)
            options.args.push_back(remaining[i]);
    }

    // Fixup: if there are no expressions, the first positional argument is a filename.
    if (options.expressions.empty() and not options.args.empty())
    {
        if (options.args[0].size() > 0 and options.args[0][0] == '-')
        {
            std::cerr<<"Error: unrecognized option '"<<options.args[0]<<"'\n";
            std::cerr<<"See --help for more info.\n";
            std::exit(1);
        }
        options.filename = options.args[0];
        options.args.erase(options.args.begin());
    }

    options.interactive = interactive ? (interactive.value()) : (options.expressions.empty() and not options.filename);

    options.error_exit = error_exit ? (error_exit.value()) : not options.interactive;

    options.echo = echo ? (echo.value()) : options.interactive;
    
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
    auto cmd_line = parse_cmd_line(argc, argv);

    if ( cmd_line.version )
    {
        std::cout << RbVersion().getVersion() << std::endl;
        exit(0);
    }

    /* Set default session properties from cmd line flags */
    auto& settings = RbSettings::userSettings();
    settings.setInteractive( cmd_line.interactive );
    settings.setEcho( cmd_line.echo );
    settings.setErrorExit( cmd_line.error_exit );

    /* Set user options from cmd line */
    for(auto& option: cmd_line.options)
    {
        std::vector<std::string> tokens;
        StringUtilities::stringSplit(option, "=", tokens);
        if (tokens.size() != 2)
        {
            std::cerr << "Option '" << option << "' must have the form key=value\n";
            std::exit(1);
        }
        else
        {
            settings.setOption( tokens[0], tokens[1], false );
        }
    }

    /* Set seed from command line */
    if (cmd_line.seed)
    {
        RevBayesCore::RandomNumberGenerator *rng = RevBayesCore::GLOBAL_RNG;
        rng->setSeed( cmd_line.seed.value() );
    }

    /* initialize environment */
    RevLanguageMain rl = RevLanguageMain(cmd_line.quiet);

    /* Set output stream */
    CommandLineOutputStream *rev_output = new CommandLineOutputStream();
    RevLanguage::UserInterface::userInterface().setOutputStream( rev_output );

    /* Any script or expressions from `-e expr` are executed here. */
    rl.startRevLanguageEnvironment(cmd_line.expressions, cmd_line.filename, cmd_line.args);

    /* Interactive or jupyter commands are executed here. */
    if ( cmd_line.jupyter )
    {
        RevClient::startJupyterInterpreter();
    }
    else if ( settings.getInteractive() )
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

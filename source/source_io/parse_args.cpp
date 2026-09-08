#include "parse_args.h"
#include "build_info.h"
#include "input_help.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>

#include "module_parameter/read_input.h"
#include "source_base/version.h"

#if defined(COMMIT_INFO)
#include "commit.h"
#endif

namespace ModuleIO
{

/**
 * @brief Format description for brief display in search results.
 *
 * Converts multi-line descriptions with list markers to single-line format:
 * - Replaces newlines with spaces
 * - Converts * list markers to -
 * - Collapses multiple spaces
 * - Truncates to max_length characters
 */
static std::string format_brief_description(const std::string& desc, size_t max_length = 60) {
    std::string result;
    result.reserve(desc.length());

    bool prev_space = false;
    bool at_line_start = true;

    for (size_t i = 0; i < desc.length(); ++i) {
        char c = desc[i];

        if (c == '\n') {
            // Replace newline with space (collapse multiple)
            if (!prev_space && !result.empty()) {
                result += ' ';
                prev_space = true;
            }
            at_line_start = true;
        } else if (c == '*' && at_line_start) {
            // Convert * list marker to -
            result += '-';
            prev_space = false;
            at_line_start = false;
        } else if (c == ' ' || c == '\t') {
            if (!prev_space && !result.empty()) {
                result += ' ';
                prev_space = true;
            }
            // Don't clear at_line_start for leading spaces
        } else {
            result += c;
            prev_space = false;
            at_line_start = false;
        }
    }

    // Trim trailing spaces
    while (!result.empty() && result.back() == ' ') {
        result.pop_back();
    }

    // Truncate if needed
    if (result.length() > max_length) {
        result = result.substr(0, max_length) + "...";
    }

    return result;
}

void print_build_info()
{
    std::cout << "ABACUS H0 build information\n";
    std::cout << "Version: " << VERSION << '\n';
#if defined(COMMIT)
    std::cout << "Git commit: " << COMMIT << '\n';
#else
    std::cout << "Git commit: unavailable\n";
#endif
    std::cout << "Build type: " << ABACUS_BUILD_TYPE << '\n';
    std::cout << "C++ compiler: " << ABACUS_CXX_COMPILER_ID << ' '
              << ABACUS_CXX_COMPILER_VERSION << '\n';
    std::cout << "MPI: " << ABACUS_MPI_STATUS << '\n';
    std::cout << "OpenMP: " << ABACUS_OPENMP_STATUS << '\n';
    std::cout << "Basis: LCAO\nXC: built-in LDA/GGA\n";
}

void parse_args(int argc, char** argv)
{
    for (int i = 1; i < argc; ++i) // Start from 1 to skip the program name
    {
        std::string arg = argv[i];
        if (arg == "--version" || arg == "-v" || arg == "-V")
        {
#if defined(VERSION)
            const char* version = VERSION;
#else
            const char* version = "unknown";
#endif
            std::cout << "ABACUS version " << version << std::endl;
            std::exit(0);
        }
        else if (arg == "--info" || arg == "-i" || arg == "-I")
        {
            print_build_info();
            std::exit(0);
        }
        else if (arg == "-h" || arg == "--help")
        {
            // Handle -h or -h <key>
            if (i + 1 < argc) {
                // Next argument exists - check if it's a parameter key
                std::string next_arg = argv[i + 1];
                if (!next_arg.empty() && next_arg[0] != '-') {
                    // Not another flag - treat as parameter key
                    ParameterHelp::initialize();
                    // When parameter is found, output to stdout; when not found, output to stderr
                    if (!ParameterHelp::show_parameter_help(next_arg, std::cout)) {
                        std::cerr << "Error: Unknown parameter '" << next_arg << "'" << std::endl;

                        // Try to find similar parameters (fuzzy matching)
                        auto suggestions = ParameterHelp::find_similar_parameters(next_arg, 5, 3);
                        if (!suggestions.empty()) {
                            std::cerr << "\nDid you mean one of these?" << std::endl;
                            for (const auto& suggestion : suggestions) {
                                std::cerr << "  - " << suggestion << std::endl;
                            }
                        }

                        std::cerr << "\nUse 'abacus -s <keyword>' to search for parameters." << std::endl;
                        std::exit(1);
                    }
                    std::exit(0);
                }
            }
            // No argument or next is a flag - show general help
            ParameterHelp::show_general_help();
            std::exit(0);
        }
        else if (arg == "-s" || arg == "--search")
        {
            // Require search query
            if (i + 1 >= argc || argv[i + 1][0] == '-') {
                std::cerr << "Error: -s requires a search query" << std::endl;
                std::exit(1);
            }

            std::string query = argv[++i];

            // Initialize help system
            ParameterHelp::initialize();

            auto results = ParameterHelp::search_parameters(query);
            if (results.empty()) {
                std::cerr << "No parameters found matching '" << query << "'" << std::endl;
                std::exit(1);
            }

            // Display results
            std::cout << "\nFound " << results.size() << " parameter(s) matching '" << query << "':\n\n";
            for (const auto& param : results) {
                auto metadata = ParameterHelp::get_metadata(param);
                std::cout << "  " << std::left << std::setw(30) << param;
                if (!metadata.name.empty() && !metadata.description.empty()) {
                    // Format description for brief display
                    std::string desc = format_brief_description(metadata.description, 60);
                    std::cout << " - " << desc;
                }
                std::cout << std::endl;
            }
            std::cout << "\nUse 'abacus -h <parameter>' for detailed help." << std::endl;
            std::exit(0);
        }
        else if (arg == "--generate-parameters-yaml")
        {
            ParameterHelp::generate_yaml(std::cout);
            std::exit(0);
        }
        else if (arg == "--check-input")
        {
            ModuleIO::ReadInput::check_mode = true;
        }
        else
        {
            // Error message goes to stderr
            std::cerr << "Error: Unknown argument: " << arg << std::endl;
            std::cerr << std::endl;

            // Display help information to stderr
            ParameterHelp::show_general_help(std::cerr);
            std::exit(1);
        }
    }
}

} // namespace ModuleIO

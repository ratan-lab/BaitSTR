#ifndef CLPARSING_H_
#define CLPARSING_H_

#include <inttypes.h>

#include "utilities.h"
#include "sllist.h"

// This reflects a CLI interface for my applications.
// A normal use-case for this is as follows:
// a) The program creates a new structure to store all the CL arguments that are
//    options:
//      commandLineArguments* arguments = NewCommandLineArguments()
//    --help and --debug are include by default. Only options beginning with
//    -- are permitted. I do not support short options beginning with a "-".
//    A single dash (-) is considered to be an argument that is not an option.
//    A double dash (--) stops the command line processing (this behavior is
//    similar to that of getopt).
// b) All the valid options should be added to the db using the function
//    addOption:
//      AddOption(&args, optString, optDefaultValue,FALSE,TRUE,helpString,func);
//    The AddOption method has six parameters (besides the argument where all
//    these options have to be added). The first parameter is a
//    string that represents the option. The second parameter is a character
//    string which depicts the default value for the option (TRUE/FALSE for
//    boolean options). The third parameter is a boolean that
//    specifies whether the option requires an argument or not. In the case of a
//    boolean option (sometimes referred to as a flag) an argument value is not
//    present so false is passed. The fourth parameter is a boolean that
//    specifies whether this option is documented or not. An undocumented option
//    is not printed in the help string for the program. The next parameter is
//    the description of the option. This description will be used in the usage
//    text of the application. The last argument is a pointer to the function
//    that must be run every time the value of this option is set using the API
//    here.
// c) The function ParseOptions should be called now to read the options from
//    the command line. It modifies the argv array so that only parameters that
//    are not options remain in it. Those can then be parsed by the program
//    separately.
//

// structure that holds a single command line argument
typedef struct CommandLineArguments_st {
    struct CommandLineArguments_st* next;
    const char* option_string;
    Bool is_value_expected;
    Bool is_option_documented;
    const char* help_string;
    char* option_default_value;
    char* option_value;
    void (*func)();
}CommandLineArguments;

// create a new structure to hold the arguments.
CommandLineArguments* NewCommandLineArguments();

// add the following options to the database. These are the options that should
// permitted to be modified from the command line.
void AddOption(CommandLineArguments** pArguments,
               const char* const option_string,
               char* const default_option_value,
               const Bool is_argument_expected,
               const Bool is_option_documented,
               const char* const help_string,
               void (*sanity_check)());

// return TRUE if this is a valid option that can  be modified using the command
// line.
Bool HasOption(CommandLineArguments* const arguments,
               const char* const option_string);

// return the boolean (TRUE/FALSE) value associated with the option. This
// function will fail if the value associated with this option is not TRUE or
// FALSE.
Bool GetOptionBoolValueOrDie(CommandLineArguments* const arguments,
                             const char* const option_string);

// return the integer value associated with this option. This function fails
// if the value associated with this option is not an integer.
int GetOptionIntValueOrDie(CommandLineArguments* const arguments,
                           const char* const option_string);

// return the double value associated with this option. This function fails
// if the value associated with this option is not an double.
double GetOptionDoubleValueOrDie(CommandLineArguments* const arguments,
                                 const char* const option_string);

// return the unsigned integer associated with this option. This function fails
// if the value associated with this option is not an uint.
uint GetOptionUintValueOrDie(CommandLineArguments* const arguments,
                             const char* const option_string);

// return the integer64 value associated with this option. This function fails
// if the value associated with this option is not an integer.
uint64_t GetOptionUint64ValueOrDie(CommandLineArguments* const arguments,
                                   const char* const option_string);

// return the string value associated with this option. It can be NULL, which
// represents an empty string
char* GetOptionStringValue(CommandLineArguments* const arguments,
                           const char* const option_string);

// set the value of this option to the following string
Bool SetOptionValue(CommandLineArguments* const arguments,
                    const char* const option_string,
                    char* const option_value);

// the work horse which is resposible to parse the command line arguments and
// uses that to modify the values of the various options. It also modifies argc
// and argv such that they only contain the arguments that are not options.
void ParseOptions(CommandLineArguments** const pArguments,
                  int* const pargc,
                  char*** pargv);

// print the usage string associated with this program.
void PrintSimpleUsageString(const CommandLineArguments* const arguments);

// free the resources used by this data structure. This is primarily so that
// this does not show up in a leak-check.
void FreeParseOptions(CommandLineArguments** pArguments, char*** pargv);

#endif  // COMMAND_LINE_PARSING_H_

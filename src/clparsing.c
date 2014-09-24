#include "clparsing.h"

extern char* program_name;
extern char* program_description;
extern char* program_use;
extern char* program_version_major;
extern char* program_version_minor;
extern char* program_revision_date;

// create a new structure to hold the arguments.
CommandLineArguments* NewCommandLineArguments() {
    CommandLineArguments* cmds = NULL;

    // by default --help, --memdebug and --debug are to be included
    CommandLineArguments* cmd = CkalloczOrDie(sizeof(CommandLineArguments));
    cmd->option_string = "help";
    cmd->is_value_expected = FALSE;
    cmd->is_option_documented = TRUE;
    cmd->option_default_value = "FALSE";
    cmd->option_value = cmd->option_default_value;
    cmd->help_string = "Print this string and quit.";
    cmd->func = NULL;
    SllAddHead(&cmds, cmd);

    cmd = CkalloczOrDie(sizeof(CommandLineArguments));
    cmd->option_string = "debug";
    cmd->is_value_expected = FALSE;
    cmd->is_option_documented = TRUE;
    cmd->option_default_value = "FALSE";
    cmd->option_value = cmd->option_default_value;
    cmd->help_string = "Print extra debug information in this run.";
    cmd->func = NULL;
    SllAddHead(&cmds, cmd);

    cmd = CkalloczOrDie(sizeof(CommandLineArguments));
    cmd->option_string = "memDebug";
    cmd->is_value_expected = FALSE;
    cmd->is_option_documented = FALSE;
    cmd->option_default_value = "FALSE";
    cmd->option_value = cmd->option_default_value;
    cmd->help_string = "Print memory related debug information in this run.";
    cmd->func = NULL;
    SllAddHead(&cmds, cmd);

    return cmds;
}

// add the following options to the database. These are the options that should
// permitted to be modified from the command line.
void AddOption(CommandLineArguments** pArguments,
               const char* const option_string,
               char* const default_option_value,
               const Bool is_argument_expected,
               const Bool is_option_documented,
               const char* const help_string,
               void (*sanity_check)()) {
    CommandLineArguments* cmd = CkalloczOrDie(sizeof(CommandLineArguments));
    cmd->option_string = option_string;
    cmd->is_value_expected = is_argument_expected;
    cmd->is_option_documented = is_option_documented;
    cmd->help_string = help_string;
    if (sanity_check != NULL) sanity_check(default_option_value);
    cmd->option_default_value = default_option_value;
    cmd->option_value = cmd->option_default_value;
    cmd->func = sanity_check;
    SllAddHead(pArguments, cmd);
    *pArguments = cmd;
}

// return TRUE if this is a valid option that can  be modified using the command
// line.
Bool HasOption(CommandLineArguments* const arguments,
               const char* const option_string) {
    CommandLineArguments* iter;
    for (iter = arguments; iter ; iter = iter->next) {
        if (SameString(iter->option_string, option_string)) {
            return TRUE;
        }
    }

    return FALSE;
}

// return the string assocaiated with this option. Return NULL if this option is
// not found in the database.
static char* GetOptionValue(CommandLineArguments* const arguments,
                            const char* const option_string) {
    CommandLineArguments* iter;
    for (iter = arguments; iter ; iter = iter->next) {
        if (SameString(iter->option_string, option_string)) {
            return iter->option_value;
        }
    }

    return NULL;
}

// return the boolean (TRUE/FALSE) value associated with the option. This
// function will fail if the value associated with this option is not TRUE or
// FALSE.
Bool GetOptionBoolValueOrDie(CommandLineArguments* const arguments,
                             const char* const option_string) {
    char* ptr = GetOptionValue(arguments, option_string);
    if (ptr == NULL) {
        PrintMessageThenDie("unknown option %s", option_string);
    }

    ForceAssert(SameString(ptr, "TRUE") || SameString(ptr, "FALSE"));
    return SameString(ptr, "TRUE");
}

// return the unsigned integer associated with this option. This function fails
// if the value associated with this option is not an uint.
uint GetOptionUintValueOrDie(CommandLineArguments* const arguments,
                             const char* const option_string) {
    char* ptr = GetOptionValue(arguments, option_string);
    if (ptr == NULL) {
        PrintMessageThenDie("unknown option %s", option_string);
    }

    uint result;
    if (sscanf(ptr, "%u", &result) != 1) {
        PrintMessageThenDie("%s: Expected unsigned int, but did not find it\n",
        option_string);
    }
    return result;
}

// return the integer value associated with this option. This function fails
// if the value associated with this option is not an integer.
int GetOptionIntValueOrDie(CommandLineArguments* const arguments,
                           const char* const option_string) {
    char* ptr = GetOptionValue(arguments, option_string);
    if (ptr == NULL) {
        PrintMessageThenDie("unknown option %s", option_string);
    }

    int result;
    if (sscanf(ptr, "%d", &result) != 1) {
        PrintMessageThenDie("%s: Expected an integer, but did not find it\n",
        option_string);
    }
    return result;
}

// return the double value associated with this option. This function fails
// if the value associated with this option is not an double.
double GetOptionDoubleValueOrDie(CommandLineArguments* const arguments,
                                 const char* const option_string) {
    char* ptr = GetOptionValue(arguments, option_string);
    if (ptr == NULL) {
        PrintMessageThenDie("unknown option %s", option_string);
    }

    double result;
    if (sscanf(ptr, "%lf", &result) != 1) {
        PrintMessageThenDie("%s: Expected an integer, but did not find it\n",
        option_string);
    }
    return result;
}

// return the uint64 value associated with this option. This function fails
// if the value associated with this option is not an integer.
uint64_t GetOptionUint64ValueOrDie(CommandLineArguments* const arguments,
                                   const char* const option_string) {
    char* ptr = GetOptionValue(arguments, option_string);
    if (ptr == NULL) {
        PrintMessageThenDie("unknown option %s", option_string);
    }

    uint64_t result;
    if (sscanf(ptr, "%"PRIu64, &result) != 1) {
        PrintMessageThenDie("%s: Expected an integer, but did not find it\n",
        option_string);
    }
    return result;
}

// return the string value associated with this option. It can be NULL, which
// represents an empty string
char* GetOptionStringValue(CommandLineArguments* const arguments,
                           const char* const option_string) {
    return GetOptionValue(arguments, option_string);
}

// set the value of this option to the following string
Bool SetOptionValue(CommandLineArguments* const arguments,
                    const char* const option_string,
                    char* const option_value) {
    CommandLineArguments* iter;
    for (iter = arguments; iter ; iter = iter->next) {
        if (SameString(iter->option_string, option_string)) {
            if (iter->func != NULL) {
                iter->func(option_value);
            }
            iter->option_value = option_value;
            return TRUE;
        }
    }

    return FALSE;
}

// the work horse which is resposible to parse the command line arguments and
// uses that to modify the values of the various options. It also modifies argc
// and argv such that they only contain the arguments that are not options.
void ParseOptions(CommandLineArguments** const pArguments,
                  int* const pargc,
                  char*** pargv) {
    int argc = *pargc;
    char** argv = *pargv;
    CommandLineArguments* arguments = *pArguments;

    // this is the modified set of arguments which does not include the options
    char** argvCopy = CkalloczOrDie(argc * sizeof(char*));

    int idx, idxForCopy;
    for (idx = 0, idxForCopy = 0; idx < argc; idx++) {
        if (CompareNames(argv[idx], "--", 2)) {
            char* ptr = strchr(argv[idx], '=');
            if (ptr == NULL) {
                // this could be a boolean option
                if (CompareNames(argv[idx], "--no", 4)) {
                    // set this option to FALSE
                    if (SetOptionValue(arguments,
                                       argv[idx] + 4,
                                       "FALSE") == FALSE) {
                        PrintMessageThenDie("unknown option %s", argv[idx] + 4);
                    }
                } else {
                    // set this option to TRUE
                    if (SetOptionValue(arguments,
                                       argv[idx] + 2,
                                      "TRUE") == FALSE) {
                        PrintMessageThenDie("unknown option %s", argv[idx] + 2);
                    }
                }
            } else {
                // this is an option of the format key=value
                *ptr = 0;
                if (SetOptionValue(arguments,
                                   argv[idx] + 2,
                                   ptr+1) == FALSE) {
                    PrintMessageThenDie("unknown option %s", argv[idx] + 2);
                }
            }
        } else {
            // this is not an option
            argvCopy[idxForCopy++] = argv[idx];
        }
    }

    // lets make help and debug as the first options
    SllReverse(pArguments);

    *pargc = idxForCopy;
    *pargv = argvCopy;
}

// print the usage string associated with this program.
void PrintSimpleUsageString(const CommandLineArguments* const arguments) {
    printf("Program : %s (%s)\n", program_name, program_description);
    printf("Version : %s.%s released %s\n\n",
    program_version_major, program_version_minor, program_revision_date);

    printf("usage:\n");
    printf("\t%s\n\n", program_use);
    printf("options:\n");

    // find the size of the longest string in the first column
    const CommandLineArguments* iter;
    uint maxLength = 0;
    for (iter = arguments; iter; iter = iter->next) {
        if (iter->is_option_documented == TRUE &&
           strlen(iter->option_string) > maxLength) {
            maxLength = strlen(iter->option_string);
        }
    }

    for (iter = arguments; iter; iter = iter->next) {
        if (iter->is_option_documented == FALSE) continue;

        printf("\t%s", iter->option_string);
        int idx = maxLength - strlen(iter->option_string);
        while (idx > 0) {
            printf(" ");
            idx--;
        }
        printf("  %s", iter->help_string);
        if (iter->is_value_expected == TRUE) {
            printf("[--%s=%s]", iter->option_string, iter->option_default_value);
        } else {
            if (SameString(iter->option_default_value, "TRUE") == TRUE) {
                printf("[--%s]", iter->option_string);
            } else {
                assert(SameString(iter->option_default_value, "FALSE") == TRUE);
                printf("[--no%s]", iter->option_string);
            }
        }
        printf("\n");
    }
}

// free the resources used by this data structure. This is primarily so that
// this does not show up in a leak-check.
void FreeParseOptions(CommandLineArguments** pArguments, char*** pArgv) {
    SllFreeList(pArguments);
    Ckfree(*pArgv);
    *pArguments = NULL;
}

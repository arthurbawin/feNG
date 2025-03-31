// Option parser copied from the MFEM library.
// Credit goes to the MFEM authors (https://mfem.org).

#ifndef _FEOPTIONSPARSER_
#define _FEOPTIONSPARSER_

#include "feNG.h"
#include "feMessage.h"
#if defined(HAVE_PETSC)
  #include "petscsys.h"
#endif

class feOptionsParser
{
public:
  enum OptionType { INT, DOUBLE, STRING, ENABLE, DISABLE };

private:
  struct Option {
    OptionType type;
    void *var_ptr;
    const char *short_name;
    const char *long_name;
    const char *description;
    bool required;

    Option() = default;

    Option(OptionType _type, void *_var_ptr, const char *_short_name, const char *_long_name,
           const char *_description, bool req)
      : type(_type), var_ptr(_var_ptr), short_name(_short_name), long_name(_long_name),
        description(_description), required(req)
    {
    }
  };

  int argc;
  char **argv;
  std::vector<Option> options;
  std::vector<int> option_check;
  // error_type can be:
  //  0 - no error
  //  1 - print help message
  //  2 - unrecognized option at argv[error_idx]
  //  3 - missing argument for the last option argv[argc-1]
  //  4 - option with index error_idx is specified multiple times
  //  5 - invalid argument in argv[error_idx] for option in argv[error_idx-1]
  //  6 - required option with index error_idx is missing
  int error_type, error_idx;
  bool _ignoreUnrecognizedOptions;

  static void WriteValue(const Option &opt, std::ostream &out);

public:
  /* Add a boolean option.
   Enable/disable tags are used to set the bool to true/false
   respectively. */
  void addOption(bool *var, const char *enable_short_name, const char *enable_long_name,
                 const char *disable_short_name, const char *disable_long_name,
                 const char *description, bool required = false)
  {
    options.push_back(
      Option(ENABLE, var, enable_short_name, enable_long_name, description, required));
    options.push_back(
      Option(DISABLE, var, disable_short_name, disable_long_name, description, required));
  }

  /// Add an integer option.
  void addOption(int *var, const char *short_name, const char *long_name, const char *description,
                 bool required = false)
  {
    options.push_back(Option(INT, var, short_name, long_name, description, required));
  }

  /// Add a double option.
  void addOption(double *var, const char *short_name, const char *long_name,
                 const char *description, bool required = false)
  {
    options.push_back(Option(DOUBLE, var, short_name, long_name, description, required));
  }

  /// Add a string (char*) option.
  void addOption(const char **var, const char *short_name, const char *long_name,
                 const char *description, bool required = false)
  {
    options.push_back(Option(STRING, var, short_name, long_name, description, required));
  }

private:
  void Parse();
#if defined(HAVE_PETSC)
  int ParsePetsc();
#endif

public:
  // Parse the command line options and return a status.
  feStatus parse();

  // Print the error message
  void PrintError(std::ostream &out) const;

  // Print the help message
  void PrintHelp(std::ostream &out) const;

  void PrintUsage(std::ostream &out) const;

  /// Construct a command line option parser with '_argc' and '_argv'.
  feOptionsParser(int _argc, char *_argv[], bool ignoreUnrecognizedOptions = false)
   : argc(_argc), argv(_argv), _ignoreUnrecognizedOptions(ignoreUnrecognizedOptions)
  {
    error_type = error_idx = 0;
  }
};

#endif
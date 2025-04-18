// Copied from the MFEM library.
// Credit goes to the MFEM authors (https://mfem.org).

#include "feOptionsParser.h"

int isValidAsInt(char *s)
{
  if(s == NULL || *s == '\0') {
    return 0; // Empty string
  }

  if(*s == '+' || *s == '-') {
    ++s;
  }

  if(*s == '\0') {
    return 0; // sign character only
  }

  while(*s) {
    if(!isdigit(*s)) {
      return 0;
    }
    ++s;
  }

  return 1;
}

int isValidAsDouble(char *s)
{
  // A valid floating point number for atof using the "C" locale is formed by
  // - an optional sign character (+ or -),
  // - followed by a sequence of digits, optionally containing a decimal-point
  //   character (.),
  // - optionally followed by an exponent part (an e or E character followed by
  //   an optional sign and a sequence of digits).

  if(s == NULL || *s == '\0') {
    return 0; // Empty string
  }

  if(*s == '+' || *s == '-') {
    ++s;
  }

  if(*s == '\0') {
    return 0; // sign character only
  }

  while(*s) {
    if(!isdigit(*s)) {
      break;
    }
    ++s;
  }

  if(*s == '\0') {
    return 1; // s = "123"
  }

  if(*s == '.') {
    ++s;
    while(*s) {
      if(!isdigit(*s)) {
        break;
      }
      ++s;
    }
    if(*s == '\0') {
      return 1; // this is a fixed point double s = "123." or "123.45"
    }
  }

  if(*s == 'e' || *s == 'E') {
    ++s;
    return isValidAsInt(s);
  } else {
    return 0; // we have encounter a wrong character
  }
}

void feOptionsParser::Parse()
{
  option_check.resize(options.size(), 0);
  for(int i = 1; i < argc;) {
    if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
      // print help message
      error_type = 1;
      return;
    }

    for(size_t j = 0; true; j++) {
      if(j >= options.size()) {
        // unrecognized option
      #if defined(HAVE_PETSC)
        // Keep parsing to accept unlisted PETSc options (e.g. -ksp_view)
        i++;
        break;
      #else
        if(_ignoreUnrecognizedOptions) {
          // Keep parsing.
          // This is used e.g. when parsing several times
          // and selecting only a subset of the options each time.
          i++;
          break;
        }
        error_type = 2;
        error_idx = i;
        return;
      #endif
      }

      if(strcmp(argv[i], options[j].short_name) == 0 ||
         strcmp(argv[i], options[j].long_name) == 0) {
        OptionType type = options[j].type;

        if(option_check[j]) {
          error_type = 4;
          error_idx = j;
          return;
        }
        option_check[j] = 1;

        i++;
        if(type != ENABLE && type != DISABLE && i >= argc) {
          // missing argument
          error_type = 3;
          error_idx = j;
          return;
        }

        int isValid = 1;
        switch(options[j].type) {
          case INT:
            isValid = isValidAsInt(argv[i]);
            *(int *)(options[j].var_ptr) = atoi(argv[i++]);
            break;
          case DOUBLE:
            isValid = isValidAsDouble(argv[i]);
            *(double *)(options[j].var_ptr) = atof(argv[i++]);
            break;
          case STRING:
            *(const char **)(options[j].var_ptr) = argv[i++];
            break;
          case ENABLE:
            *(bool *)(options[j].var_ptr) = true;
            option_check[j + 1] = 1; // Do not allow the DISABLE Option
            break;
          case DISABLE:
            *(bool *)(options[j].var_ptr) = false;
            option_check[j - 1] = 1; // Do not allow the ENABLE Option
            break;
        }

        if(!isValid) {
          error_type = 5;
          error_idx = i;
          return;
        }

        break;
      }
    }
  }

  // check for missing required options
  for(size_t i = 0; i < options.size(); i++)
    if(options[i].required &&
       (option_check[i] == 0 || (options[i].type == ENABLE && option_check[++i] == 0))) {
      error_type = 6; // required option missing
      error_idx = i; // for a boolean option i is the index of DISABLE
      return;
    }

  error_type = 0;
}

#if defined(HAVE_PETSC)
int feOptionsParser::ParsePetsc()
{
  // Actually parse the command line options
  Parse();

  // Continue parsing with PETSc if the parse went well,
  // that is, if no error or if an unknown parameter was met
  // (this is OK as we want to include all PETSc options
  // which are not listed).
  // In particular, this allows to exit if a recognized option
  // was given an incompatible type.
  if(error_type != 0 && error_type != 2)
    return 1;

  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "", "");
  {
    // Always add the help to the PETSc database
    char dummyHelp[64];
    PetscCall(PetscOptionsString("-h", "--help", "Display the available options", dummyHelp, dummyHelp, sizeof(dummyHelp), NULL));

    // Store each option in a dummy variable, but signal to PETSc that the option exists
    // and has been populated (in the Parse() call above)
    for(auto &option : options) {
      switch(option.type)
      {
        case INT:
        {
          int foo = 0;
          PetscCall(PetscOptionsInt(option.short_name, option.long_name, option.description, foo, &foo, NULL));
          break;
        }
        case DOUBLE:
        {
          double foo = 0.;
          PetscCall(PetscOptionsReal(option.short_name, option.long_name, option.description, foo, &foo, NULL));
          break;
        }
        case STRING:
        {
          char foo[64];
          PetscCall(PetscOptionsString(option.short_name, option.long_name, option.description, foo, foo, sizeof(foo), NULL));
          break;
        }
        case ENABLE:
        {
          PetscBool foo = PETSC_TRUE;
          PetscBool flg;
          PetscCall(PetscOptionsBool(option.short_name, option.long_name, option.description, foo, &flg, NULL));
          break;
        }
        case DISABLE:
        {
          PetscBool foo = PETSC_FALSE;
          PetscBool flg;
          PetscCall(PetscOptionsBool(option.short_name, option.long_name, option.description, foo, &flg, NULL));
          break;
        }
      }
    }
  }
  PetscOptionsEnd();
  return 0;
}
#endif

feStatus feOptionsParser::parse()
{
#if defined(HAVE_PETSC)
  ParsePetsc();
  // error_type == 2 is OK to accept unlisted PETSc options (e.g. -ksp_view)
  if(error_type == 0 || error_type == 2) {
    return FE_STATUS_OK;
  } else if(error_type == 1) {
    PrintUsage(std::cout);
    return FE_PRINT_HELP;
  } else {
    PrintUsage(std::cout);
    return feErrorMsg(FE_STATUS_ERROR, "Parse error");
  }
#else
  Parse();
  if(error_type == 0) {
    return FE_STATUS_OK;
  } else if(error_type == 1) {
    PrintUsage(std::cout);
    return FE_PRINT_HELP;
  } else {
    PrintUsage(std::cout);
    return feErrorMsg(FE_STATUS_ERROR, "Parse error");
  }
#endif
}

void feOptionsParser::WriteValue(const Option &opt, std::ostream &out)
{
  switch(opt.type) {
    case INT:
      out << *(int *)(opt.var_ptr);
      break;

    case DOUBLE:
      out << *(double *)(opt.var_ptr);
      break;

    case STRING:
      out << *(const char **)(opt.var_ptr);
      break;

    default: // provide a default to suppress warning
      break;
  }
}

void feOptionsParser::PrintError(std::ostream &out) const
{
  static const char *line_sep = "";

  out << line_sep;
  switch(error_type) {
    case 2:
      out << "Unrecognized option: " << argv[error_idx] << '\n' << line_sep;
      break;

    case 3:
      out << "Missing argument for the last option: " << argv[argc - 1] << '\n' << line_sep;
      break;

    case 4:
      if(options[error_idx].type == ENABLE)
        out << "Option " << options[error_idx].long_name << " or "
            << options[error_idx + 1].long_name << " provided multiple times\n"
            << line_sep;
      else if(options[error_idx].type == DISABLE)
        out << "Option " << options[error_idx - 1].long_name << " or "
            << options[error_idx].long_name << " provided multiple times\n"
            << line_sep;
      else
        out << "Option " << options[error_idx].long_name << " provided multiple times\n"
            << line_sep;
      break;

    case 5:
      out << "Wrong option format: " << argv[error_idx - 1] << " " << argv[error_idx] << '\n'
          << line_sep;
      break;

    case 6:
      out << "Missing required option: " << options[error_idx].long_name << '\n' << line_sep;
      break;
  }
  out << std::endl;
}

void feOptionsParser::PrintHelp(std::ostream &out) const
{
  static const char *indent = "   ";
  static const char *seprtr = ", ";
  static const char *descr_sep = "\n\t";
  static const char *line_sep = "";
  static const char *types[] = {" <int>", " <double>",   " <string>",     "",
                                "",       " '<int>...'", " '<double>...'"};

  out << indent << "-h" << seprtr << "--help" << descr_sep << "Print this help message and exit.\n"
      << line_sep;
  for(size_t j = 0; j < options.size(); j++) {
    OptionType type = options[j].type;

    out << indent << options[j].short_name << types[type] << seprtr << options[j].long_name
        << types[type] << seprtr;
    if(options[j].required) {
      out << "(required)";
    } else {
      if(type == ENABLE) {
        j++;
        out << options[j].short_name << types[type] << seprtr << options[j].long_name << types[type]
            << seprtr << "current option: ";
        if(*(bool *)(options[j].var_ptr) == true) {
          out << options[j - 1].long_name;
        } else {
          out << options[j].long_name;
        }
      } else {
        out << "current value: ";
        WriteValue(options[j], out);
      }
    }
    out << descr_sep;

    if(options[j].description) {
      out << options[j].description << '\n';
    }
    out << line_sep;
  }
}

void feOptionsParser::PrintUsage(std::ostream &out) const
{
  static const char *line_sep = "";

  PrintError(out);
  out << "Usage: " << argv[0] << " [options] ...\n" << line_sep << "Options:\n" << line_sep;
  PrintHelp(out);
}
/*
  Copied and adapted from the Hxt project
  See https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/contrib/hxt/core/src
  Original author : Celestin Marot
*/

#ifndef _FEMESSAGE_
#define _FEMESSAGE_

/* Verbosity level :
  - 0 : No information messages, only print warnings and errors
  - 1 (default) : Moderate information messages
  - 2 : All information messages */
#define VERBOSE_NONE 0
#define VERBOSE_MODERATE 1
#define VERBOSE_ALL 2

typedef enum {
  FE_MSGLEVEL_INFO = 0,
  FE_MSGLEVEL_INFO_COLLECTIVE = 5,
  FE_MSGLEVEL_DEBUG = 1,
  FE_MSGLEVEL_WARNING = 2,
  FE_MSGLEVEL_ERROR = 3,
  FE_MSGLEVEL_TRACE = 4
} feLevel;

/* MESSAGE */
typedef struct {
  /* the message content */
  const char *string; // lifetime = time of callback function

  /* information about the location of the code which sent the message */
  const char *func; // lifetime = forever
  const char *file; // lifetime = forever
  const char *line; // lifetime = forever
  int threadId; // the thread which sent the message
  int numThreads; // the number of threads

  feLevel level;

} feMessage;

typedef enum {
  // Positive values mean a success => feCheck does nothing for positive values
  FE_STATUS_OK = 0,
  FE_STATUS_TRUE = 0,
  FE_STATUS_FALSE = 1,
  FE_PRINT_HELP = 2,
  // Fatal Errors => feCheck gives trace message and return
  FE_STATUS_ERROR = -1,
  FE_STATUS_FAILED = -2,
  FE_STATUS_ASSERTION_FAILED = -3,
  FE_STATUS_OUT_OF_MEMORY = -4,
  FE_STATUS_FILE_CANNOT_BE_OPENED = -5,
  FE_STATUS_POINTER_ERROR = -6,
  FE_STATUS_READ_ERROR = -7,
  FE_STATUS_WRITE_ERROR = -8,
  FE_STATUS_RANGE_ERROR = -9,
  FE_STATUS_FORMAT_ERROR = -10,

} feStatus;

#define STR(x) #x
#define STRINGIFY(x) STR(x)

#define feInfo(...) feMessageInfo(__func__, __FILE__, STRINGIFY(__LINE__), ##__VA_ARGS__)
#define feInfoCond(cond, ...) ((cond) ? feInfo(__VA_ARGS__) : 0)
#define feInfoCollective(...) feMessageInfoCollective(__func__, __FILE__, STRINGIFY(__LINE__), ##__VA_ARGS__)
#define feWarning(...) feMessageWarning(__func__, __FILE__, STRINGIFY(__LINE__), ##__VA_ARGS__)

/* print an error message corresponding to a status */
#define feErrorMsg(status, ...)                                                                    \
  feMessageError(status, __func__, __FILE__, STRINGIFY(__LINE__), ##__VA_ARGS__)
#define feError(status) feErrorMsg(status, NULL)

/* give trace message if status is not OK, but does not return */
#define FE_TRACE_MSG(status, ...)                                                                  \
  feMessageTraceError(status, __func__, __FILE__, STRINGIFY(__LINE__), ##__VA_ARGS__)
#define FE_TRACE(status) FE_TRACE_MSG(status, NULL)

/* give trace message and return the error status if not ok */
#define FE_CHECK_MSG(status, ...)                                                                  \
  do {                                                                                             \
    feStatus _tmp_ = status;                                                                       \
    if(_tmp_ < 0) {                                                                                \
      if(_tmp_ < 0) {                                                                              \
        FE_TRACE_MSG(_tmp_, ##__VA_ARGS__);                                                        \
      }                                                                                            \
      return _tmp_;                                                                                \
    }                                                                                              \
    if(_tmp_ == 2) {                                                                               \
      return _tmp_;                                                                                \
    }                                                                                              \
  } while(0)
#define feCheck(status) FE_CHECK_MSG(status, NULL)
#define feCheckReturn(status) if(status != FE_STATUS_OK) return status;

void setVerbose(int level);

// Prints info on MPI rank 0 if MPI is enabled
feStatus feMessageInfo(const char *func, const char *file, const char *line, const char *fmt, ...);
// Prints info on all MPI ranks
feStatus feMessageInfoCollective(const char *func, const char *file, const char *line, const char *fmt, ...);

feStatus feMessageWarning(const char *func, const char *file, const char *line, const char *fmt,
                          ...);
feStatus feMessageError(feStatus status, const char *func, const char *file, const char *line,
                        const char *fmt, ...);
feStatus feMessageTraceError(feStatus status, const char *func, const char *file, const char *line,
                             const char *fmt, ...);

#endif
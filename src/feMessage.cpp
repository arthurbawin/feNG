/*
  Copied and adapted from the Hxt project
  See https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/contrib/hxt/core/src
  Original author : Celestin Marot
*/

#include "feNG.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#if defined(HAVE_OMP)
#include "omp.h"
#endif
#if defined(HAVE_MPI)
#include "mpi.h"
#endif
#include "feMessage.h"

int FE_VERBOSE = VERBOSE_MODERATE;

void setVerbose(int level) { FE_VERBOSE = level; }

feStatus defaultMessageCallback(feMessage *msg);

static feStatus (*msgCallback)(feMessage *msg) =
  &defaultMessageCallback; // the message callback is declared as a global variable
static const char *encodingIssue = "~~~~ (encoding issue) ~~~~";

const char *feGetStatusString(feStatus status)
{
  switch(status) {
    case FE_STATUS_OK:
      return "no error";
      break;
    case FE_STATUS_FAILED:
      return "function failed";
      break;
    case FE_STATUS_OUT_OF_MEMORY:
      return "out of memory";
      break;
    case FE_STATUS_ERROR:
      return "error";
      break;
    case FE_STATUS_FILE_CANNOT_BE_OPENED:
      return "file cannot be opened";
      break;
    case FE_STATUS_ASSERTION_FAILED:
      return "assertion failed";
      break;
    case FE_STATUS_POINTER_ERROR:
      return "wrong pointer given";
      break;
    case FE_STATUS_READ_ERROR:
      return "read error";
      break;
    case FE_STATUS_WRITE_ERROR:
      return "write error";
      break;
    case FE_STATUS_RANGE_ERROR:
      return "number out of range";
      break;
    case FE_STATUS_FORMAT_ERROR:
      return "wrong format";
      break;
    default:
      if(status < 0)
        return "unknown error";
      else
        return "positive return value (no error)";
      break;
  }
}

feStatus defaultMessageCallback(feMessage *msg)
{
  if(msg->level == FE_MSGLEVEL_INFO)
    fprintf(stdout, "Info : %s\n", msg->string);
  else if(msg->level == FE_MSGLEVEL_INFO_COLLECTIVE) {
#if defined(HAVE_MPI)
      int rank, size;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      char processor_name[MPI_MAX_PROCESSOR_NAME];
      int name_len;
      MPI_Get_processor_name(processor_name, &name_len);
      fprintf(stdout, "Info %s - %d/%d : %s\n", processor_name, rank, size, msg->string);
#endif
  }
  else if(msg->level == FE_MSGLEVEL_ERROR)
    fprintf(stderr, "\n= X = Error : %s   \n in %s -> %s:%s\n", msg->string, msg->func, msg->file,
            msg->line);
  else if(msg->level == FE_MSGLEVEL_TRACE)
    fprintf(stderr, "  - trace -   %s -> %s:%s  \t(%s)\n", msg->func, msg->file, msg->line,
            msg->string);
  else if(msg->level == FE_MSGLEVEL_WARNING)
    fprintf(stderr, "/!\\ Warning : %s\n", msg->string);
  else if(msg->level == FE_MSGLEVEL_DEBUG)
    fprintf(stderr, "Debug : %s   \t(in %s -> %s:%s)\n", msg->string, msg->func, msg->file,
            msg->line);
  return FE_STATUS_OK;
}

feStatus feSetMessageCallback(feStatus (*newMsgCallback)(feMessage *msg))
{
  if(newMsgCallback == NULL)
    msgCallback = defaultMessageCallback;
  else
    msgCallback = newMsgCallback;
  return FE_STATUS_OK;
}

static void feMessageGeneral(int messageLevel, const char *func, const char *file, const char *line,
                             const char *encodingIssueString, const char *alternativeString,
                             const char *fmt, va_list args)
{
  feMessage msg;

  char str[4096];
  if(fmt != NULL) {
    int err = vsnprintf(str, sizeof(str), fmt, args);

    if(err >= 0)
      msg.string = str;
    else
      msg.string = encodingIssueString;
  } else {
    msg.string = alternativeString;
  }

  msg.func = func;
  msg.file = file;
  msg.line = line;
  msg.level = (feLevel)messageLevel;

#if defined(HAVE_OMP)
  msg.threadId = omp_get_thread_num();
  msg.numThreads = omp_get_num_threads();
#endif
  // #pragma omp critical
  {
    msgCallback(&msg);
  }
}

feStatus feMessageInfo(const char *func, const char *file, const char *line, const char *fmt, ...)
{
#if defined(HAVE_MPI)
  int initialized;
  MPI_Initialized(&initialized);
  if(initialized) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
#endif
    {
      va_list args;
      va_start(args, fmt);
      feMessageGeneral(FE_MSGLEVEL_INFO, func, file, line, encodingIssue, "", fmt, args);
      va_end(args);
    }
#if defined(HAVE_MPI)
  }
#endif
  return FE_STATUS_OK;
}

feStatus feMessageInfoCollective(const char *func, const char *file, const char *line, const char *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  feMessageGeneral(FE_MSGLEVEL_INFO_COLLECTIVE, func, file, line, encodingIssue, "", fmt, args);
  va_end(args);
  return FE_STATUS_OK;
}

feStatus feMessageWarning(const char *func, const char *file, const char *line, const char *fmt,
                          ...)
{
  va_list args;
  va_start(args, fmt);
  feMessageGeneral(FE_MSGLEVEL_WARNING, func, file, line, encodingIssue, "", fmt, args);
  return FE_STATUS_OK;
}

feStatus feMessageError(feStatus status, const char *func, const char *file, const char *line,
                        const char *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  feMessageGeneral(FE_MSGLEVEL_ERROR, func, file, line, feGetStatusString(status),
                   feGetStatusString(status), fmt, args);
  va_end(args);
  return status;
}

feStatus feMessageTraceError(feStatus status, const char *func, const char *file, const char *line,
                             const char *fmt, ...)
{
  if(status == FE_STATUS_OK) return FE_STATUS_OK;
  va_list args;
  va_start(args, fmt);
  feMessageGeneral(FE_MSGLEVEL_TRACE, func, file, line, encodingIssue, feGetStatusString(status),
                   fmt, args);
  va_end(args);
  return status;
}

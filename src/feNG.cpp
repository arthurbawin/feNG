#include "feNG.h"

using namespace std::chrono_literals;

double tic(int mode)
{
  static std::chrono::steady_clock::time_point begin;

  if(mode == 0){
    begin = std::chrono::steady_clock::now();
    return 0.;
  }
  else {
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> d = end - begin;
    double timeInSeconds  = d / 1.0s;
    return timeInSeconds;
  }
}

double toc() { return tic(1); }
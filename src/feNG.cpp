#include "feNG.h"

void tic(int mode)
{
  static std::chrono::_V2::system_clock::time_point t_start;

  if(mode == 0)
    t_start = std::chrono::high_resolution_clock::now();
  else {
    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout << "Elapsed time is " << (t_end - t_start).count() * 1E-9 << " seconds\n";
  }
}

void toc() { tic(1); }
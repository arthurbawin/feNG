#ifndef _FENG_TESTS_
#define _FENG_TESTS_

int compareOutputFiles(const std::string &testRoot, const std::stringstream &resultBuffer)
{
  std::string expectedFile = testRoot + ".output";

	// Read result and expected output and compare both buffers
  std::ifstream refFile(expectedFile);
  std::stringstream expectedBuffer;
  expectedBuffer << refFile.rdbuf();
  std::cout << std::endl;
  std::cout << "===============================================" << std::endl;
  std::cout << "Expected output is:" << std::endl;
  std::cout << "===============================================" << std::endl;

  if(refFile) {
    std::cout << expectedBuffer.str() << std::endl;
  } else {
    std::string tmpOutputFile = testRoot + ".tmp_output";

    std::cout << std::endl;
    feWarning("Output file %s was not found! Test is likely to fail.", expectedFile.data());
    feWarning("Writing temporary output file : %s", tmpOutputFile.data());
    feWarning("If this file matches the expected output, rename it as %s", expectedFile.data());
    feWarning("Attention! This will set this result as the new reference test output.");
    std::cout << std::endl;

    std::ofstream tmpOutput;
    tmpOutput.open(tmpOutputFile);
    if(!tmpOutput) {
      std::cout << std::endl;
      feWarning("!!! Could not open temporary output file : %s !!!", tmpOutputFile.data());
      std::cout << std::endl;
    } else {
      tmpOutput << resultBuffer.rdbuf();
      tmpOutput.close();
    }
  }
  std::cout << "===============================================" << std::endl;
  std::cout << "End of expected output" << std::endl;
  std::cout << "===============================================" << std::endl;

  std::cout << std::endl;
  std::cout << "===============================================" << std::endl;
  std::cout << "Test output is:" << std::endl;
  std::cout << "===============================================" << std::endl;
  std::cout << resultBuffer.str() << std::endl;
  std::cout << "===============================================" << std::endl;
  std::cout << "End of test output" << std::endl;
  std::cout << "===============================================" << std::endl;

  return expectedBuffer.str().compare(resultBuffer.str());
}

void computeAndPrintConvergence(const int dim,
                                const int nMesh,
                                const std::vector<double> &error,
                                const std::vector<int> &nElm,
                                std::stringstream &resultBuffer)
{
  std::vector<double> rate(nMesh, 0.);
  for(int i = 1; i < nMesh; ++i)
  {
    rate[i] = -log(error[i] / error[i-1]) / log( pow((double) nElm[i] / (double) nElm[i-1], 1./(double) dim) );
  }
  printf("%12s \t %12s \t %12s \n", "nElm", "||E||", "rate");
  resultBuffer
    << std::setw(16) << std::right << "nElm"
    << std::setw(16) << std::right << "error"
    << std::setw(16) << std::right << "rate" << std::endl;
  for(int i = 0; i < nMesh; ++i) {
    printf("%12d \t %12.6e \t %12.6e\n", nElm[i], error[i], rate[i]);
    resultBuffer
      << std::scientific
      << std::setw(16) << std::right
      << std::setprecision(6) << nElm[i] << std::setw(16) << std::right
      << std::setprecision(6) << error[i] << std::setw(16) << std::right
      << std::setprecision(6) << rate[i]
      << std::endl;
  }
  resultBuffer << std::endl;
}

#endif
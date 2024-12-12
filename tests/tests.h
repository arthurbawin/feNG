#ifndef _FENG_TESTS_
#define _FENG_TESTS_

extern int my_argc;
extern char** my_argv;

int compareOutputFiles(std::string &testRoot, std::stringstream &resultBuffer)
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

#endif
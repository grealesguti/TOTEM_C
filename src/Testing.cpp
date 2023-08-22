#include "Testing.h"

using namespace arma;

Testing::Testing() {
    // Initialize a member variable (if needed)
}


int Testing::testArmadillo() {
  arma::mat A = {{1, 2}, {3, 4}};
  arma::mat B = {{5, 6}, {7, 8}};
  arma::mat result = A + B;

  // Assert the expected result
  arma::mat expected = {{6, 8}, {10, 12}};
  if (arma::approx_equal(result, expected, "absdiff", 1e-6)) {
    std::cout << "Armadillo integration test passed!" << std::endl;
    return 0;
  } else {
    std::cout << "Armadillo integration test failed!" << std::endl;
  }
}
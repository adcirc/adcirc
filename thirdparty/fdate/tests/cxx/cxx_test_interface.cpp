//
// This creates a test interface for the C++ code to check that
// date parsing and operations match a python implementation.
// Then, the python can be called with random data and we can
// see if our implementation matches
//

#include <iostream>
#include <string>

#include "DateTime.hpp"
#include "TimeDelta.hpp"

int parse_test(int argc, char *argv[]) {
  // Check for sanity
  if (argc < 3) {
    std::cerr << "ERROR: Invalid command line arguments" << std::endl;
    exit(1);
  }

  // Convert the argv to a string
  const auto input_string = std::string(argv[2]);

  // Create a DateTime object
  const auto dt = DateTime::strptime(input_string);

  // Check if the DateTime object is valid
  if (!dt.valid()) {
    std::cerr << "Invalid DateTime string: " << input_string << std::endl;
    return 1;
  }

  // Write out the number of seconds since the epoch
  std::cout << dt.timestamp() << std::endl;

  return 0;
}

int arithmetic_test(int argc, char *argv[]) {
  // Check for sanity
  if (argc < 7) {
    std::cerr << "ERROR: Invalid command line arguments" << std::endl;
    exit(1);
  }

  // Convert the argv to a string
  const auto input_string = std::string(argv[2]);

  // Create a DateTime object
  const auto start_time = DateTime::strptime(input_string);

  // Check if the DateTime object is valid
  if (!start_time.valid()) {
    std::cerr << "Invalid DateTime string: " << input_string << std::endl;
    return 1;
  }

  // Get the time delta
  const int days = std::stoi(argv[3]);
  const int hours = std::stoi(argv[4]);
  const int minutes = std::stoi(argv[5]);
  const int seconds = std::stoi(argv[6]);
  const auto dt = TimeDelta(days, hours, minutes, seconds, 0);

  // Add
  const auto result = start_time + dt;

  // Write out the result
  std::cout << result.strftime("%Y-%m-%dT%H:%M:%S") << std::endl;

  return 0;
}

int main(int argc, char *argv[]) {
  // Check that we have the right number of arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0]
              << " <test_type> <datetime_string> OR <test_type> "
                 "<datetime_string> <days> <hours> <minutes> <seconds>"
              << std::endl;
    return 1;
  }

  const auto type = std::string(argv[1]);
  if (type == "parse") {
    if (argc != 3) {
      std::cerr << "Invalid input" << std::endl;
      return 1;
    }
    return parse_test(argc, argv);
  } else if (type == "arithmetic") {
    if (argc != 7) {
      std::cerr << "Invalid input" << std::endl;
      return 1;
    }
    return arithmetic_test(argc, argv);
  }
}

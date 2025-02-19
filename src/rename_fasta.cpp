#include <Rcpp.h>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

using namespace Rcpp;

// [[Rcpp::export]]
void rename_fasta(std::string input_file, std::string output_file, std::string prefix) {
  std::ifstream infile(input_file);
  std::ofstream outfile(output_file);
  std::string line;
  int count = 1;

  if (!infile.is_open()) {
    stop("Could not open input file.");
  }
  if (!outfile.is_open()) {
    stop("Could not open output file.");
  }

  while (std::getline(infile, line)) {
    if (line[0] == '>') {
      std::ostringstream new_header;
      new_header << ">" << prefix << "_" << std::setw(9) << std::setfill('0') << count;
      outfile << new_header.str() << std::endl;
      count++;
    } else {
      outfile << line << std::endl;
    }
  }

  infile.close();
  outfile.close();
}
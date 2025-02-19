#include <Rcpp.h>
#include <fstream>
#include <string>
#include <unordered_map>
#include <sstream>
#include <cstdio> // for std::remove and std::rename

using namespace Rcpp;

// [[Rcpp::export]]
void seq2tab(std::string input_table_file, std::string fasta_file, std::string id_column) {
  // Create a map to store sequences from the FASTA file
  std::unordered_map<std::string, std::string> sequences;

  // Read the FASTA file
  std::ifstream infile(fasta_file);
  if (!infile.is_open()) {
    stop("Could not open FASTA file.");
  }

  std::string line, current_id, current_sequence;
  while (std::getline(infile, line)) {
    if (line[0] == '>') {
      if (!current_id.empty()) {
        sequences[current_id] = current_sequence;
      }
      current_id = line.substr(1); // Remove '>'
      current_sequence.clear();
    } else {
      current_sequence += line;
    }
  }
  if (!current_id.empty()) {
    sequences[current_id] = current_sequence;
  }
  infile.close();

  // Read the input table and write the output table with the new sequence column
  std::ifstream table_infile(input_table_file);
  std::string temp_table_file = input_table_file + ".tmp";
  std::ofstream table_outfile(temp_table_file);
  if (!table_infile.is_open()) {
    stop("Could not open input table file.");
  }
  if (!table_outfile.is_open()) {
    stop("Could not open temporary output table file.");
  }

  std::string header;
  std::getline(table_infile, header);
  table_outfile << header << "\tsequence\n"; // Use tab as the delimiter

  std::string row;
  while (std::getline(table_infile, row)) {
    std::istringstream row_stream(row);
    std::string cell, id;
    std::getline(row_stream, id, '\t'); // Assuming the ID column is the first column and tab-separated

    std::string sequence = (sequences.find(id) != sequences.end()) ? sequences[id] : "NA";
    table_outfile << row << "\t" << sequence << "\n"; // Use tab as the delimiter
  }

  table_infile.close();
  table_outfile.close();

  // Replace the original file with the temporary file
  std::remove(input_table_file.c_str());
  std::rename(temp_table_file.c_str(), input_table_file.c_str());
}
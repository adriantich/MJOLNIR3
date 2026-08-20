#include "seq2tab_core.h"

#include <fstream>
#include <string>
#include <unordered_map>
#include <sstream>
#include <cstdio> // for std::remove and std::rename

/*
In order to be able to create a standalone binary, the seq2tab.cpp in now 
splitted in different files. There is a core of the function and then a 
wrapper for R and another for the CLI.
The core is in src/ folder: seq2tab_core.h and seq2tab_core.cpp. 
The wrapper for R is in src/ too: seq2tab.cpp
The wrapper for the CLI is in tools/: seq2tab_cli.cpp
*/ 

void seq2tab_core(const std::string& input_table_file,
                  const std::string& fasta_file,
                  const std::string& id_column) {
  // Create a map to store sequences from the FASTA file
  std::unordered_map<std::string, std::string> sequences;

  // Read the FASTA file
  std::ifstream infile(fasta_file);
  if (!infile.is_open()) {
    // stop("Could not open FASTA file.");
    throw std::runtime_error("Could not open FASTA file.");
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
  if (!table_infile.is_open()) {
    // stop("Could not open input table file.");
    throw std::runtime_error("Could not open input table file.");
  }

  std::string temp_table_file = input_table_file + ".tmp";
  std::ofstream table_outfile(temp_table_file);
  if (!table_outfile.is_open()) {
    // stop("Could not open temporary output table file.");
    throw std::runtime_error("Could not open temporary output table file.");
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
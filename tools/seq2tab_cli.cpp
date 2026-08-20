#include "../src/seq2tab_core.h"

#include <exception>
#include <iostream>
#include <string>

/*
In order to be able to create a standalone binary, the seq2tab.cpp in now 
splitted in different files. There is a core of the function and then a 
wrapper for R and another for the CLI.
The core is in src/ folder: seq2tab_core.h and seq2tab_core.cpp. 
The wrapper for R is in src/ too: seq2tab.cpp
The wrapper for the CLI is in tools/: seq2tab_cli.cpp
*/ 

int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cerr << "Usage: seq2tab <input_table> <fasta_file> <id_column>\n";
    return 1;
  }

  try {
    seq2tab_core(argv[1], argv[2], argv[3]);
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }
}
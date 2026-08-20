#pragma once

#include <string>

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
                  const std::string& id_column);
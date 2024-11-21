#ifndef amino_acid_sequence_reader
#define amino_acid_sequence_reader

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>


std::string read_fasta(std::string& filepath);


class input_parameters{
    public:
        input_parameters(std::string filepath);
        std::string job;
        std::vector<std::vector<int>> box_size;
        int steps;
        double starting_temperature;
        bool anneal, anneal_to_zero;
        double rate;
        double max_distance;
        std::string sequence;


};

#endif // amino_acid_sequence_reader

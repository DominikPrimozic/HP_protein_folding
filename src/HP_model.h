#ifndef polar_nonpolar_model
#define polar_nonpolar_model

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <utility>
#include <functional>
#include <sequence_reader.h>
#include <set>



using coordinates = std::vector<std::vector<int>>;
//hashing cooridnates into unordered map is faster for lookup

class HP_model_2D : public input_parameters{ //deriving class from input
    private:
        void initial_();
        double conformatio_energy(coordinates& new_protein);
        bool Metropolis_algorithm(double newE, double oldE, double T=300);
        coordinates move(coordinates protein, int& movable);
        void move_residue(coordinates& protein, int& movable);
        double compactness(coordinates& new_protein);
        double annealing(int& step, bool& isTrue,double rate,  bool toZero);


        std::vector<std::vector<int>> lattice; //lattice only has [minx,miny],[maxx,maxy]
        coordinates protein;
        std::random_device rd;
        std::mt19937 gen;
        std::uniform_int_distribution<> distrib;
        int aa_size;
        


    public:
        HP_model_2D(std::string filepath);
        void print_to_file(std::string path,double E = std::nan(""),double C = std::nan(""),double T = std::nan(""));
        void run_simulation();
        double max_compactness, min_energy;
        int best_compactness, best_energy;
};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class HP_model_3D : public input_parameters{ //deriving class from input
    private:
        void initial_();
        double conformatio_energy(coordinates& new_protein);
        bool Metropolis_algorithm(double newE, double oldE, double T=300);
        coordinates move(coordinates protein, int& movable);
        void move_residue(coordinates& protein, int& movable);
        double compactness(coordinates& new_protein);
        double annealing(int& step, bool& isTrue,double rate,  bool toZero);

        bool is_valid_move_residue(coordinates& new_protein,std::vector<std::vector<int>>& lattice, std::vector<int>& new_position, int& residue, double distance);
        int neighbours(coordinates& protein, int& current,std::string sequence );
        bool is_valid_move(coordinates& new_protein,std::vector<std::vector<int>>& lattice);

        std::vector<std::vector<int>> lattice; //lattice only has [minx,miny, minz],[maxx,maxy, maxz]
        coordinates protein;
        std::random_device rd;
        std::mt19937 gen;
        std::uniform_int_distribution<> distrib;
        int aa_size;



    public:
        HP_model_3D(std::string filepath);
        void print_to_file(std::string path,double E = std::nan(""),double C = std::nan(""),double T = std::nan(""));
        void run_simulation();
        double max_compactness, min_energy;
        int best_compactness, best_energy;
};
#endif // polar_nonpolar_model

#include <HP_model.h>
#include <sequence_reader.h>
#include <filesystem>

namespace fs = std::filesystem;
void cleanTxtFiles(const std::string& directoryPath) { //I used chatGPT for this because I have no idea how to program this
    try {
        // Check if the directory exists
        if (fs::exists(directoryPath) && fs::is_directory(directoryPath)) {
            for (const auto& entry : fs::directory_iterator(directoryPath)) {
                // Check if the entry is a file and has a .txt extension
                if (entry.is_regular_file() && entry.path().extension() == ".txt") {
                    std::cout << "Removing file: " << entry.path() << '\n';
                    fs::remove(entry.path()); // Remove the file
                }
            }
        } else {
            std::cerr << "The directory does not exist or is not a directory.\n";
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << '\n';
    } catch (const std::exception& e) {
        std::cerr << "General error: " << e.what() << '\n';
    }
}

std::string trim_main(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return ""; 
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, (last - first + 1));
}

int main(){

    std::string directory = "output/conformations";
    cleanTxtFiles(directory);
    
    std::string filepath = "input/input.txt";

    std::ifstream inputFile(filepath);
    if (!inputFile) {
        std::cerr << "Error: Could not open input file.\n";
        return 1;
    }
    std::string line;
    int mode;
    if (std::getline(inputFile, line)) {
        line = trim_main(line);
        
        if (line == "3D") {
            mode=1;
        } else if (line == "2D") {
            mode=0;
        } else {
            std::cout << "Unknown dimensionality specified in the file.\n";
        }
    } else {
        std::cerr << "Error: Could not read the first line of the file.\n";
    }

    inputFile.close();

    if (mode==1){
        HP_model_3D hp(filepath);
        hp.run_simulation();
        std::ofstream outfile("output/best.txt");
        outfile << "Minimal energy " << hp.min_energy << " at configuration: " << hp.best_energy << std::endl;
        outfile << "Maximal compactness " << hp.max_compactness << " at configuration: " << hp.best_compactness << std::endl;
        outfile.close();
    
    }
    else if (mode==0){
        HP_model_2D hp(filepath);
        hp.run_simulation();
        std::ofstream outfile("output/best.txt");
        outfile << "Minimal energy " << hp.min_energy << " at configuration: " << hp.best_energy << std::endl;
        outfile << "Maximal compactness " << hp.max_compactness << " at configuration: " << hp.best_compactness << std::endl;
        outfile.close();
    }
    else{
        std::cerr << "Error: Could not run.\n";
    }
    
    return 0;
}


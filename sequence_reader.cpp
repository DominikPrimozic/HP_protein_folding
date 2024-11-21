#include <sequence_reader.h>

std::string read_fasta(std::string& filepath){
    std::ifstream fastaFile(filepath);

    std::string line;
    std::string sequence;

    while (std::getline(fastaFile, line)) {
        if (line.empty()) continue;
        else if (line[0] == '>') continue;
        else{
            sequence +=line;
        }
    }
    fastaFile.close();
    return sequence;
}





std::string trim(const std::string& s) {
    auto start = s.begin();
    while (start != s.end() && std::isspace(*start)) start++;
    
    auto end = s.end();
    do {
        end--;
    } while (std::distance(start, end) > 0 && std::isspace(*end));
    
    return std::string(start, end + 1);
}

// Remove comments and trim whitespace
std::string removeComment(const std::string& line) {
    size_t pos = line.find("//");
    std::string cleanLine = (pos != std::string::npos) ? line.substr(0, pos) : line;
    return trim(cleanLine);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

input_parameters::input_parameters(std::string filepath){
    std::ifstream file(filepath);
    std::string line;

    std::getline(file, job);

    while (std::getline(file, line)){

        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);

        if (line.rfind("//", 0) == 0) {
        line = line.substr(2);
            }
        
        if (line.empty() || line[0] == '#') continue;

        else if (line.find("box size") != std::string::npos) {
            box_size.clear();  

            while (std::getline(file, line) && !line.empty()) {  
                std::istringstream boxStream(line);
                std::vector<int> dimensions;
                int value;

                
                while (boxStream >> value) {
                    dimensions.push_back(value);
                    if (boxStream.peek() == ',' || boxStream.peek() == ' ') {
                        boxStream.ignore();  
                    }
                }

                
                if (!dimensions.empty()) {
                    box_size.push_back(dimensions);
                }

                
                if (box_size.size() == 2) {
                    break;
                }
            }
        }
        else if (line.find("steps") != std::string::npos) {
            std::getline(file, line);
            steps = std::stoi(line);
        }
        else if (line.find("starting temperature") != std::string::npos) {
            std::getline(file, line);
            starting_temperature = std::stod(line);
        }
        else if (line.find("anneal") != std::string::npos) {
            std::getline(file, line);
            line=removeComment(line);
            anneal = (line == "true");

            std::getline(file, line);
            line=removeComment(line);
            anneal_to_zero = (line == "true");

            if (std::getline(file, line)) {
                
                if (!line.empty() && line.find("//") != std::string::npos) {
                    rate = std::stod(line);  
                }
            }
        }
        else if (line.find("max distance") != std::string::npos) {
            std::getline(file, line);
            max_distance = std::stod(line);
        }
        
        else if (line.find("sequence") != std::string::npos) {
            std::string seq_line;
            while (std::getline(file, seq_line) && !seq_line.empty()) {
                if (seq_line[0] == '>') continue;
                else {
                    sequence += seq_line;
            }
        }
    }
    }

    file.close();
        
}
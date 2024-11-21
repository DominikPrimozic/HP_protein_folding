#include <HP_model.h>


void printVectorOfVectors(const std::vector<std::vector<int>>& vec) {
     std::cout <<"vector is "<< std::endl;
    for (const auto& pair : vec) {
        std::cout << "(" << pair[0] << ", " << pair[1] << ") ";
    }
    std::cout << std::endl;
}

int HP_translator(char& aa){ 
    std::unordered_map<char, int> amino_acid_map = {
        {'A', 0}, // Alanine - Hydrophobic
        {'R', 1}, // Arginine - Polar
        {'N', 1}, // Asparagine - Polar
        {'D', 1}, // Aspartic acid - Polar
        {'C', 0}, // Cysteine - Hydrophobic
        {'Q', 1}, // Glutamine - Polar
        {'E', 1}, // Glutamic acid - Polar
        {'G', 0}, // Glycine - Hydrophobic
        {'H', 1}, // Histidine - Polar
        {'I', 0}, // Isoleucine - Hydrophobic
        {'L', 0}, // Leucine - Hydrophobic
        {'K', 1}, // Lysine - Polar
        {'M', 0}, // Methionine - Hydrophobic
        {'F', 0}, // Phenylalanine - Hydrophobic
        {'P', 0}, // Proline - Hydrophobic
        {'S', 1}, // Serine - Polar
        {'T', 1}, // Threonine - Polar
        {'W', 0}, // Tryptophan - Hydrophobic
        {'Y', 1}, // Tyrosine - Polar
        {'V', 0}  // Valine - Hydrophobic
    };
    return amino_acid_map[aa];
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


HP_model_2D::HP_model_2D(std::string  filepath) : input_parameters(filepath), gen(rd()), distrib(0,3){
    aa_size=sequence.size();
    initial_();
}

bool is_valid_move(coordinates& new_protein,std::vector<std::vector<int>>& lattice){
    
    int min_x = lattice[0][0];
    int min_y = lattice[0][1];
    int max_x = lattice[1][0];
    int max_y = lattice[1][1];

    std::set<std::vector<int>> collected_residues;
    for (auto& residue : new_protein){
        if (residue[0] < min_x || residue[0] > max_x || residue[1] < min_y || residue[1]> max_y) {
            return false; 
        }

        if (collected_residues.find(residue) != collected_residues.end()) {
            return false;
        }
        collected_residues.insert(residue);
        
    }
    return true;
}

coordinates up_fold(coordinates protein, int section_start){ //90 anticlock
    std::vector<int> center =protein[section_start];
    for (int i=section_start+1;i<protein.size();i++){
        std::vector<int> current_protein=protein[i];
        int rotatedx=center[0]-(current_protein[1]-center[1]);
        int rotatedy=center[1]+(current_protein[0]-center[0]);
        protein[i]={rotatedx,rotatedy};
    }
    return protein;
}

coordinates down_fold(coordinates protein, int section_start){ //90 clocl
    std::vector<int> center =protein[section_start];
    for (int i=section_start+1;i<protein.size();i++){
        std::vector<int> current_protein=protein[i];
        int rotatedx=center[0]+(current_protein[1]-center[1]);
        int rotatedy=center[1]-(current_protein[0]-center[0]);
        protein[i]={rotatedx,rotatedy};
    }
    return protein;
}

coordinates flip_horizontal_fold(coordinates protein, int section_start){
    std::vector<int> center =protein[section_start];
    for (int i=section_start+1;i<protein.size();i++){
        std::vector<int> current_protein=protein[i];
        int mirroredx=2*center[0]-current_protein[0];
        protein[i]={mirroredx,current_protein[1]};
    }
    return protein;
}
coordinates flip_vertical_fold(coordinates protein, int section_start){
    std::vector<int> center =protein[section_start];
    for (int i=section_start+1;i<protein.size();i++){
        std::vector<int> current_protein=protein[i];
        int mirroredy=2*center[1]-current_protein[1];
        protein[i]={current_protein[0],mirroredy};
    }
    return protein;
}

coordinates HP_model_2D::move(coordinates protein, int& movable){ //we are going to pick at random a monomer, then we will fold the short end of the protein and check if we get any overlaps or too long distances
    int validity_counter=0;
    int direction; 
    bool valid=false;
    while (!valid){
        validity_counter++;
        direction = distrib(gen);
        coordinates new_protein;
        switch (direction) {
                case 0:  // Move left
                    new_protein=up_fold(protein,movable);
                    break;
                case 1:  // Move up
                    new_protein=down_fold(protein,movable);
                    break;
                case 2:  // Move right
                    new_protein=flip_horizontal_fold(protein,movable);
                    break;
                case 3:  // Move down
                    new_protein=flip_vertical_fold(protein,movable);
                    break;

            }
            if (is_valid_move(new_protein,box_size)){
                valid=true;
                return new_protein;
            } 
            if (validity_counter>100) break;
    }
    return protein;
}

double neighbour_distance(std::vector<int>& a,std::vector<int>& b){
    double distance_ab=0;
    for (int i=0; i<b.size();i++){
        distance_ab+=(b[i]-a[i])*(b[i]-a[i]);
    }
    distance_ab = std::sqrt(distance_ab);
    

    return distance_ab;
}

bool is_valid_move_residue(coordinates& new_protein,std::vector<std::vector<int>>& lattice, std::vector<int>& new_position, int& residue, double distance){
    int min_x = lattice[0][0];
    int min_y = lattice[0][1];
    int max_x = lattice[1][0];
    int max_y = lattice[1][1];
    if (new_position[0] < min_x || new_position[0] > max_x || new_position[1] < min_y || new_position[1]> max_y) {
        return false; 
    }

    auto it = std::find(new_protein.begin(), new_protein.end(), new_position);
    if (it != new_protein.end()) {
        return false; 
    }

    if (residue==0){
        if ((neighbour_distance(new_protein[residue+1],new_position)>distance)) return false;
    }
    else if (residue==new_protein.size()-1){
        if ((neighbour_distance(new_protein[residue-1],new_position)>distance)) return false;
    }
    else {
        if ((neighbour_distance(new_protein[residue-1],new_position)>distance) || (neighbour_distance(new_protein[residue+1],new_position)>distance)) return false;
    }
    return true;
}

void HP_model_2D::move_residue(coordinates& protein, int& movable){ 
    int validity_counter=0;
    int direction; 
    bool valid=false;
    while (!valid){
        validity_counter++;
        std::vector<int>  new_position=protein[movable];

        direction = distrib(gen);
        
        switch (direction) {
                case 0:  // Move left
                    new_position[0]-=1;
                    break;
                case 1:  // Move up
                    new_position[1]+=1;
                    break;
                case 2:  // Move right
                    new_position[0]+=1;
                    break;
                case 3:  // Move down
                    new_position[1]-=1;
                    break;
            }
            if (is_valid_move_residue(protein,box_size,new_position, movable, max_distance)){
                valid=true;
                protein[movable]=new_position;
            } 
            if (validity_counter>30){
                valid=true;//maybe if it fails a lot it should reinitilize
            }
            
    }
}

void HP_model_2D::initial_(){
    int direction; 
    bool place=false;
     std::discrete_distribution<int> distrib_bias({1, 2, 1, 1});
    std::vector<int> pos={(box_size[0][0]+box_size[1][0])/2,(-box_size[0][1]+box_size[1][1])/10};
    protein.push_back(pos);
    for (int aa=1;aa<aa_size;aa++){
        int attempt_count = 0;
        place=false;
        while(!place){
            if (attempt_count >= 30) { 
                protein.clear();
                pos = {(box_size[0][0] + box_size[1][0]) / 2, (-box_size[0][1] + box_size[1][1]) / 10};
                protein.push_back(pos);
                aa = 0; 
                attempt_count = 0;
                break; 
            }

            direction = distrib_bias(gen);
            std::vector<int> try_pos=pos;
            switch (direction) {
                    case 0:  // Move left
                        try_pos[0]-=1;
                        break;
                    case 1:  // Move up
                        try_pos[1]+=1;
                        break;
                    case 2:  // Move right
                        try_pos[0]+=1;
                        break;
                    case 3:  // Move 
                        try_pos[1]-=1;
                        break;
                }
            auto it = std::find(protein.begin(), protein.end(), try_pos);
            if (it == protein.end()) {
                place=true;
                protein.push_back(try_pos);
            }
            else {
                attempt_count++; 
            }
        }
        pos=protein[aa];
    }
    double E=conformatio_energy(protein);
    double Comp=compactness(protein);
    print_to_file("output/initial.txt", E, Comp, starting_temperature);
}



void HP_model_2D::print_to_file(std::string path, double E,double C, double T){
     std::ofstream outFile(path);

    outFile << "Energy: " << E << std::endl;
    outFile << "Compactness: " << C << std::endl;
    outFile << "Temperature: " << T << std::endl;

     for (size_t i = 0; i < protein.size(); ++i) {
        const auto& pair = protein[i]; 
        int x = pair[0];
        int y = pair[1];

        if (i < sequence.size()) {
            char letter = sequence[i]; 
            outFile << letter << "\t\t" << x << "\t" << y << std::endl;
        }
    }

    
    outFile.close();
}

int neighbours(coordinates& protein, int& current,std::string sequence ){
    int valid=0;
    int x=protein[current][0];
    int y=protein[current][1];

    int current_HP=HP_translator(sequence[current]);

    std::vector<std::vector<int>> directions = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
    for (auto& direction : directions){
        int checkX=x+direction[0];
        int checkY=y+direction[1];
        std::vector<int> check ={checkX,checkY};
        auto it = std::find(protein.begin(), protein.end(),check );
        if (it != protein.end()) {
            int index = std::distance(protein.begin(), it);
            if ((index != current + 1) && (index != current - 1) && current_HP == 0) { //count hydrophobic
                    valid++;
                }
        }
    }
    return valid;
}

double  HP_model_2D::conformatio_energy(coordinates& new_protein){
    double bond_cout=0;
    for (int i=0;i<protein.size();i++){
        bond_cout+=neighbours(new_protein,i,sequence);
    }
    return -bond_cout/2.0;
}

bool  HP_model_2D::Metropolis_algorithm(double newE, double oldE, double T){
    double deltaE=newE-oldE;
    if (deltaE<0) {
        return true;}
    else{
        double probability = std::exp(-(deltaE)/(T));
        std::uniform_real_distribution<> prob_dis(0.0, 1.0);
        double acceptance=prob_dis(gen);
        return acceptance<probability;
    }

}

void HP_model_2D::run_simulation(){ 
    
    std::uniform_int_distribution<> aa_picker(0,aa_size-1);
    double E,E1;
    int change_counter=0;
    E=conformatio_energy(protein);
    double Comp=compactness(protein);
    min_energy=E; max_compactness=Comp;
    for (int i=0;i<steps;i++){
        int aa=aa_picker(gen);
        move_residue(protein,aa);//this can fail and nothing happens, so it should,'t be a problem too much, otherwise i need to return it here already, which is annoying, but then move can be by reference
        coordinates new_protein=move(protein, aa);
        E1=conformatio_energy(new_protein);
        double temeprature=annealing(i,anneal, rate, anneal_to_zero);
        bool accept=Metropolis_algorithm(E1,E,temeprature);
        if (accept) {
            protein=new_protein;
            E=E1;
            Comp=compactness(new_protein);
            std::string filepath = "output/conformations/" + std::to_string(change_counter) + ".txt";
            print_to_file(filepath, E, Comp,temeprature);
            if (E<min_energy){min_energy=E; best_energy=change_counter;}
            if (Comp>max_compactness){max_compactness=Comp; best_compactness=change_counter;}
            change_counter++;
        }

        
    }
}

double HP_model_2D::compactness(coordinates& new_protein){
    double compactness=0;
    std::string all_are_valid(protein.size(), 'A'); 
     for (int i=0;i<protein.size();i++){
        compactness+=neighbours(new_protein,i,all_are_valid);
    }
    return compactness/2.0;
}

double HP_model_2D::annealing(int& step, bool& isTrue, double rate, bool toZero){ 
    if (!isTrue) return starting_temperature;

    if (toZero){
        rate=(starting_temperature-1e-5)/steps; //maybe dumb to compute rate everytime but ehh
    }

    return starting_temperature-step*rate;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

HP_model_3D::HP_model_3D(std::string  filepath) : input_parameters(filepath), gen(rd()), distrib(0,5){
    aa_size=sequence.size();
    initial_();
}

bool HP_model_3D::is_valid_move(coordinates& new_protein,std::vector<std::vector<int>>& lattice){
    
    int min_x = lattice[0][0];
    int min_y = lattice[0][1];
    int min_z = lattice[0][2];
    int max_x = lattice[1][0];
    int max_y = lattice[1][1];
    int max_z = lattice[1][2];
    

    std::set<std::vector<int>> collected_residues;
    for (auto& residue : new_protein){
        if (residue[0] < min_x || residue[0] > max_x || residue[1] < min_y || residue[1]> max_y || residue[2]< min_z || residue[2]> max_z) {
            return false; 
        }

        if (collected_residues.find(residue) != collected_residues.end()) {
            return false;
        }
        collected_residues.insert(residue);
        
    }
    return true;
}

coordinates z_clock_wise(coordinates protein, int section_start){
    std::vector<int> center =protein[section_start];
    for (int i=section_start+1;i<protein.size();i++){
        std::vector<int> current_protein=protein[i];
        int rotatedx=center[0]+(current_protein[1]-center[1]);
        int rotatedy=center[1]-(current_protein[0]-center[0]);
        protein[i]={rotatedx,rotatedy, current_protein[2]};
    }
    return protein;
}

coordinates z_anticlock_wise(coordinates protein, int section_start){ 
    std::vector<int> center =protein[section_start];
    for (int i=section_start+1;i<protein.size();i++){
        std::vector<int> current_protein=protein[i];
        int rotatedx=center[0]-(current_protein[1]-center[1]);
        int rotatedy=center[1]+(current_protein[0]-center[0]);
        protein[i]={rotatedx,rotatedy, current_protein[2]};
    }
    return protein;
}

coordinates x_clock_wise(coordinates protein, int section_start){ 
    std::vector<int> center =protein[section_start];
    for (int i=section_start+1;i<protein.size();i++){
        std::vector<int> current_protein=protein[i];
        int rotatedy=center[1] + (current_protein[2] - center[2]);
        int rotatedz=center[2] - (current_protein[1] - center[1]);
        protein[i]={current_protein[0],rotatedy, rotatedz};
    }
    return protein;
}

coordinates x_anticlock_wise(coordinates protein, int section_start){ 
    std::vector<int> center =protein[section_start];
    for (int i=section_start+1;i<protein.size();i++){
        std::vector<int> current_protein=protein[i];
        int rotatedy=center[1] - (current_protein[2] - center[2]);
        int rotatedz=center[2] + (current_protein[1] - center[1]);
        protein[i]={current_protein[0],rotatedy, rotatedz};
    }
    return protein;
}
coordinates y_clock_wise(coordinates protein, int section_start){ 
    std::vector<int> center =protein[section_start];
    for (int i=section_start+1;i<protein.size();i++){
        std::vector<int> current_protein=protein[i];
        int rotatedX = center[0] + (current_protein[2] - center[2]);
        int rotatedZ = center[2] - (current_protein[0] - center[0]);
        protein[i] = {rotatedX, current_protein[1], rotatedZ};
    }
    return protein;
}

coordinates y_anticlock_wise(coordinates protein, int section_start){ 
    std::vector<int> center =protein[section_start];
    for (int i=section_start+1;i<protein.size();i++){
        std::vector<int> current_protein=protein[i];
        int rotatedX = center[0] - (current_protein[2] - center[2]);
        int rotatedZ = center[2] + (current_protein[0] - center[0]);
        protein[i] = {rotatedX, current_protein[1], rotatedZ};
    }
    return protein;
}

coordinates flip_yz(coordinates protein, int section_start){
    std::vector<int> center =protein[section_start];
    for (int i=section_start+1;i<protein.size();i++){
        std::vector<int> current_protein=protein[i];
        int mirroredx=2*center[0]-current_protein[0];
        protein[i] = {mirroredx, current_protein[1], current_protein[2]};
    }
    return protein;
}

coordinates flip_xz(coordinates protein, int section_start){
    std::vector<int> center =protein[section_start];
    for (int i=section_start+1;i<protein.size();i++){
        std::vector<int> current_protein=protein[i];
        int mirroredy = 2 * center[1] - current_protein[1];
        protein[i] = {current_protein[0], mirroredy, current_protein[2]}; 
    }
    return protein;
}

coordinates flip_xy(coordinates protein, int section_start){
    std::vector<int> center =protein[section_start];
    for (int i=section_start+1;i<protein.size();i++){
        std::vector<int> current_protein=protein[i];
        int mirroredz=2*center[2]-current_protein[2];
        protein[i] = {current_protein[0], current_protein[1], mirroredz};
    }
    return protein;
}

coordinates HP_model_3D::move(coordinates protein, int& movable){ //we are going to pick at random a monomer, then we will fold the short end of the protein and check if we get any overlaps or too long distances
    int validity_counter=0;
    int direction; 
    bool valid=false;
    while (!valid){
        validity_counter++;
        std::uniform_int_distribution<> distrib9(0,8);
        direction = distrib9(gen);
        coordinates new_protein;
        switch (direction) {
                case 0: 
                    new_protein=x_clock_wise(protein,movable);
                    break;
                case 1: 
                    new_protein=x_anticlock_wise(protein,movable);
                    break;
                case 2:  
                    new_protein=y_clock_wise(protein,movable);
                    break;
                case 3:  
                    new_protein=y_anticlock_wise(protein,movable);
                    break;
                case 4:  
                    new_protein=z_clock_wise(protein,movable);
                    break;
                case 5:  
                    new_protein=z_anticlock_wise(protein,movable);
                    break;
                case 6:  
                    new_protein=flip_yz(protein,movable);
                    break;
                case 7:  
                    new_protein=flip_xz(protein,movable);
                    break;
                case 8:  
                    new_protein=flip_xy(protein,movable);
                    break;

            }
            if (is_valid_move(new_protein,box_size)){
                valid=true;
                return new_protein;
            } 
            if (validity_counter>100) break;
    }
    return protein;
}



bool HP_model_3D::is_valid_move_residue(coordinates& new_protein,std::vector<std::vector<int>>& lattice, std::vector<int>& new_position, int& residue, double distance){
    int min_x = lattice[0][0];
    int min_y = lattice[0][1];
    int min_z = lattice[0][2];
    int max_x = lattice[1][0];
    int max_y = lattice[1][1];
    int max_z = lattice[1][2];
    
    if (new_position[0] < min_x || new_position[0] > max_x || new_position[1] < min_y || new_position[1]> max_y || new_position[2]<min_z || new_position[2]> max_z ) {
        return false; 
    }

    auto it = std::find(new_protein.begin(), new_protein.end(), new_position);
    if (it != new_protein.end()) {
        return false; 
    }

    if (residue==0){
        if ((neighbour_distance(new_protein[residue+1],new_position)>distance)) return false;
    }
    else if (residue==new_protein.size()-1){
        if ((neighbour_distance(new_protein[residue-1],new_position)>distance)) return false;
    }
    else {
        if ((neighbour_distance(new_protein[residue-1],new_position)>distance) || (neighbour_distance(new_protein[residue+1],new_position)>distance)) return false;
    }
    return true;
}

void HP_model_3D::move_residue(coordinates& protein, int& movable){ 
    int validity_counter=0;
    int direction; 
    bool valid=false;
    while (!valid){
        validity_counter++;
        std::vector<int>  new_position=protein[movable];

        direction = distrib(gen);
        
        switch (direction) {
                case 0:  // Move left
                    new_position[0]-=1;
                    break;
                case 1:  // Move up
                    new_position[1]+=1;
                    break;
                case 2:  // Move right
                    new_position[0]+=1;
                    break;
                case 3:  // Move down
                    new_position[1]-=1;
                    break;
                case 4: //move out
                    new_position[2]+=1;
                    break;
                case 5: //move in
                    new_position[2]-=1;
                    break;
            }
            if (is_valid_move_residue(protein,box_size,new_position, movable, max_distance)){
                valid=true;
                protein[movable]=new_position;
            } 
            if (validity_counter>30){
                valid=true;
            }
            
    }
}

void HP_model_3D::initial_(){
    int direction; 
    bool place=false;
    std::discrete_distribution<int> distrib_bias({1, 2, 1, 1,1,1});
    std::vector<int> pos={(box_size[0][0]+box_size[1][0])/2,(-box_size[0][1]+box_size[1][1])/10 , (box_size[0][2]+box_size[1][2])/2};
    protein.push_back(pos);
    for (int aa=1;aa<aa_size;aa++){
        int attempt_count = 0;
        place=false;
        while(!place){
            if (attempt_count >= 30) { 
                protein.clear();
                pos = {(box_size[0][0] + box_size[1][0]) / 2, (-box_size[0][1] + box_size[1][1]) / 10, (box_size[0][2]+box_size[1][2])/2};
                protein.push_back(pos);
                aa = 0; 
                attempt_count = 0;
                break; 
            }

            direction = distrib_bias(gen);
            std::vector<int> try_pos=pos;
            switch (direction) {
                    case 0:  
                        try_pos[0]-=1;
                        break;
                    case 1: 
                        try_pos[1]+=1;
                        break;
                    case 2: 
                        try_pos[0]+=1;
                        break;
                    case 3:  
                        try_pos[1]-=1;
                        break;
                    case 4:  
                        try_pos[2]+=1;
                        break;
                    case 5:  
                        try_pos[2]-=1;
                        break;
                }
            auto it = std::find(protein.begin(), protein.end(), try_pos);
            if (it == protein.end()) {
                place=true;
                protein.push_back(try_pos);
            }
            else {
                attempt_count++; 
            }
        }
        pos=protein[aa];
    }
    double E=conformatio_energy(protein);
    double Comp=compactness(protein);
    print_to_file("output/initial.txt", E, Comp, starting_temperature);
}



void HP_model_3D::print_to_file(std::string path, double E,double C, double T){
     std::ofstream outFile(path);

    outFile << "Energy: " << E << std::endl;
    outFile << "Compactness: " << C << std::endl;
    outFile << "Temperature: " << T << std::endl;

     for (size_t i = 0; i < protein.size(); ++i) {
        const auto& triplet = protein[i];
        int x = triplet[0];
        int y = triplet[1];
        int z = triplet[2];

        
        if (i < sequence.size()) {
            char letter = sequence[i]; //
            outFile << letter << "\t\t" << x << "\t" << y<<"\t"<<z << std::endl;
        }
    }

   
    outFile.close();
}

int HP_model_3D::neighbours(coordinates& protein, int& current,std::string sequence ){
    int valid=0;
    int x=protein[current][0];
    int y=protein[current][1];
    int z=protein[current][2];

    int current_HP=HP_translator(sequence[current]);

    std::vector<std::vector<int>> directions = {{-1, 0,0}, {1, 0,0}, {0, -1,0}, {0, 1,0}, {0,0,1},{0,0,-1}};
    for (auto& direction : directions){
        int checkX=x+direction[0];
        int checkY=y+direction[1];
        int checkZ=z+direction[2];
        std::vector<int> check ={checkX,checkY,checkZ};
        auto it = std::find(protein.begin(), protein.end(),check );
        if (it != protein.end()) {
            int index = std::distance(protein.begin(), it);
            if ((index != current + 1) && (index != current - 1) && current_HP == 0) { //count hydrophobic
                    valid++;
                }
        }
    }
    return valid;
}

double  HP_model_3D::conformatio_energy(coordinates& new_protein){
    double bond_cout=0;
    for (int i=0;i<protein.size();i++){
        bond_cout+=neighbours(new_protein,i,sequence);
    }
    return -bond_cout/2.0;
}

bool  HP_model_3D::Metropolis_algorithm(double newE, double oldE, double T){
    double deltaE=newE-oldE;
    if (deltaE<0) {
        return true;}
    else{
        double probability = std::exp(-(deltaE)/(T));
        std::uniform_real_distribution<> prob_dis(0.0, 1.0);
        double acceptance=prob_dis(gen);
        return acceptance<probability;
    }

}

void HP_model_3D::run_simulation(){ 
    
    std::uniform_int_distribution<> aa_picker(0,aa_size-1);
    double E,E1;
    int change_counter=0;
    E=conformatio_energy(protein);
    double Comp=compactness(protein);
    min_energy=E; max_compactness=Comp;
    for (int i=0;i<steps;i++){
        int aa=aa_picker(gen);
        move_residue(protein,aa);
        coordinates new_protein=move(protein, aa);
        E1=conformatio_energy(new_protein);
        double temperature=annealing(i,anneal, rate, anneal_to_zero);
        bool accept=Metropolis_algorithm(E1,E,temperature);
        if (accept) {
            protein=new_protein;
            E=E1;
            Comp=compactness(new_protein);
            std::string filepath = "output/conformations/" + std::to_string(change_counter) + ".txt"; 
            print_to_file(filepath, E, Comp, temperature);
            if (E<min_energy){min_energy=E; best_energy=change_counter;}
            if (Comp>max_compactness){max_compactness=Comp; best_compactness=change_counter;}
            change_counter++;
        }

        
    }
}

double HP_model_3D::compactness(coordinates& new_protein){
    double compactness=0;
    std::string all_are_valid(protein.size(), 'A'); 
     for (int i=0;i<protein.size();i++){
        compactness+=neighbours(new_protein,i,all_are_valid);
    }
    return compactness/2.0;
}

double HP_model_3D::annealing(int& step, bool& isTrue, double rate, bool toZero){ 
    if (!isTrue) return starting_temperature;

    if (toZero){
        rate=(starting_temperature-1e-5)/steps; //maybe dumb to compute rate everytime but ehh
    }

    return starting_temperature-step*rate;

}

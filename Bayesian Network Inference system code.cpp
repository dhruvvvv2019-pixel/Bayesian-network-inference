#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <map>
#include <iomanip>
#include <cmath>

using namespace std;

class Graph_Node {
private:
    string Node_Name;
    vector<int> Children;
    vector<string> Parents;
    int nvalues;
    vector<string> values;
    vector<float> CPT;

public:
    Graph_Node(string name, int n, vector<string> vals) {
        Node_Name = name;
        nvalues = n;
        values = vals;
    }

    string get_name() {
        return Node_Name;
    }

    vector<int> get_children() {
        return Children;
    }

    vector<string> get_Parents() {
        return Parents;
    }

    vector<float> get_CPT() {
        return CPT;
    }

    int get_nvalues() {
        return nvalues;
    }

    vector<string> get_values() {
        return values;
    }

    void set_CPT(vector<float> new_CPT) {
        CPT.clear();
        CPT = new_CPT;
    }

    void set_Parents(vector<string> Parent_Nodes) {
        Parents.clear();
        Parents = Parent_Nodes;
    }

    int add_child(int new_child_index) {
        for(int i = 0; i < Children.size(); i++) {
            if(Children[i] == new_child_index)
                return 0;
        }
        Children.push_back(new_child_index);
        return 1;
    }
};


class EMStructure{
    public:
    map<string,float> count_value;
 
    map<string, map<string, float>> count_conditional;


};

class network {
    list<Graph_Node> Pres_Graph;

    public:
    vector<int> missing;
    vector<vector<string>> data;
    vector<double> weights;

      // Constructor
    network() {
        missing.clear();
        data.clear();
        weights.clear();
    }

public:
    int addNode(Graph_Node node) {
        Pres_Graph.push_back(node);
        return 0;
    }

    list<Graph_Node>::iterator getNode(int i) {
        int count = 0;
        list<Graph_Node>::iterator listIt;
        for(listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++) {
            if(count++ == i)
                break;
        }
        return listIt;
    }

    int netSize() {
        return Pres_Graph.size();
    }

    int get_index(string val_name) {
        list<Graph_Node>::iterator listIt;
        int count = 0;
        for(listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++) {
            if(listIt->get_name().compare(val_name) == 0)
                return count;
            count++;
        }
        return -1;
    }

    list<Graph_Node>::iterator get_nth_node(int n) {
        list<Graph_Node>::iterator listIt;
        int count = 0;
        for(listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++) {
            if(count == n)
                return listIt;
            count++;
        }
        return listIt;
    }

    list<Graph_Node>::iterator search_node(string val_name) {
        list<Graph_Node>::iterator listIt;
        for(listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++) {
            if(listIt->get_name().compare(val_name) == 0)
                return listIt;
        }
        cout << "node not found: " << val_name << "\n";
        return listIt;
    }

    //function to initialze cpt

    void initialize_cpt(vector<vector<string>> &data){
       
        for(auto& present_node:Pres_Graph){

            vector<float> cpt=present_node.get_CPT();
            int numValues=present_node.get_nvalues();

            int n=cpt.size();

            if(cpt.empty() || cpt[0]==-1){
                cpt.clear();

                int parent_combination=1;

                vector<string> parent=present_node.get_Parents();

                for(auto& parent_node:parent){
                    auto parent_present = search_node(parent_node);
                    int size=parent_present->get_nvalues();
                    parent_combination *= size;
                }


                for(int i=0;i<parent_combination;i++){

                    for(int j=0;j<numValues;j++){

                        cpt.push_back(1.0f/numValues);
                    }

                }

                present_node.set_CPT(cpt);

            }

        }


    }


    //deepseek ----->

int get_CPT_Index(vector<int> &vals, vector<int> &sizes) {

    if(vals.size() == 0)
    { 
        return 0;
    }
    
    int index = 0;
    int b = 1;

    int M = sizes.size();

    for(int i = M-1; i >= 0; i--) {
        index = index + b * vals[i];
        b = b * sizes[i];
    }

    return index;

}

    //Function to implement Expectation ----->

void expectation(vector<vector<string>> &data,vector<EMStructure>& ems){
       
       //clear weights

        weights.clear();
    
        int N=data.size();

        for(int i=0;i<N;i++){

            int miss=-1;

            for(int r=0;r<data[i].size();r++){
                if(data[i][r]=="?"){
                    miss=r;
                    break;
                }
            }

            if(miss==-1){
                //No missing value case 
                weights.push_back(1.0);

                for(int j=0;j<data[i].size();j++){

                    string node_val=data[i][j];
                    auto node=get_nth_node(j);
                    vector<string> parents=node->get_Parents();

                    if(parents.empty()){
                        ems[j].count_value[node_val] += 1.0;
                    }else{

                        string parent_config= "";

                        for(auto &parent_name : parents){

                            int parent_idx=get_index(parent_name);

                            if(!parent_config.empty()) parent_config+="|";
                            parent_config+=data[i][parent_idx];

                        }
                        
                        ems[j].count_conditional[node_val][parent_config] += 1.0;

                    }

                }

            }else{

                //missing value case

                auto missing_node=get_nth_node(miss);
                vector<string> missing_values=missing_node->get_values();
                int num_vals=missing_values.size();

                double denominator=0.0;

                vector<double> all_nums;

                for(int t=0;t<num_vals;t++){

                    double numerator=1.0;
                    vector<string> temp=data[i];
                    temp[miss] = missing_values[t];

                    vector<int> temp_indices;

                    for(int k=0;k<temp.size();k++){
                        auto node_k=get_nth_node(k);

                        vector<string> node_values=node_k->get_values();
                        int idx=-1;

                        for(int m=0;m<node_values.size();m++){
                            
                            if(node_values[m] == temp[k]){
                                idx=m;
                                break;
                            }
                        }

                        temp_indices.push_back(idx);

                    }

                    vector<int> children = missing_node->get_children();

                    for(int child_idx:children){

                        auto child_node=get_nth_node(child_idx);
                        vector<int> vals,sizes;


                        vals.push_back(temp_indices[child_idx]);
                        sizes.push_back(child_node->get_nvalues());

                        vector<string> child_parents = child_node->get_Parents();

                        for(auto & parent_name:child_parents){

                            int parent_idx=get_index(parent_name);
                            vals.push_back(temp_indices[parent_idx]);
                            sizes.push_back(get_nth_node(parent_idx)->get_nvalues());

                        }

                        //calculate P(child / parents) from cpt

                        int cpt_idx= get_CPT_Index(vals,sizes);
                        numerator *= child_node->get_CPT()[cpt_idx];


                    }

                    vector<int> vals,sizes;

                    vals.push_back(t);
                    sizes.push_back(num_vals);

                    vector<string> misssing_parents = missing_node->get_Parents();

                    for(auto & parent_name:misssing_parents){

                        int parent_idx=get_index(parent_name);
                        vals.push_back(temp_indices[parent_idx]);
                        sizes.push_back(get_nth_node(parent_idx)->get_nvalues());

                    }

                    int cpt_idx=get_CPT_Index(vals,sizes);
                    numerator *= missing_node->get_CPT()[cpt_idx];

                    denominator += numerator;
                    all_nums.push_back(numerator);

                }


                for(int t=0;t<num_vals;t++){
                    double w=all_nums[t]/denominator;
                    weights.push_back(w);


                    vector<string> temp2 =data[i];
                    temp2[miss] = missing_values[t];

                    for(int j=0;j<data[i].size();j++){
                        string node_val=temp2[j];

                        auto node=get_nth_node(j);

                        vector<string> node_parents = node->get_Parents();

                        if(node_parents.empty()){
                            ems[j].count_value[node_val] += w;
                        }else{
                            string parent_config="";

                            for(auto &parent_name : node_parents){
                                int parent_idx=get_index(parent_name);
                                if(!parent_config.empty())  parent_config += "|";
                                parent_config += temp2[parent_idx];
                            }
                            ems[j].count_conditional[node_val][parent_config] += w;
                            // cout<<"praetn--->"<<parent_config<<endl;
                        }
                    }
                }




            }


        }



    }



//Function for maximization step ------>

void maximization(bool& maximization_happening,vector<EMStructure> &ems){

        int index=0;

        for(auto& node:Pres_Graph){

            auto parent_node=node.get_Parents();

            if(parent_node.empty()){
                //No parents exist  ------>
                
                float total_count=0.0;
                for(auto node_values:node.get_values()){

                    total_count+=ems[index].count_value[node_values];

                }

                int N=node.get_nvalues();
                vector<float> cpt=node.get_CPT();

                int cpt_index=0;
                for(auto node_values:node.get_values()){
                    cpt[cpt_index]=ems[index].count_value[node_values]/total_count;
                    cpt_index++;
                }

                node.set_CPT(cpt);

            }else {

    // Here Parents Exist
    auto parents = node.get_Parents();
    
    
    int parent_combinations = 1;

    vector<int> parent_sizes;
    for(auto parent : parents) {
        int parent_idx = get_index(parent);
        int parent_nvalues = get_nth_node(parent_idx)->get_nvalues();
        parent_sizes.push_back(parent_nvalues);
        parent_combinations *= parent_nvalues;
    }
    
    vector<float> new_cpt;
    
    for(int comb = 0; comb < parent_combinations; comb++) {
     
vector<string> parent_config_values(parents.size());
int tmp = comb;


for (int p = (int)parents.size() - 1; p >= 0; --p) {
    int parent_val_idx = tmp % parent_sizes[p];
    tmp /= parent_sizes[p];

    int parent_idx = get_index(parents[p]);
    auto parent_node = get_nth_node(parent_idx);
    parent_config_values[p] = parent_node->get_values()[parent_val_idx];
}


        string parent_config_str = "";
        for(auto parent_val : parent_config_values) {
            if(!parent_config_str.empty()) parent_config_str += "|";
            parent_config_str += parent_val;
        }

        // cout<<"ocnfig----->"<<parent_config_str<<endl;
        
        float total_for_this_config = 0.0;
        for(auto node_val : node.get_values()) {
            total_for_this_config += ems[index].count_conditional[node_val][parent_config_str];
        }


        for(auto node_val : node.get_values()) {
            float count = ems[index].count_conditional[node_val][parent_config_str];
            float prob = (count + 0.001) / (total_for_this_config + 0.001 * node.get_nvalues());
            new_cpt.push_back(prob);
        }
    }
    
    node.set_CPT(new_cpt);
}

            index++;

        }

    }

};


// Helper function to trim whitespace
string trim(const string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (string::npos == first) {
        return str;
    }
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, (last - first + 1));
}



network read_network(const char* filename) {
    network BayesNet;
    string line;
    ifstream myfile(filename);
    
    if (!myfile.is_open()) {
        cout << "Error: Could not open file " << filename << endl;
        return BayesNet;
    }

    while (getline(myfile, line)) {
        line = trim(line);
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') {
            continue;
        }

        stringstream ss(line);
        string token;
        ss >> token;

        // Parse variable declarations
        if (token == "variable") {
            string var_name;
            ss >> var_name;
            
            // Read the next line with type and values
            getline(myfile, line);
            stringstream ss2(line);
            
            string type_keyword, discrete_keyword, bracket, equals;
            int num_values;
            ss2 >> type_keyword >> discrete_keyword >> bracket >> num_values >> bracket >> equals >> bracket;
            
            // Read values within brackets
            vector<string> values;
            string value;
            while (ss2 >> value) {
                if (value == "};") break;
                // Remove trailing comma if present
                if (value.back() == ',') {
                    value = value.substr(0, value.length() - 1);
                }
                values.push_back(value);
            }
            
            Graph_Node new_node(var_name, num_values, values);
            BayesNet.addNode(new_node);
        }
        // Parse probability tables
        else if (token == "probability") {
            string paren, node_name;
            ss >> paren >> node_name;
            
            // Handle multi-line probability declarations
            string full_prob_line = line;
            while (full_prob_line.find('{') == string::npos) {
                string next_line;
                if (!getline(myfile, next_line)) break;
                full_prob_line += " " + trim(next_line);
            }
            
            // Extract node name and parents from the probability declaration
            size_t start_paren = full_prob_line.find('(');
            size_t end_paren = full_prob_line.find(')');
            size_t pipe_pos = full_prob_line.find('|');
            
            if (start_paren == string::npos || end_paren == string::npos) {
                continue;
            }
            
            string prob_content = full_prob_line.substr(start_paren + 1, end_paren - start_paren - 1);
            stringstream prob_ss(prob_content);
            
            prob_ss >> node_name;
            
            list<Graph_Node>::iterator listIt = BayesNet.search_node(node_name);
            int index = BayesNet.get_index(node_name);
            
            vector<string> parents;
            
            // Check if there are parent nodes (conditional probability)
            if (pipe_pos != string::npos && pipe_pos < end_paren) {
                string parents_str = full_prob_line.substr(pipe_pos + 1, end_paren - pipe_pos - 1);
                stringstream parent_ss(parents_str);
                string parent;
                
                while (parent_ss >> parent) {
                    // Remove trailing comma if present
                    if (parent.back() == ',') {
                        parent = parent.substr(0, parent.length() - 1);
                    }
                    parents.push_back(parent);
                    
                    list<Graph_Node>::iterator parentIt = BayesNet.search_node(parent);
                    parentIt->add_child(index);
                }
            }
            
            listIt->set_Parents(parents);
            
            // Read CPT values - everything between { and };
            vector<float> cpt;
            bool reading_cpt = false;
            
            while (getline(myfile, line)) {
                line = trim(line);
                if (line == "};") break;
                if (line.empty()) continue;
                
                // Parse the line for probability values
                size_t close_paren = line.find(')');
                string prob_part;
                
                if (close_paren != string::npos) {
                    // Line has conditional values like "(A) 0.5, 0.5;"
                    prob_part = line.substr(close_paren + 1);
                } else if (line.find("table") != string::npos) {
                    // Line has "table" keyword
                    size_t table_pos = line.find("table");
                    prob_part = line.substr(table_pos + 5);
                } else {
                    prob_part = line;
                }
                
                // Extract all probability values from the line
                stringstream ss_prob(prob_part);
                string token;
                while (ss_prob >> token) {
                    // Remove trailing punctuation
                    while (!token.empty() && (token.back() == ',' || token.back() == ';')) {
                        token = token.substr(0, token.length() - 1);
                    }
                    
                    // Check if it's a valid number
                    if (!token.empty() && (isdigit(token[0]) || token[0] == '.' || token[0] == '-')) {
                        cpt.push_back(atof(token.c_str()));
                    }
                }
            }
            
            listIt->set_CPT(cpt);
        }
    }
    
    myfile.close();
    return BayesNet;
}


void write_network(const char* filename, network& BayesNet) {
    ofstream outfile(filename);

    if (!outfile.is_open()) {
        cout << "Error: Could not open file " << filename << " for writing" << endl;
        return;
    }

    outfile << "// Bayesian Network" << endl << endl;

    int N = BayesNet.netSize();

    // Write all nodes first
    for (int i = 0; i < N; i++) {
        auto node = BayesNet.get_nth_node(i);

        outfile << "variable " << node->get_name() << " {" << endl;
        outfile << "  type discrete [ " << node->get_nvalues() << " ] = { ";

        vector<string> vals = node->get_values();
        for (int j = 0; j < (int)vals.size(); j++) {
            outfile << vals[j];
            if (j < (int)vals.size() - 1) outfile << ", ";
        }
        outfile << " };" << endl;
        outfile << "}" << endl;
    }

    // Write probability tables
    outfile << std::fixed << std::setprecision(6); // consistent float formatting
    for (int i = 0; i < N; i++) {
        auto node = BayesNet.get_nth_node(i);
        vector<string> parents = node->get_Parents();
        vector<string> values = node->get_values();
        vector<float> cpt = node->get_CPT();

        // Header
        outfile << "probability ( " << node->get_name();
        if (!parents.empty()) {
            outfile << " | ";
            for (int j = 0; j < (int)parents.size(); j++) {
                outfile << parents[j];
                if (j < (int)parents.size() - 1) outfile << ", ";
            }
        }
        outfile << " ) {" << endl;

        // compute parent radices (number of values for each parent)
        vector<int> radices;
        radices.reserve(parents.size());
        for (auto &pname : parents) {
            auto pnode = BayesNet.search_node(pname);
            radices.push_back(pnode->get_nvalues());
        }

        // product of radices
        int parent_combinations = 1;
        for (int r : radices) parent_combinations *= r;

        // Write the CPT values. Keep a single cpt_index that we advance in the same order read_network used
        int cpt_index = 0;

        if (parents.empty()) {
            // unconditional: use table
            outfile << "    table ";
            for (int k = 0; k < (int)values.size(); k++) {
                if (cpt_index < (int)cpt.size()) outfile << cpt[cpt_index++];
                else outfile << "-1";
                if (k < (int)values.size() - 1) outfile << ", ";
            }
            outfile << ";" << endl;
        } else {
            // conditional: for each parent combination write "( v1, v2, ... ) p1, p2, ...;"
            for (int comb = 0; comb < parent_combinations; comb++) {
                // compute mixed-radix digits (indices for each parent),
                // with the rightmost parent cycling fastest (so we fill indices from right->left)
                vector<int> idx(parents.size(), 0);
                int tmp = comb;
                for (int p = (int)parents.size() - 1; p >= 0; p--) {
                    idx[p] = tmp % radices[p];
                    tmp /= radices[p];
                }

                // print parent values in the original parent order
                outfile << "    ( ";
                for (int p = 0; p < (int)parents.size(); p++) {
                    auto pnode = BayesNet.search_node(parents[p]);
                    auto pvals = pnode->get_values();
                    int vidx = idx[p];
                    outfile << pvals[vidx];
                    if (p < (int)parents.size() - 1) outfile << ", ";
                }
                outfile << " ) ";

                // print the probabilities for this parent combination (one probability per child value)
                for (int k = 0; k < (int)values.size(); k++) {
                    if (cpt_index < (int)cpt.size()) outfile << cpt[cpt_index++]; // round to last non zero decimal valu
                    else outfile << "-1";
                    if (k < (int)values.size() - 1) outfile << ", ";
                }
                outfile << ";" << endl;
            }
        }

        // ***** THIS IS THE FIX *****
        outfile << "};" << endl << endl;
    }

    outfile.close();
    cout << "Network written to file: " << filename << endl;
}



//This function is used to read from record.dat and return vector<vector<string>> data

vector<vector<string>> read_data(string filename) {

    ifstream myfile(filename);
    string line;
    vector<vector<string>> data;

    if (!myfile.is_open()) {
        cout << "Error : Cannot Open File " << filename << endl;
        return data;
    }

    while (getline(myfile, line)) {
       
        if (line.empty()) continue;
        
        stringstream ss(line);
        string val;
        vector<string> record;

        while (getline(ss, val, ',')) {
          
            val = trim(val);
            
            if (val.size() >= 2 && val[0] == '"' && val[val.size()-1] == '"') {
                val = val.substr(1, val.size() - 2);
            }  
            record.push_back(val);
        }
        
        data.push_back(record);
    }


    myfile.close();


    cout << "Loaded " << data.size() << " records from " << filename << endl;
    return data;
}


#ifndef BN_LIB
int main(int argc, char* argv[]) {
     if(argc < 3) {
        cout << "Usage: " << argv[0] << " <network.bif> <data.csv>" << endl;
        return 1;
    }

    cout << "Reading network from: " << argv[1] << endl;
    network BayesNet = read_network(argv[1]);
    
    cout << "Network loaded successfully!" << endl;
    cout << "Number of nodes: " << BayesNet.netSize() << endl;

    if(BayesNet.netSize() == 0) {
        cout << "ERROR: No nodes loaded from network file!" << endl;
        return 1;
    }

    cout << "Reading data from: " << argv[2] << endl;
    vector<vector<string>> data = read_data(argv[2]);
    cout << "Loaded " << data.size() << " records" << endl;

    if(data.empty()) {
        cout << "ERROR: No data loaded!" << endl;
        return 1;
    }


   BayesNet.missing.clear();
    for(int i = 0; i < data.size(); i++) {
        int miss_idx = -1;
        for(int j = 0; j < data[i].size(); j++) {
            if(data[i][j] == "?") {
                miss_idx = j;
                break;
            }
        }
        BayesNet.missing.push_back(miss_idx);
    }

     // Initialize CPTs if needed
    BayesNet.initialize_cpt(data);

    // Run EM algorithm
    vector<EMStructure> ems(BayesNet.netSize());
    bool maximization_happening = true;
    
    for(int iteration = 0; iteration < 10; iteration++) {
        cout << "Iteration " << iteration + 1 << ":" << endl;
        
for (auto &e : ems) {
    e.count_value.clear();
    e.count_conditional.clear();
}

     //Expectation step

        BayesNet.expectation(data, ems);
        cout << "E-step completed" << endl;
        
        // Maximization step  
        BayesNet.maximization(maximization_happening, ems);
        cout << "M-step completed" << endl;
        
       
    }

    // Write final network
    string output_file = "solved_" + string(argv[1]);
    write_network(output_file.c_str(), BayesNet);
    
    cout << "EM completed! Output written to " << output_file << endl;
    
    return 0;


}
#endif // BN_LIB




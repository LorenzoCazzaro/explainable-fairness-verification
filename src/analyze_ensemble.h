#include "analyze_ensemble_utils.h"
#include <iterator>
#include <queue>
#include <thread>
#include <mutex> 
#include <chrono>
#include <string.h>
#define SYMATKCOMPARATOR priority_pairs_cut_ind
#define MIN_SYMATK_TO_BE_SPLITTED 10

void execute_iters_analyze_ensemble(set<SymbolicAttack, SYMATKCOMPARATOR>& C, set<SymbolicAttack>& E, Attacker A, Ensemble &T, size_t max_iter) {
    
    size_t iter = 0;

    while (C.size() && iter < max_iter) {         
        size_t k = 0;
        if(C.size() > MIN_SYMATK_TO_BE_SPLITTED)
            k = std::max(MIN_SYMATK_TO_BE_SPLITTED, int(C.size()*0.05));
        else k = C.size();

        vector<SymbolicAttack> to_be_splitted;
        auto C_iter = C.begin();
        size_t j = 0;
        while(j < k) {
            to_be_splitted.push_back(*C_iter);
            C_iter++;
            j++;
        }
      
        C.erase(C.begin(), C_iter);
        
        size_t i = 0;

        while (i < k) {
            const SymbolicAttack u = to_be_splitted[i];
            std::vector<Tree*> uT;
            validate_s2(u, T, &uT);
            std::set<SymbolicAttack> u_splitted;
            split_sa(u, u_splitted, uT, A);

            if (u_splitted.size()) {
                for(const SymbolicAttack& u_splt : u_splitted) {
                    size_t type = validate_s2(u_splt, T);
                    
                    if (type==1)
                        E.insert(u_splt);
                    else if (type > 1)
                        C.insert(u_splt); 
                }
            } else {
                E.insert(u);
            }
            i++;
        }
        ++iter;
    }
}

void execute_iters_analyze_ensemble_with_metrics_computations(set<SymbolicAttack, SYMATKCOMPARATOR>& C, Ensemble &T, Attacker& A, size_t max_iter, VV_f& data, V_f& labels, vector<bool>& robustness_x_instance, vector<vector<bool>>& resilience_x_instance, vector<HyperRectangle>& Ns, double& analyzer_time) {
    size_t iter = 0;
    size_t n_iter_clear_E = 1;

    set<SymbolicAttack> E;

    while (C.size() && iter < max_iter) {
        auto start = std::chrono::high_resolution_clock::now();
        size_t k;
        if(C.size() > MIN_SYMATK_TO_BE_SPLITTED)
            k = std::max(MIN_SYMATK_TO_BE_SPLITTED, int(C.size()*0.05));
        else k = C.size();

        vector<SymbolicAttack> to_be_splitted;
        auto C_iter = C.begin();
        size_t j = 0;
        while(j < k) {
            to_be_splitted.push_back(*C_iter);
            C_iter++;
            j++;
        }
      
        C.erase(C.begin(), C_iter);
        
        size_t i = 0;

        while (i < k) {
            const SymbolicAttack u = to_be_splitted[i];
            std::vector<Tree*> uT;
            validate_s2(u, T, &uT);
            std::set<SymbolicAttack> u_splitted;
            split_sa(u, u_splitted, uT, A);

            if (u_splitted.size()) {
                for(const SymbolicAttack& u_splt : u_splitted) {
                    size_t type = validate_s2(u_splt, T);
                    if (type==1)
                        E.insert(u_splt);
                    else if (type > 1)
                        C.insert(u_splt); 
                }
            } else {
                E.insert(u);
            }
            i++;
        }
        ++iter;
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> dur = end - start;
        analyzer_time += dur.count();

        if(iter % n_iter_clear_E == 0) {
           for(size_t i = 0; i < data.size(); i++) {
                if(robustness_x_instance[i]) {
                    float pred = T.predict(data[i]);
                    robustness_x_instance[i] = (pred == labels[i]) && !check_in2(data[i], E, T, labels[i]);
                    for(size_t k = 0; k < Ns.size(); k++) {
                        if(resilience_x_instance[k][i]) {
                            HyperRectangle N_data = data[i] + Ns[k];
                            resilience_x_instance[k][i] = robustness_x_instance[i] && check_inter(N_data, E);
                        }
                    }
                }
            }
            E.clear();   
        }
    }
}

template <class Comp>
void compute_rob_res_data(std::set<SymbolicAttack, Comp>& Cp, VV_f& data, V_f& labels, vector<HyperRectangle>& Ns, Ensemble &T, vector<bool>& robustness_x_instance, vector<vector<bool>>& resilience_x_instance) {
    
    std::thread::id this_id = std::this_thread::get_id();

    for(size_t i = 0; i < data.size(); i++){

        if(robustness_x_instance[i]) {
            float pred = T.predict(data[i]);
            robustness_x_instance[i] = (pred == labels[i]) && !check_in2(data[i], Cp, T, labels[i]);
            for(size_t k = 0; k < Ns.size(); k++) {
                if(resilience_x_instance[k][i]) {
                    HyperRectangle N_data = data[i] + Ns[k];
                    resilience_x_instance[k][i] = robustness_x_instance[i] && check_inter(N_data, Cp);
                }
            }
        }
    }
}

char* print_time() {
    std::chrono::system_clock::time_point today = std::chrono::system_clock::now();
    time_t tt = std::chrono::system_clock::to_time_t ( today );
    char* time_string = ctime(&tt);
    time_string[strlen(time_string)-1] = '\0';
    return time_string;
}


void analyze_ensemble(Ensemble &T, Attacker& A, size_t n_threads, size_t max_iter, size_t iterations_x_log, set<SymbolicAttack>& U){

    std::cout << "[" << print_time() << "]\n" << "Analyze tree ensemble stability" <<  std::endl;
    
    ens = T;
    std::set<SymbolicAttack> C, E, C2;

    for (Tree* t : T.get_trees())
        analyze_tree(C, *t, A);
    
    std::cout << "C size after analyzying trees: " << C.size() << std::endl; 

    for (const SymbolicAttack u : C) {
        size_t type = validate_s2(u, T);
        if (type==1)
            E.insert(u);   
        if (type > 1)
            C2.insert(u);
    }
    
    C = C2;
    std::cout << "C size after first cleaning: " << C.size() << std::endl;
    C2.clear();

    std::vector<std::set<SymbolicAttack, SYMATKCOMPARATOR>> Cs(n_threads);
    std::vector<std::set<SymbolicAttack>> Es(n_threads);

    bool C_not_empty = C.size() != 0;

    size_t i = 0;
    for(auto iter = C.begin(); iter != C.end(); iter++) {
        i %= n_threads;
        Cs[i].insert(*iter);
        i++;
    } 
    C.clear();
    
    for(auto iter = E.begin(); iter != E.end(); iter++) {
        i %= n_threads;
        Es[i].insert(*iter);
        i++;
    } 

    size_t iter = 0;

    cout << "ITERATIONS: " << iter << endl;
    for(size_t i = 0; i < n_threads; i++){
        cout << "THREAD " << i << ": " << "size C: " << Cs[i].size() << ", size E: " << Es[i].size() << endl; 
    }
    cout << endl;
        
    while(iter < max_iter && C_not_empty) {

        std::vector<std::thread> threads(n_threads);

        for(size_t i = 0; i < n_threads; i++) {
            threads[i] = thread(execute_iters_analyze_ensemble, std::ref(Cs[i]), std::ref(Es[i]), A, std::ref(T), iterations_x_log);
        }
        
        for(size_t i = 0; i < n_threads; i++)
            threads[i].join();

        iter += iterations_x_log;

        cout << "ITERATIONS: " << iter << endl;
        for(size_t i = 0; i < n_threads; i++){
            cout << "THREAD " << i << ": " << "size C: " << Cs[i].size() << ", size E: " << Es[i].size() << endl; 
        }
        cout << endl;

        size_t i = 0;
        C_not_empty = false;
        while(!C_not_empty && i < n_threads) {
            if(Cs[i].size())
                C_not_empty = true;
            i++;
        }
    }

    for(size_t i = 0; i < n_threads; i++) {
        U.insert(Cs[i].begin(), Cs[i].end());
        U.insert(Es[i].begin(), Es[i].end());
    }
}

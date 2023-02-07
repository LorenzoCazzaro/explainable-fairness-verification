#include <cstdlib>
#include <iterator>
#include <queue>
#include <algorithm>
#include <thread>
#include <mutex> 
#include <chrono>
#include <string.h>
#define SYMATKCOMPARATOR priority_pairs_cut_ind
#define MIN_SYMATK_TO_BE_SPLITTED 10
#define LOGGING_OUTPUT 0
#define LOGGING_OUTPUT_DUMP_END_THREAD 0

void analyze_ensemble_with_priority_queue(std::set<SymbolicAttack, SYMATKCOMPARATOR>& U, std::set<SymbolicAttack>& U_ended, Attacker& A, Ensemble &T, size_t max_iter, size_t& S_elem_counter, double& analyzer_time)
{
    size_t iter = 0;
    std::thread::id this_id = std::this_thread::get_id();
    set<SymbolicAttack> S;

    while (U.size() && iter < max_iter) {
        auto start = std::chrono::high_resolution_clock::now();
        #if LOGGING_OUTPUT == 1
            if(iter%1 == 0) 
            std::cout << "[" << print_time() << "]" << " THREAD ID: " << this_id << ", ITER: " << iter << " (" << U.size() << ", " << U_ended.size() << ", " << S_elem_counter << ")" << std::endl;  
        #endif

        size_t k;
        if(U.size() > MIN_SYMATK_TO_BE_SPLITTED)
            k = std::max(MIN_SYMATK_TO_BE_SPLITTED, int(U.size()*0.05));
        else k = U.size();

        vector<SymbolicAttack> to_be_splitted;
        auto U_iter = U.begin();
        size_t j = 0;
        while(j < k) {
            to_be_splitted.push_back(*U_iter);
            U_iter++;
            j++;
        }
      
        U.erase(U.begin(), U_iter);
        
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
                    if (type == 0)
                        S.insert(u_splt);
                    if (type==1)
                        U_ended.insert(u_splt);
                    else U.insert(u_splt); 
                }
                //cout << endl;
            } else {
                U_ended.insert(u);
            }
            i++;
        }
        ++iter;
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> dur = end - start;
        analyzer_time += dur.count();
        S_elem_counter += S.size();
        S.clear();
    }

    #if LOGGING_OUTPUT_DUMP_END_THREAD==1
        std::cout << "[" << print_time() << "]" << " THREAD ID: " << this_id << ", END ITER: " << iter << " (" << U.size() << ", " << U_ended.size() << ", " << S_elem_counter << ")" << std::endl;  
    #endif
}

void analyze(set<SymbolicAttack>& U, Attacker& A, Ensemble &T, size_t n_threads, size_t max_iter, size_t iterations_x_dump, ofstream* ofs){

    VM_fst ths_occ_aly = T.get_ths_occ();

    size_t n_f_v_p = 0;
    size_t n_feat = ths_occ_aly.size();
    std::set<float> v_set;
    
    for (auto fv_set : ths_occ_aly) {
        for (auto& v_occ : fv_set) {
            float vv = v_occ.first;
            size_t oc = v_occ.second;
            n_f_v_p = n_f_v_p + oc;
            //std::cout << vv << std::endl;
            v_set.insert(vv);
        }
    }

    size_t n_ths = v_set.size();

    std::cout << "#atk dim: " << A.get_dim() << std::endl;
    std::cout << "#features: " << n_feat << std::endl;
    std::cout << "#thresholds: " << n_ths << std::endl;
    std::cout << "#feat thres pairs: " << n_f_v_p << std::endl;

    char* buffer = new char[20000];

    std::cout << "[" << print_time() << "]" << " analyze_ensemble_priority_queue_for_splitting with dumps" <<  std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();

    double total_time_ms = 0;
    ens = T;
    size_t iter = 0;

    std::set<SymbolicAttack> Up;
    std::set<SymbolicAttack> U_ended;
    std::set<SymbolicAttack> S;

    for (Tree* t : T.get_trees())
        analyze_tree(U, *t, A);
    std::cout << "U size after analyze_tree: " << U.size() << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> dur = end - start;
    total_time_ms += dur.count();
    cout << "[" << print_time() << "] " << "Tempo per analyze single trees: " << dur.count() << endl;
    sprintf(buffer, "[%s] Tempo per analyze single trees: %f\n", print_time(), dur.count());
    ofs->write(buffer, strlen(buffer));
    ofs->flush();

    start = std::chrono::high_resolution_clock::now();
    for (const SymbolicAttack u : U) {
        size_t type = validate_s2(u, T);
        if (type==1)
            U_ended.insert(u);   
        if (type > 1)
            Up.insert(u);
    }
    
    U = Up;
    std::cout << "U size after first cleaning: " << U.size() << std::endl;
    std::cerr << Up.size() << " UP1" << std::endl;
    Up.clear();
    std::cerr << Up.size() << " UP2" << std::endl;

    end = std::chrono::high_resolution_clock::now();
    dur = end - start;
    total_time_ms += dur.count();

    cout << endl << endl;

    cout << endl << endl;
    cout << "MAX ITER: " << max_iter << endl;
    cout << "ITERATIONS X DUMP: " << iterations_x_dump << endl;
    cout << "INITIAL U DIM: " << U.size() << endl;
    cout << "INITIAL U_ENDED DIM: " << U_ended.size() << endl;
    cout << "INITIAL S DIM: " << S.size() << endl;
    sprintf(buffer, "[%s] INITIAL U DIM: %lu\n", print_time(), U.size());
    cout << buffer << endl;
    ofs->write(buffer, strlen(buffer));
    ofs->flush();
    sprintf(buffer, "[%s] INITIAL U_ENDED DIM: %lu\n", print_time(), U_ended.size());
    cout << buffer << endl;
    ofs->write(buffer, strlen(buffer));
    ofs->flush();
    sprintf(buffer, "[%s] INITIAL S DIM: %lu\n", print_time(), S.size());
    cout << buffer << endl;
    ofs->write(buffer, strlen(buffer));
    ofs->flush();
    sprintf(buffer, "[%s] Tot iter: %zu\n", print_time(), iter);
    cout << buffer << endl;
    ofs->write(buffer, strlen(buffer));
    ofs->flush();
    cout << "Total time: " << total_time_ms << endl;
    sprintf(buffer, "[%s] Total time: %f\n", print_time(), total_time_ms);
    ofs->write(buffer, strlen(buffer));
    ofs->flush();

    start = std::chrono::high_resolution_clock::now();
    std::vector<std::set<SymbolicAttack, SYMATKCOMPARATOR>> Us(n_threads);
    std::vector<std::set<SymbolicAttack>> U_endeds(n_threads);

    size_t i = 0;
    bool U_not_empty = U.size() != 0;
	cout << "U è vuoto? La size è " << U.size() << endl;

    for(auto iter = U.begin(); iter != U.end(); iter++) {
        Us[i].insert(*iter);
        i++;
        i = i%n_threads;
    }

    U.clear();

    i = 0;
    for(auto iter = U_ended.begin(); iter != U_ended.end(); iter++) {
        U_endeds[i].insert(*iter);
        i++;
        i = i%n_threads;
    } 
    U_ended.clear();

    end = std::chrono::high_resolution_clock::now();
    dur = end - start;
    total_time_ms += dur.count();

    while(iter < max_iter && U_not_empty) {
        size_t summ_cums = 0;
        for (size_t ixx = 0 ; ixx < n_threads; ++ixx) {
            summ_cums += Us[ixx].size();
        }
        std::cout << "U size after iter: " << iter << " " << summ_cums << std::endl;

        start = std::chrono::high_resolution_clock::now();
        std::vector<std::thread> threads(n_threads);

        cout << "ITERATIONS OF THE PARALLEL ALG: " << iter << endl;

        vector<double> total_times_analyzer_x_threads(n_threads, 0);
        vector<size_t> S_elem_counter_x_threads(n_threads, 0);

        end = std::chrono::high_resolution_clock::now();
        dur = end - start;
        total_time_ms += dur.count();

        for(size_t i = 0; i < n_threads; i++) {
            threads[i] = thread(analyze_ensemble_with_priority_queue, std::ref(Us[i]), std::ref(U_endeds[i]), std::ref(A), std::ref(T), iterations_x_dump, std::ref(S_elem_counter_x_threads[i]), std::ref(total_times_analyzer_x_threads[i]));
        }
        
        for(size_t i = 0; i < n_threads; i++)
            threads[i].join();

        iter += iterations_x_dump;

        total_time_ms += *(max_element(total_times_analyzer_x_threads.begin(), total_times_analyzer_x_threads.end()));

        U_not_empty = false;
        while(!U_not_empty && i < n_threads) {
            if(Us[i].size())
                U_not_empty = true;
            i++;
        }
    }

    for(size_t i = 0; i < n_threads; i++) {
        U.insert(Us[i].begin(), Us[i].end());
        U.insert(U_endeds[i].begin(), U_endeds[i].end());
    }
}

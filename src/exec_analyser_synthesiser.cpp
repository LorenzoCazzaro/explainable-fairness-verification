#include "synthetiser.hpp"
#include "loader.h"
//#include "analyzer.h"
#include "analyze_ensemble.h"
//#include "utils.hpp"

#define UNION 0
#define CONDITIONS_OUTPUT 0
#define REDUCED_FEATURE_SPACE 1
#define COMPACT_CONDITIONS 1

using namespace jsoncons;

void single_dim_hyperrectangle_union(unordered_set<HyperRectangle> &hyperrects)
{

    auto iter = hyperrects.begin();
    while (iter != hyperrects.end())
    {
        auto iter_union = iter;
        iter_union++;
        bool unified = false;
        bool result;
        HyperRectangle res;
        while (!unified && iter_union != hyperrects.end())
        {
            result = iter->single_dim_union(*iter_union, res);
            if (result) {
                unified = true;
            }
            else
                iter_union++;
        }
        if (unified)
        {
            hyperrects.erase(iter);
            hyperrects.erase(iter_union);
            hyperrects.insert(res);
            iter = hyperrects.begin();
        }
        else
            iter++;
    }
}

void single_dim_hyperrectangle_union(set<HyperRectangle> &hyperrects)
{
    auto iter = hyperrects.begin();
    while (iter != hyperrects.end())
    {
        auto iter_union = iter;
        iter_union++;
        bool unified = false;
        bool result;
        HyperRectangle res;
        while (!unified && iter_union != hyperrects.end())
        {
            result = iter->single_dim_union(*iter_union, res);
            if (result)
                unified = true;
            else
                iter_union++;
        }
        if (unified)
        {
            hyperrects.erase(iter);
            hyperrects.erase(iter_union);
            hyperrects.insert(res);
            iter = hyperrects.begin();
        }
        else
            iter++;
    }
}

void remove_prot_attribute_interval_hyperrect(size_t protected_attribute_index, HyperRectangle &hyper, unordered_set<HyperRectangle> &set_to_fill)
{
    hyper[protected_attribute_index] = Interval(-INF, INF, 0, 0);
    set_to_fill.insert(hyper);
}

void remove_prot_attribute_interval_hyperrects(size_t protected_attribute_index, vector<SymbolicAttack>::iterator hypers_to_convert_begin, vector<SymbolicAttack>::iterator hypers_to_convert_end, unordered_set<HyperRectangle> &set_to_fill)
{
    while (hypers_to_convert_begin != hypers_to_convert_end)
    {
        HyperRectangle pre_to_insert = hypers_to_convert_begin->get_pre();
        remove_prot_attribute_interval_hyperrect(protected_attribute_index, pre_to_insert, set_to_fill);
        hypers_to_convert_begin++;
    }
}

bool check_hyperrect_validity_x_categorical_features(HyperRectangle &hyper_to_check, vector<vector<size_t>> &categorical_values_columns_index)
{
    size_t i = 0;
    bool valid_hyper = true;
    while (valid_hyper && i < categorical_values_columns_index.size())
    {
        size_t j = 0;
        ushort n_values_1 = 0;
        ushort n_values_0 = 0;
        while (j < categorical_values_columns_index[i].size() && n_values_1 < 2)
        {
            const Interval &intv = hyper_to_check[categorical_values_columns_index[i][j]];
            if (intv.get_left() == 0.5)
                ++n_values_1;
            else if (intv.get_right() == 0.5)
                ++n_values_0;
            j++;
        }
        if (n_values_1 > 1) {
            valid_hyper = false;
            //cout << "CONDIZIONE FORTUNATA!!!!" << endl;
        }
        else if(n_values_0 == categorical_values_columns_index[i].size()) {
            //cout << "CONDIZIONE SFIGATA!!!!" << endl;
            valid_hyper = false;
        }
        i++;
    }

    return valid_hyper;
}

void complete_condition(HyperRectangle cond_to_complete, Condition& result, size_t feature_space_dims, map<size_t, size_t>& numerical_values_mapping, map<size_t, size_t>& categorical_values_mapping, vector<vector<size_t>>& categorical_indexes_groups)
{
    const vector<Interval>& intvs = cond_to_complete.get_intvs_ref();
    //assume result is an empty condition
    result = Condition(feature_space_dims);
    size_t index = 0;
    while(index < intvs.size()) {
        bool is_left_inf = intvs[index].get_left() == -INF;
        bool is_right_inf = intvs[index].get_right() == INF;
        if(!is_left_inf || !is_right_inf) {
            if(numerical_values_mapping.find(index) != numerical_values_mapping.end()) {
                result.insert_atomic_condition(numerical_values_mapping[index], Interval((is_left_inf) ? 0 : intvs[index].get_left(), (is_right_inf) ? 1 : intvs[index].get_right(), (is_left_inf) ? 0 : intvs[index].get_left_type(), (is_right_inf) ? 0 : intvs[index].get_right_type()));
                index++;
            } else {
                bool found = false;
                auto iter = categorical_indexes_groups.begin();
                bool is_1_value_in = false;
                size_t index_is_1_value_in = 0;
                while(!found && iter != categorical_indexes_groups.end()) {
                    if(find(iter->begin(), iter->end(), categorical_values_mapping[index]) != iter->end()) {
                        found = true;
                    } else iter++;
                }
                if(intvs[index].get_left() == 0.5) {
                    index_is_1_value_in = index;
                    is_1_value_in = true;
                }
                vector<size_t> indexes_of_same_cat_variable;
                indexes_of_same_cat_variable.push_back(index);
                index++;
                while(index < intvs.size() && find(iter->begin(), iter->end(), categorical_values_mapping[index]) != iter->end()) {
                    if(intvs[index].get_left() != -INF || intvs[index].get_right() != INF) {
                        if(!is_1_value_in && intvs[index].get_left() == 0.5) {
                            is_1_value_in = true;
                            index_is_1_value_in = index;
                        }   
                        indexes_of_same_cat_variable.push_back(index);
                    }
                    index++;
                }
                if(is_1_value_in) {
                    result.insert_atomic_condition(categorical_values_mapping[index_is_1_value_in], Interval(0.5, 1, 1, 0));
                    for(auto cat_index_in_group = iter->begin(); cat_index_in_group != iter->end(); cat_index_in_group++) {
                        if(*cat_index_in_group != categorical_values_mapping[index_is_1_value_in]) {
                            result.insert_atomic_condition(*cat_index_in_group, Interval(0, 0.5, 0, 0));
                        }
                    }
                } else {
                    
                    for(auto indexes = indexes_of_same_cat_variable.begin(); indexes != indexes_of_same_cat_variable.end(); indexes++) {
                        result.insert_atomic_condition(categorical_values_mapping[*indexes], Interval(0, 0.5, 0, 0));
                    }
                }
            }
        } else {
            index++;
        }
    }
}

void complete_condition_compact(HyperRectangle cond_to_complete, Condition& result, size_t feature_space_dims, map<size_t, size_t>& numerical_values_mapping, map<size_t, size_t>& categorical_values_mapping, vector<vector<size_t>>& categorical_indexes_groups) {
    const vector<Interval>& intvs = cond_to_complete.get_intvs_ref();
    //assume result is an empty condition
    result = Condition(feature_space_dims);
    size_t index = 0;
    while(index < intvs.size()) {
        bool is_left_inf = intvs[index].get_left() == -INF;
        bool is_right_inf = intvs[index].get_right() == INF;
        if(!is_left_inf || !is_right_inf) {
            if(numerical_values_mapping.find(index) != numerical_values_mapping.end()) {
                result.insert_atomic_condition(numerical_values_mapping[index], Interval((is_left_inf) ? 0 : intvs[index].get_left(), (is_right_inf) ? 1 : intvs[index].get_right(), (is_left_inf) ? 0 : intvs[index].get_left_type(), (is_right_inf) ? 0 : intvs[index].get_right_type()));
                index++;
            } else {
                bool found = false;
                auto iter = categorical_indexes_groups.begin();
                bool is_1_value_in = false;
                size_t index_is_1_value_in = 0;
                while(!found && iter != categorical_indexes_groups.end()) {
                    if(find(iter->begin(), iter->end(), categorical_values_mapping[index]) != iter->end()) {
                        found = true;
                    } else iter++;
                }
                vector<size_t> indexes_of_same_cat_variable_with_0;
                if(intvs[index].get_left() == 0.5) {
                    index_is_1_value_in = index;
                    is_1_value_in = true;
                } else indexes_of_same_cat_variable_with_0.push_back(categorical_values_mapping[index]);
                index++;
                while(index < intvs.size() && find(iter->begin(), iter->end(), categorical_values_mapping[index]) != iter->end()) {
                    if(intvs[index].get_left() != -INF || intvs[index].get_right() != INF) {
                        if(!is_1_value_in && intvs[index].get_left() == 0.5) {
                            is_1_value_in = true;
                            index_is_1_value_in = index;
                        } else indexes_of_same_cat_variable_with_0.push_back(categorical_values_mapping[index]);
                    }
                    index++;
                }
                if(is_1_value_in) {
                    result.insert_atomic_condition(categorical_values_mapping[index_is_1_value_in], Interval(0.5, 1, 1, 0));
                } else {
                    for(auto iter_indexes_cat_group = iter->begin(); iter_indexes_cat_group != iter->end(); iter_indexes_cat_group++) {
                        if(find(indexes_of_same_cat_variable_with_0.begin(), indexes_of_same_cat_variable_with_0.end(), *iter_indexes_cat_group) == indexes_of_same_cat_variable_with_0.end())
                            result.insert_atomic_condition(*iter_indexes_cat_group, Interval(0.5, 1, 1, 0));
                    }
                }
            }
        } else {
            index++;
        }
    }
}

// pre: non ci son duplicati tra hypers_to_filter_begin e hypers_to_filter_end
void check_hyperrects_validity_x_categorical_features(vector<HyperRectangle>::iterator hypers_to_filter_begin, vector<HyperRectangle>::iterator hypers_to_filter_end, unordered_set<HyperRectangle> &set_to_fill, vector<vector<size_t>> &categorical_values_columns_index)
{
    while (hypers_to_filter_begin != hypers_to_filter_end)
    {
        if (check_hyperrect_validity_x_categorical_features(*hypers_to_filter_begin, categorical_values_columns_index))
            set_to_fill.insert(*hypers_to_filter_begin);
        hypers_to_filter_begin++;
    }
}

void remove_protected_attr(HyperRectangle &hyper, size_t protected_attribute_index) {
    hyper[protected_attribute_index] = Interval(-INF, INF, 0, 0);
}

void check_validity_and_remove_prot_attribute_interval_hyperrect(size_t protected_attribute_index, HyperRectangle &hyper, unordered_set<HyperRectangle> &set_to_fill, vector<vector<size_t>> &categorical_values_columns_index)
{
    remove_protected_attr(hyper, protected_attribute_index);
    if (check_hyperrect_validity_x_categorical_features(hyper, categorical_values_columns_index))
        set_to_fill.insert(hyper);
}

void check_validity_and_remove_prot_attribute_interval_hyperrects(size_t protected_attribute_index, vector<SymbolicAttack>::iterator hypers_to_convert_begin, vector<SymbolicAttack>::iterator hypers_to_convert_end, unordered_set<HyperRectangle> &set_to_fill, vector<vector<size_t>> &categorical_values_columns_index)
{
    while (hypers_to_convert_begin != hypers_to_convert_end)
    {
        HyperRectangle pre_to_insert = hypers_to_convert_begin->get_pre();
        check_validity_and_remove_prot_attribute_interval_hyperrect(protected_attribute_index, pre_to_insert, set_to_fill, categorical_values_columns_index);
        hypers_to_convert_begin++;
    }
}

size_t compute_complexity_cond(const Condition& cond, vector<size_t>& numerical_features_indexes, vector<vector<size_t>>& categorical_features_indexes) {

    const dynamic_bitset<>& constr_feats = cond.get_constrained_features_ref();
    size_t index = constr_feats.find_first();
    size_t n_constr_feats = 0;
    while(index != dynamic_bitset<>::npos) {
        if(find(numerical_features_indexes.begin(), numerical_features_indexes.end(), index) != numerical_features_indexes.end()){
            n_constr_feats++;
            index = constr_feats.find_next(index);
        }
        else {
            short i = 0;
            while(find(categorical_features_indexes[i].begin(), categorical_features_indexes[i].end(), index) == categorical_features_indexes[i].end())
                i++;
            while(index != dynamic_bitset<>::npos && find(categorical_features_indexes[i].begin(), categorical_features_indexes[i].end(), index) != categorical_features_indexes[i].end())
                index = constr_feats.find_next(index);
            n_constr_feats++;
        }
    }
    return n_constr_feats;
}

void save_hyperrects(string filename, unordered_set<Condition> parity_not_valid_conditions)
{
    std::ofstream output_cond_stream;
    output_cond_stream.open(filename, std::ifstream::out);

    compact_json_stream_encoder encoder(output_cond_stream);

    encoder.begin_array();

    for (Condition c : parity_not_valid_conditions)
    {
        const dynamic_bitset<> &fs = c.get_constrained_features_ref();
        encoder.begin_object();
        encoder.key("dim");
        encoder.uint64_value(c.get_constrained_features_ref().size());
        for (size_t index = fs.find_first(); index != dynamic_bitset<>::npos; index = fs.find_next(index))
        {
            encoder.key(to_string(index));
            float left = (c(index).get_left() == -INF) ? 0 : c(index).get_left();
            float right = (c(index).get_right() == INF) ? 1 : c(index).get_right();
            bool left_type = (c(index).get_left() == -INF) ? 0 : c(index).get_left_type();
            bool right_type = (c(index).get_right() == INF) ? 0 : c(index).get_right_type();
            encoder.begin_array();
            encoder.double_value(left);
            encoder.double_value(right);
            encoder.uint64_value(left_type);
            encoder.uint64_value(right_type);
            encoder.end_array();
        }
        encoder.end_object();
    }

    encoder.end_array();
    encoder.flush();
    output_cond_stream.close();
}

// EXEC BOTH ANALYZER AND SYNTHESISER
int main(int argc, char *argv[])
{
    bool flag_print = false;

    /* check input parameters */
    if (argc < 13)
    {
        std::cout << "./executable <ensemble_filename.json> <columns_filename.json> <categorical_columns_names.json> <numerical_columns_indexes.json> <categorical_columns_indexes.json> <protected_attribute> <n_iter_analyze> <n_threads> <filename_dump_hyperrects.json> <n_iter_synthesiser> <output_conditions_filename_base.json> <normalizations_columns_conditions.json>" << std::endl;
        return -1;
    }

    #if REDUCED_FEATURE_SPACE == 1
        Ensemble *T = load_model(argv[1], true);
        V_st used_features = T->get_f_idx();
    #else
        Ensemble *T = load_model(argv[1]);
    #endif

    string columns_string = read_file(argv[2], "columns");
    vector<string> used_columns;
    auto loaded = jc::decode_json<vector<string>>(columns_string);

    #if REDUCED_FEATURE_SPACE == 1
        for (auto iter = used_features.begin(); iter != used_features.end(); iter++)
            used_columns.push_back(loaded[*iter]);
    #else
        used_columns = loaded;
    #endif

    vector<string> columns = loaded;
    size_t feature_space_dim = columns.size();

    string numerical_columns_indexes_string = read_file(argv[4], "numerical_columns_indexes");
    vector<size_t> numerical_columns_indexes = jc::decode_json<vector<size_t>>(numerical_columns_indexes_string);

    string categorical_columns_indexes_string = read_file(argv[5], "categorical_columns_indexes");
    vector<vector<size_t>> categorical_columns_indexes = jc::decode_json<vector<vector<size_t>>>(categorical_columns_indexes_string);

    map<size_t, size_t> numerical_columns_reduced_map;
    map<size_t, size_t> categorical_columns_reduced_map;

    size_t index = 0;
    while(index < used_features.size()) {
        if(find(numerical_columns_indexes.begin(), numerical_columns_indexes.end(), used_features[index]) != numerical_columns_indexes.end())
            numerical_columns_reduced_map[index] = used_features[index];
        else
            categorical_columns_reduced_map[index] = used_features[index];
        index++;
    }

    string categorical_columns_string = read_file(argv[3], "categorical_columns_names");
    vector<vector<size_t>> categorical_columns_indexes_used;
    vector<vector<string>> categorical_columns_names_used;

    auto loaded2 = jc::decode_json<vector<vector<string>>>(categorical_columns_string);
    for (vector<string> &names : loaded2)
    {
        auto iter = names.begin();
        vector<size_t> categorical_columns_subvector_indexes_used;
        vector<string> categorical_columns_subvector_names_used;
        while (iter != names.end())
        {
            auto found_name_pointer = find(used_columns.begin(), used_columns.end(), *iter);
            if (found_name_pointer != used_columns.end())
            {
                categorical_columns_subvector_indexes_used.push_back(found_name_pointer - used_columns.begin());
                categorical_columns_subvector_names_used.push_back(*found_name_pointer);
            }
            iter++;
        }
        categorical_columns_indexes_used.push_back(categorical_columns_subvector_indexes_used);
        categorical_columns_names_used.push_back(categorical_columns_subvector_names_used);
    }

    string protected_attribute = string(argv[6]);
    auto iter_pa = find(used_columns.begin(), used_columns.end(), protected_attribute);
    size_t protected_attribute_index = (iter_pa != used_columns.end()) ? iter_pa - used_columns.begin() : -1;

    size_t i = 0;
    vector<Interval> atk_v;
    V_f c_atk;
    for (size_t i; i < used_columns.size(); i++)
    {
        if (i != protected_attribute_index)
        {
            atk_v.emplace_back(0, 0, 1, 1);
            c_atk.push_back(10000);
        }
        else
        {
            atk_v.emplace_back(-2, 2, 0, 0); //suppose the feature values to be normalized in 0-1
            c_atk.push_back(1);
        }
    }
    size_t b = 1;
    HyperRectangle atk(atk_v);
    Attacker A(atk, c_atk, b);

    size_t n_iter;
    sscanf(argv[7], "%zu", &n_iter);

    size_t n_threads;
    sscanf(argv[8], "%zu", &n_threads);

    string filename_dump_hyperrects = argv[9];

    set<SymbolicAttack> U;
    set<HyperRectangle> U_preimgs;
    analyze_ensemble(*T, A, n_threads, n_iter, n_iter, U);

    unordered_set<HyperRectangle> parity_not_valid(0);
    unordered_set<HyperRectangle> parity_not_valid_partial(0);
    size_t counter = 0;
    for (auto iter = U.begin(); iter != U.end(); iter++)
    {
        HyperRectangle hyper = iter->get_pre();
        remove_protected_attr(hyper, protected_attribute_index);
        parity_not_valid_partial.insert(hyper);
    }

    for (auto iter = parity_not_valid_partial.begin(); iter != parity_not_valid_partial.end(); iter++)
    {
        HyperRectangle hyper = *iter;
        if (check_hyperrect_validity_x_categorical_features(hyper, categorical_columns_indexes_used))
            parity_not_valid.insert(hyper);
    }

    #if UNION == 1
        single_dim_hyperrectangle_union(parity_not_valid);
    #endif

    // CONVERTING HYPERS TO PARTIAL FUNCTIONS
    unordered_set<Condition> parity_not_valid_conditions(0);

    for (auto iter = parity_not_valid.begin(); iter != parity_not_valid.end(); iter++)
    {   
        Condition condition_complete;
        complete_condition(*iter, condition_complete, feature_space_dim, numerical_columns_reduced_map, categorical_columns_reduced_map, categorical_columns_indexes);
        parity_not_valid_conditions.emplace(condition_complete);
    }

    size_t n_iter_synthesiser;
    sscanf(argv[10], "%zu", &n_iter_synthesiser);

    string output_conditions_filename_base = argv[11];

    list<Condition> parity_holds_conditions;

    // SYNTHESIS
    synthetise(parity_not_valid_conditions, parity_holds_conditions, n_iter_synthesiser, output_conditions_filename_base, numerical_columns_indexes, categorical_columns_indexes);

    #if COMPACT_CONDITIONS == 1
        post_proc_condition_catvar(parity_holds_conditions, numerical_columns_indexes, categorical_columns_indexes);

        vector<list<Condition>> n_conds_x_complexities(n_iter_synthesiser);
        list<Condition> list_printed_conds;
        for(auto iter = parity_holds_conditions.begin(); iter != parity_holds_conditions.end(); iter++) {
            n_conds_x_complexities[iter->get_complexity()-1].push_back(*iter);
        }
        for(unsigned short i = 0; i < n_conds_x_complexities.size(); i++) {
            list_printed_conds.insert(list_printed_conds.end(), n_conds_x_complexities[i].begin(), n_conds_x_complexities[i].end());
            print_conditions(list_printed_conds, output_conditions_filename_base + "_COMPACT_" + to_string(i+1) + "iter.json");
        }
        list_printed_conds.clear();
        n_conds_x_complexities.clear();
    #endif

    #if CONDITIONS_OUTPUT == 1
        string norm_info_string = read_file(argv[12], "normalization");
        vector<vector<float>> normalization_conditions;
        auto loaded1 = jc::decode_json<vector<vector<float>>>(norm_info_string);
        for (vector<float> v : loaded1)
        {
            normalization_conditions.push_back(v);
        }

        cout << "###############OUTPUT CONDITIONS###############\n"
            << endl;
        cout << "N OUTPUT CONDITIONS: " << parity_holds_conditions.size() << "\n"
            << endl;

        sprintf(buffer, "###############OUTPUT CONDITIONS###############\n");
        sprintf(buffer, "N OUTPUT CONDITIONS: %zu\n", parity_holds_conditions.size());
        synthetiser_log_stream.write(buffer, strlen(buffer));
        synthetiser_log_stream.flush();
        for (Condition c : parity_holds_conditions)
        {
            const dynamic_bitset<> &fs = c.get_constrained_features_ref();
            for (size_t index = fs.find_first(); index != dynamic_bitset<>::npos; index = fs.find_next(index))
            {
                // cout << used_columns[iter->first] << ": ";
                float left = c(index).get_left();
                float right = c(index).get_right();
                string left_par = (c(index).get_left_type()) ? "(" : "[";
                string right_par = (c(index).get_right_type()) ? ")" : "]";
                // cout << ((iter->second.get_left_type()) ? "(" : "[") << left << ", " << right << ((iter->second.get_right_type()) ? ")" : "]") << endl;
                sprintf(buffer, "FEATURE %zu: %s%f, %f%s\n", index, left_par.c_str(), left, right, right_par.c_str());
                synthetiser_log_stream.write(buffer, strlen(buffer));
                synthetiser_log_stream.flush();
            }
            sprintf(buffer, "\n");
            synthetiser_log_stream.write(buffer, strlen(buffer));
            synthetiser_log_stream.flush();

            for (size_t index = fs.find_first(); index != dynamic_bitset<>::npos; index = fs.find_next(index))
            {
                // cout << used_columns[iter->first] << ": ";
                float left = c(index).get_left();
                float right = c(index).get_right();
                left = (left != -INF) ? (left * (normalization_conditions[index][1] - normalization_conditions[index][0])) + normalization_conditions[index][0] : normalization_conditions[index][0];
                right = (right != INF) ? (right * (normalization_conditions[index][1] - normalization_conditions[index][0])) + normalization_conditions[index][0] : normalization_conditions[index][1];
                string left_par = (c(index).get_left_type()) ? "(" : "[";
                string right_par = (c(index).get_right_type()) ? ")" : "]";
                // cout << ((iter->second.get_left_type()) ? "(" : "[") << left << ", " << right << ((iter->second.get_right_type()) ? ")" : "]") << endl;
                sprintf(buffer, "%s: %s %f, %f %s\n", columns[index].c_str(), left_par.c_str(), left, right, right_par.c_str());
                synthetiser_log_stream.write(buffer, strlen(buffer));
                synthetiser_log_stream.flush();
            }
            sprintf(buffer, "\n\n");
            synthetiser_log_stream.write(buffer, strlen(buffer));
            synthetiser_log_stream.flush();
        }
    #endif
}
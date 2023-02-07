#include "condition.hpp"
#include <string.h>
#include <unordered_set>
#include <jsoncons/encode_json.hpp>
#include <jsoncons/json.hpp>
#include <iostream>
#include <list>
#define HYPERS_CHECK 1

using namespace jsoncons;

bool is_intersection_empty(const Condition& left, const Condition& right) {
    const dynamic_bitset<>& left_const_feats = left.get_constrained_features_ref();
    const dynamic_bitset<>& right_const_feats = right.get_constrained_features_ref();
    dynamic_bitset<> intersection = left_const_feats & right_const_feats;
    size_t index = intersection.find_first();
    bool intersection_empty = true;
    while(intersection_empty && index != dynamic_bitset<>::npos) {
        if((left(index).intersection(right(index))).is_empty()) {
            intersection_empty = false;
        } 
        index = intersection.find_next(index);
    }
    return intersection_empty;
}

void gen_atom_from_unfair_condition(const Condition& unfair_condition, set<Condition>& conditions, const unordered_set<Condition>& U, uint16_t& base_id, vector<size_t>& numerical_feature_indexes, vector<vector<size_t>>& categorical_features_indexes) {
    const dynamic_bitset<>& constr_feats = unfair_condition.get_constrained_features_ref();
    Interval base_int = Interval(0, 1, 0, 0); 
    size_t index = constr_feats.find_first();
    Interval category_possession = Interval(0, 0.5, 0, 0);
    while(index != dynamic_bitset<>::npos) { 
        vector<Interval> setminus_result = base_int.setminus(unfair_condition(index));
        for(auto iter_res = setminus_result.begin(); iter_res != setminus_result.end(); iter_res++) {
            if(find(numerical_feature_indexes.begin(), numerical_feature_indexes.end(), index) != numerical_feature_indexes.end() || category_possession != *iter_res) {
                Condition condition(constr_feats.size(), U.size());
                condition.insert_atomic_condition(index, *iter_res);
                if(conditions.find(condition) == conditions.end()) {
                    size_t id = 0;
                    for(auto iter_hyper = U.begin(); iter_hyper != U.end(); iter_hyper++) {
                        if(!is_intersection_empty(condition, *iter_hyper)){
                            condition.insert_in_hypers_excluded(id);
                        }
                        id++;
                    }
                    conditions.insert(condition);
                }
            } 
        }
        index = constr_feats.find_next(index);
    }
}

void gen_atoms(const unordered_set<Condition>& syms_pre, set<Condition>& new_conditions, vector<vector<size_t>>& categorical_features_indexes, vector<size_t>& numerical_feature_indexes) {
    uint16_t base_id = 1;
    set<Condition> new_conditions_temp;
    for(auto iter = syms_pre.begin(); iter != syms_pre.end(); iter++) {
        gen_atom_from_unfair_condition(*iter, new_conditions_temp, syms_pre, base_id,  numerical_feature_indexes, categorical_features_indexes);
    }
    for(auto iter = new_conditions_temp.begin(); iter != new_conditions_temp.end(); iter++) {
        Condition c = *iter;
        c.set_id({base_id});
        base_id++;
        new_conditions.insert(c);
    }
}

void post_proc_condition_catvar(list<Condition>& cond_list, vector<size_t>& numerical_features_indexes, vector<vector<size_t>>& categorical_features_indexes) {
    auto prev_size = 0;
    dynamic_bitset<> num_feats_bitmap(cond_list.begin()->get_constr_feats_size());
    for(auto iter_num = numerical_features_indexes.begin(); iter_num != numerical_features_indexes.end(); iter_num++)
        num_feats_bitmap[*iter_num] = 1;

    while(cond_list.size() != prev_size) {

        prev_size = cond_list.size();
        auto iter = cond_list.begin();
        bool union_done = false;

        while (iter != cond_list.end()) {
            Condition union_conditions_catvar = *iter;
            const dynamic_bitset<>& iter_constr_feats = iter->get_constrained_features_ref();
            dynamic_bitset<> feats_num_constrained(iter_constr_feats.size());
            size_t index_mask = iter_constr_feats.find_first();
            while(find(numerical_features_indexes.begin(), numerical_features_indexes.end(), index_mask) != numerical_features_indexes.end()) {
                feats_num_constrained[index_mask] = 1;
                index_mask = iter_constr_feats.find_next(index_mask);
            }
            auto iter2 = iter;
            iter2++;
            
            while(!union_done && iter2 != cond_list.end()) {
                
                const dynamic_bitset<>& iter_constr_feats2 = iter2->get_constrained_features_ref();
                dynamic_bitset<> feats_num_constrained2 = num_feats_bitmap & iter_constr_feats2;
                if(feats_num_constrained2 == feats_num_constrained) { //check su stesse feature numeriche vincolate
                    dynamic_bitset<> different_constr_feats = iter_constr_feats ^ iter_constr_feats2; //ci deve essere almeno una diff nelle variabili cat, altrimenti sarebbe un duplicato
                    size_t diff_index = different_constr_feats.find_first();
                    auto iter_cat = categorical_features_indexes.begin();
                    if(diff_index != dynamic_bitset<>::npos) {
                        while(find(iter_cat->begin(), iter_cat->end(), diff_index) == iter_cat->end())
                            iter_cat++;
                        diff_index = different_constr_feats.find_next(diff_index);
                        while(diff_index != dynamic_bitset<>::npos && find(iter_cat->begin(), iter_cat->end(), diff_index) != iter_cat->end()) {
                            diff_index = different_constr_feats.find_next(diff_index);
                        }
                        auto cat_group_bitmap = dynamic_bitset<>(iter->get_constr_feats_size());
                        for(auto cat_group_iter = iter_cat->begin(); cat_group_iter != iter_cat->end(); cat_group_iter++) 
                            cat_group_bitmap[*cat_group_iter] = 1;
                        if(diff_index == dynamic_bitset<>::npos && different_constr_feats.count() == ((iter_constr_feats & cat_group_bitmap).count() + (iter_constr_feats2 & cat_group_bitmap).count())) { //tutte le differenze in feats vincolate son relative alla stessa var categoriale
                            size_t num_index = feats_num_constrained.find_first();
                            while(num_index != dynamic_bitset<>::npos && iter->operator()(num_index) == iter2->operator()(num_index)) 
                                num_index = feats_num_constrained.find_next(num_index);
                            if(num_index == dynamic_bitset<>::npos){ //i vincoli alle var numeriche sono effettivamente gli stessi
                                auto constr_to_insert = iter_constr_feats2 & cat_group_bitmap;
                                size_t index_constr_to_insert = constr_to_insert.find_first();
                                while(index_constr_to_insert != dynamic_bitset<>::npos){
                                    union_conditions_catvar.insert_atomic_condition(index_constr_to_insert, iter2->operator()(index_constr_to_insert));
                                    index_constr_to_insert = constr_to_insert.find_next(index_constr_to_insert);
                                }
                                union_done = true;
                            }
                        }
                    }
                }
                if(union_done) {
                    cond_list.erase(iter);
                    cond_list.erase(iter2);
                    cond_list.push_back(union_conditions_catvar);
                } else iter2++;
            }
            if(union_done)
                iter = cond_list.end();
            else
                iter++;
        }
        
    }
}

bool filter(const Condition& condition, unordered_set<Condition>& unfair_conditions){
    return condition.get_hypers_excluded_count() == unfair_conditions.size();
}

void gen_and_filter_formulas(map<u16string, list<Condition>>& conditions, map<u16string, list<Condition>>& longer_conditions, list<Condition>& fair_conditions, unordered_set<Condition>& unfair_conditions, vector<size_t>& numerical_values_columns_index, vector<vector<size_t>>& categorical_values_columns_index, bool delete_new_conds, size_t& counter) {
    
    dynamic_bitset<> numerical_feats_mask(conditions.begin()->second.begin()->get_constr_feats_size());
    for(auto iter = numerical_values_columns_index.begin(); iter != numerical_values_columns_index.end(); iter++) 
        numerical_feats_mask[*iter] = 1;
    counter = 0;
    for(auto iter_set_same_prefix = conditions.begin(); iter_set_same_prefix != conditions.end(); iter_set_same_prefix++) {
        auto iter_set_same_prefix_end = iter_set_same_prefix->second.end();
        for(auto iter_cond1 = iter_set_same_prefix->second.begin(); iter_cond1 != iter_set_same_prefix_end; iter_cond1++) {
            const Condition& cond1 = *iter_cond1;
            u16string prefix = cond1.get_id();
            list<Condition> same_prefix_conditions_list;
            auto iter_cond2 = iter_cond1;
            iter_cond2++;
            bool can_subsume = true;
            while(iter_cond2 != iter_set_same_prefix_end) {
                const Condition& cond2 = *iter_cond2;
                if(can_subsume)
                    can_subsume &= subsumes(cond1, cond2, true, numerical_values_columns_index, categorical_values_columns_index, numerical_feats_mask);
                if(!can_subsume) {
                    Condition meet_res;
                    bool defined = meet(cond1, cond2, meet_res, numerical_values_columns_index, categorical_values_columns_index);
                    if(defined) {
                        auto iter_fair_cond = fair_conditions.begin();
                        bool subsumed = false;
                        while(!subsumed && iter_fair_cond != fair_conditions.end()) {
                            subsumed = subsumes(*iter_fair_cond, meet_res, false, numerical_values_columns_index, categorical_values_columns_index, numerical_feats_mask);
                            iter_fair_cond++;
                        }
                        if(!subsumed){
                            if(filter(meet_res, unfair_conditions)) {
                                fair_conditions.push_back(meet_res);
                                post_proc_condition_catvar(fair_conditions, numerical_values_columns_index, categorical_values_columns_index);
                                
                            }
                            else {
                                if(!delete_new_conds){
                                    same_prefix_conditions_list.push_back(meet_res);
                                }
                                else 
                                    counter++;
                            }
                        }
                    }
                }
                iter_cond2++;
            }
            if(same_prefix_conditions_list.size() > 1) {
                longer_conditions[prefix] = std::move(same_prefix_conditions_list);
            }
        }
    }
}

void print_conditions(list<Condition>& cs, string base_name) {
    
    std::ofstream output_cond_stream; 
    output_cond_stream.open(base_name+".json", std::ofstream::out);

    compact_json_stream_encoder encoder(output_cond_stream);

    encoder.begin_array();
    
    for(auto iter=cs.begin(); iter != cs.end(); iter++) {
        Condition c = *iter;
        const dynamic_bitset<>& fs = c.get_constrained_features_ref();
        encoder.begin_object();
        encoder.key("dim");
        encoder.uint64_value(c.get_constrained_features_ref().size());
        for(size_t index = fs.find_first(); index != dynamic_bitset<>::npos; index = fs.find_next(index)){
            encoder.key(to_string(index));
            float left = (c(index).get_left() == -INF) ? 0 : c(index).get_left();
            float right = (c(index).get_right() == INF) ? 1: c(index).get_right();
            bool left_type = (c(index).get_left() == -INF) ? 0: c(index).get_left_type();
            bool right_type = (c(index).get_right() == INF) ? 0: c(index).get_right_type();
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

void print_conditions(unordered_set<Condition>& cs, string base_name) {
    
    std::ofstream output_cond_stream; 
    output_cond_stream.open(base_name+".json", std::ofstream::out);

    compact_json_stream_encoder encoder(output_cond_stream);

    encoder.begin_array();
    
    for(Condition c: cs) {
        const dynamic_bitset<>& fs = c.get_constrained_features_ref();
        encoder.begin_object();
        encoder.key("dim");
        encoder.uint64_value(c.get_constrained_features_ref().size());
        for(size_t index = fs.find_first(); index != dynamic_bitset<>::npos; index = fs.find_next(index)){
            encoder.key(to_string(index));
            float left = (c(index).get_left() == -INF) ? 0 : c(index).get_left();
            float right = (c(index).get_right() == INF) ? 1: c(index).get_right();
            bool left_type = (c(index).get_left() == -INF) ? 0: c(index).get_left_type();
            bool right_type = (c(index).get_right() == INF) ? 0: c(index).get_right_type();
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

//sintetizza le condizioni a partire dagli iperrettangoli contenuti in parity_not_holds
void synthetise(unordered_set<Condition>& U, list<Condition>& F, size_t n_iter, string base_name_cs_output, vector<size_t>& numerical_values_columns_index, vector<vector<size_t>>& categorical_values_columns_index) {

    set<Condition> atomic_conds;
    
    gen_atoms(U, atomic_conds, categorical_values_columns_index, numerical_values_columns_index);

    map<u16string, list<Condition>> Ci_minus_1;

    for(const Condition& c: atomic_conds) {
        if(filter(c, U))
            F.push_back(c);
        else {
            u16string s = c.get_id();
            Ci_minus_1[s.substr(0, s.size()-1)].push_back(c);
        }
    }

    size_t i = 1;
    print_conditions(F, base_name_cs_output + "_" + to_string(i) + "iterS");

    i++;

    while(i-1 != n_iter && Ci_minus_1.size() != 0) {

        map<u16string, list<Condition>> Ci;

        size_t count = 0;

        if(i != n_iter)
            gen_and_filter_formulas(Ci_minus_1, Ci, F, U, numerical_values_columns_index, categorical_values_columns_index, false, count);
        else 
            gen_and_filter_formulas(Ci_minus_1, Ci, F, U, numerical_values_columns_index, categorical_values_columns_index, true, count);

        print_conditions(F, base_name_cs_output + "_" + to_string(i) + "iterS");

        i++;
        Ci_minus_1 = std::move(Ci);
        Ci.clear();
    }
}
#include <jsoncons/encode_json.hpp>
#include <jsoncons/json_cursor.hpp>
#include <jsoncons/json.hpp>
#include <vector>

using namespace std;
using namespace jsoncons;

bool is_instance_inside_hyperrects(vector<float>& instance, unordered_set<HyperRectangle>& hypers) {
    bool contained = false;
    auto iter = hypers.begin();
    while(!contained && iter != hypers.end()) {
        if(instance > (*iter))
            contained = true;
        iter++;
    }
    return contained;
}

size_t count_contained_in(vector<vector<float>>& instances, unordered_set<HyperRectangle>& hypers) {
	size_t counter = 0;
    for(size_t i = 0; i < instances.size(); i++) {
        if(is_instance_inside_hyperrects(instances[i], hypers))
            counter++;
    }
	return counter;
}

/*size_t count_contained_in(vector<vector<float>>& instances, set<Partial_Function>& pfs) {
	size_t counter = 0;
    for(size_t i = 0; i < instances.size(); i++) {
        bool contained = false;
        auto iter = pfs.begin();
        while(!contained && iter != pfs.end()) {
            if(instances[i] > (*iter))
                contained = true;
            iter++;
        }
        if(contained)
            counter++;
    }
	return counter;
}

size_t count_contained_in(vector<vector<float>>& instances, unordered_set<Partial_Function>& pfs) {
	size_t counter = 0;
    for(size_t i = 0; i < instances.size(); i++) {
        bool contained = false;
        auto iter = pfs.begin();
        while(!contained && iter != pfs.end()) {
            if(instances[i] > (*iter))
                contained = true;
            iter++;
        }
        if(contained)
            counter++;
    }
	return counter;
}*/

bool is_instance_inside_hyperrects(vector<float>& instance, vector<Condition>& hypers) {
    bool contained = false;
    auto iter = hypers.begin();
    while(!contained && iter != hypers.end()) {
        if(instance > (*iter))
            contained = true;
        iter++;
    }
    return contained;
}

size_t count_instances_inside_cond(vector<vector<float>>& instances, Condition& pf) {
	size_t counter = 0;
    for(size_t i = 0; i < instances.size(); i++) {
        if(instances[i] > pf)
            counter++;
    }
	return counter;
}

size_t count_contained_in(vector<vector<float>>& instances, vector<Condition>& pfs) {
	size_t counter = 0;
    for(size_t i = 0; i < instances.size(); i++) {
        if(is_instance_inside_hyperrects(instances[i], pfs))
            counter++;
    }
	return counter;
}

bool is_instance_inside_hyperrects(vector<float>& instance, list<Condition>& hypers) {
    bool contained = false;
    auto iter = hypers.begin();
    while(!contained && iter != hypers.end()) {
        if(instance > (*iter))
            contained = true;
        iter++;
    }
    return contained;
}

size_t count_contained_in(vector<vector<float>>& instances, list<Condition>& pfs) {
	size_t counter = 0;
    for(size_t i = 0; i < instances.size(); i++) {
        if(is_instance_inside_hyperrects(instances[i], pfs))
            counter++;
        /*else {
            cout << "INSTANCE NOT INSIDE: " << endl;
            for(short j = 0; j < instances[j].size(); j++) 
                cout << instances[i][j] << ", ";
            cout << endl;
        }*/
    }
	return counter;
}

bool is_instance_inside_hyperrects(vector<float>& instance, unordered_set<Condition>& hypers) {
    bool contained = false;
    auto iter = hypers.begin();
    while(!contained && iter != hypers.end()) {
        if(instance > (*iter))
            contained = true;
        iter++;
    }
    return contained;
}

size_t count_contained_in(vector<vector<float>>& instances, unordered_set<Condition>& pfs) {
	size_t counter = 0;
    for(size_t i = 0; i < instances.size(); i++) {
        if(is_instance_inside_hyperrects(instances[i], pfs))
            counter++;
    }
	return counter;
}

bool is_instance_inside_hyperrects(vector<float>& instance, set<HyperRectangle>& hypers) {
    bool contained = false;
    auto iter = hypers.begin();
    while(!contained && iter != hypers.end()) {
        if(instance > (*iter))
            contained = true;
        iter++;
    }
    return contained;
}

size_t count_contained_in(vector<vector<float>>& instances, set<HyperRectangle>& pfs) {
	size_t counter = 0;
    for(size_t i = 0; i < instances.size(); i++) {
        if(is_instance_inside_hyperrects(instances[i], pfs))
            counter++;
    }
	return counter;
}

bool is_instance_inside_condition(vector<float> instance, Condition condition, vector<size_t> numerical_features_indexes, vector<vector<size_t>> categorical_features_indexes) {
    const dynamic_bitset<> constr_feats = condition.get_constrained_features_ref();
    size_t index = constr_feats.find_first();
    bool contained = true;
    while(contained && index != dynamic_bitset<>::npos){
        if(find(numerical_features_indexes.begin(), numerical_features_indexes.end(), index) != numerical_features_indexes.end()) {
            contained = instance[index] > condition(index);
            index = constr_feats.find_next(index);
        } else {
            short i = 0;
            while(find(categorical_features_indexes[i].begin(), categorical_features_indexes[i].end(), index) == categorical_features_indexes[i].end())
                i++;
            bool included = false;
            while(index != dynamic_bitset<>::npos && find(categorical_features_indexes[i].begin(), categorical_features_indexes[i].end(), index) != categorical_features_indexes[i].end()) {
                if(!included) {
                    included = instance[index] > condition(index);
                }
                index = constr_feats.find_next(index);
            }
            if(!included)
                contained = false;
        }
    }
    return contained;
}

bool is_instance_inside_conditions(vector<float> instance, list<Condition> condition, vector<size_t> numerical_features_indexes, vector<vector<size_t>> categorical_features_indexes) {
    bool is_in_conditions = false;
    auto iter = condition.begin();
    while(!is_in_conditions && iter != condition.end()) {
        is_in_conditions = is_instance_inside_condition(instance, *iter, numerical_features_indexes, categorical_features_indexes);
        /*if(is_in_conditions) {
            cout << "INSTANCE ";
            for(short j = 0; j < instance.size(); j++) 
                cout << instance[j] << ", ";
            cout << endl;
            cout << "IN " << endl;
            iter->print();
            cout << endl;
        }*/
        iter++;
    }
    return is_in_conditions;
}

bool is_instance_inside_conditions(vector<float> instance, vector<Condition> condition, vector<size_t> numerical_features_indexes, vector<vector<size_t>> categorical_features_indexes) {
    bool is_in_conditions = false;
    auto iter = condition.begin();
    while(!is_in_conditions && iter != condition.end()) {
        is_in_conditions = is_instance_inside_condition(instance, *iter, numerical_features_indexes, categorical_features_indexes);
        iter++;
    }
    return is_in_conditions;
}

bool is_instance_inside_conditions(vector<float> instance, unordered_set<Condition> condition, vector<size_t> numerical_features_indexes, vector<vector<size_t>> categorical_features_indexes) {
    bool is_in_conditions = false;
    auto iter = condition.begin();
    while(!is_in_conditions && iter != condition.end()) {
        is_in_conditions = is_instance_inside_condition(instance, *iter, numerical_features_indexes, categorical_features_indexes);
        iter++;
    }
    return is_in_conditions;
}

size_t count_instances_inside_cond(vector<vector<float>>& instances, Condition& pf, vector<size_t> numerical_features_indexes, vector<vector<size_t>> categorical_features_indexes) {
	size_t counter = 0;
    for(size_t i = 0; i < instances.size(); i++) {
        if(is_instance_inside_condition(instances[i], pf, numerical_features_indexes, categorical_features_indexes))
            counter++;
    }
	return counter;
}

size_t count_contained_in_numcatvar(vector<vector<float>>& instances, list<Condition>& pfs, vector<size_t>& numerical_features_indexes, vector<vector<size_t>>& categorical_features_indexes) {
	size_t counter = 0;
    for(size_t i = 0; i < instances.size(); i++) {
        if(is_instance_inside_conditions(instances[i], pfs, numerical_features_indexes, categorical_features_indexes)) {
            counter++;
        }
            
    }
	return counter;
}

size_t count_contained_in_numcatvar(vector<vector<float>>& instances, vector<Condition>& pfs, vector<size_t> numerical_features_indexes, vector<vector<size_t>> categorical_features_indexes) {
	size_t counter = 0;
    for(size_t i = 0; i < instances.size(); i++) {
        if(is_instance_inside_conditions(instances[i], pfs, numerical_features_indexes, categorical_features_indexes))
            counter++;
    }
	return counter;
}

size_t count_contained_in_numcatvar(vector<vector<float>>& instances, unordered_set<Condition>& pfs, vector<size_t> numerical_features_indexes, vector<vector<size_t>> categorical_features_indexes) {
	size_t counter = 0;
    for(size_t i = 0; i < instances.size(); i++) {
        if(is_instance_inside_conditions(instances[i], pfs, numerical_features_indexes, categorical_features_indexes))
            counter++;
    }
	return counter;
}

void parse_hyperrectangles_json(string path_file, unordered_set<Condition>& cs) {
	
	std::ifstream is(path_file);
	json_stream_cursor cursor(is);
	Condition c;
    Condition* pc = NULL;
    size_t dim;
	float left;
	float right;
    bool left_set = false;
	size_t n_left_type;
	size_t n_right_type;
    bool left_type;
	bool right_type;
    bool left_type_set = false;
	size_t feature = 0;
	Interval intv;
    string s;
    std::stringstream sstream;
	while(!(cursor.done())) {
		const auto& event = cursor.current();
        switch (event.event_type())
        {
            case staj_event_type::begin_array:
                break;
            case staj_event_type::end_array:
                if(pc != NULL) {
                    intv = Interval(left, right, left_type, right_type);
                    c.insert_atomic_condition(feature, intv);
                }
                break;
            case staj_event_type::begin_object:
                c = Condition();
                pc = &c;
                break;
            case staj_event_type::end_object:
                cs.insert(c);
                pc = NULL;
                break;
            case staj_event_type::key:
                s = event.get<string>();
                if(s != "dim") {
                    sstream = stringstream(s);
                    sstream >> feature;
                }
                break;
            case staj_event_type::uint64_value:
                if (s == "dim") {
                    dim = event.get<uint64_t>();
                    c.set_n_feats(dim);
                } else {
                    if(left_type_set) {
                        right_type = event.get<uint64_t>();
                        left_type_set = false;
                    } else {
                        left_type = event.get<uint64_t>();
                        left_type_set = true;
                    }
                }
                break;
            case staj_event_type::double_value:
                if(left_set) {
                    right = event.get<double>();
                    left_set = false;
                } else {
                    left = event.get<double>();
                    left_set = true;
                }
                break;
            default:
                std::cout << "Unhandled event type: " << event.event_type() << " " << "\n";
                break;
        }
		cursor.next();
	}

	is.close();
}

void parse_hyperrectangles_json(string path_file, vector<Condition>& cs) {
	
	std::ifstream is(path_file);
	json_stream_cursor cursor(is);
	Condition c;
    Condition* pc = NULL;
    size_t dim;
	float left;
	float right;
    bool left_set = false;
	size_t n_left_type;
	size_t n_right_type;
    bool left_type;
	bool right_type;
    bool left_type_set = false;
	size_t feature = 0;
	Interval intv;
    string s;
    std::stringstream sstream;
	while(!(cursor.done())) {
		const auto& event = cursor.current();
        switch (event.event_type())
        {
            case staj_event_type::begin_array:
                break;
            case staj_event_type::end_array:
                if(pc != NULL) {
                    intv = Interval(left, right, left_type, right_type);
                    c.insert_atomic_condition(feature, intv);
                }
                break;
            case staj_event_type::begin_object:
                c = Condition();
                pc = &c;
                break;
            case staj_event_type::end_object:
                cs.push_back(c);
                pc = NULL;
                break;
            case staj_event_type::key:
                s = event.get<string>();
                if(s != "dim") {
                    sstream = stringstream(s);
                    sstream >> feature;
                }
                break;
            case staj_event_type::uint64_value:
                if (s == "dim") {
                    dim = event.get<uint64_t>();
                    c.set_n_feats(dim);
                } else {
                    if(left_type_set) {
                        right_type = event.get<uint64_t>();
                        left_type_set = false;
                    } else {
                        left_type = event.get<uint64_t>();
                        left_type_set = true;
                    }
                }
                break;
            case staj_event_type::double_value:
                if(left_set) {
                    right = event.get<double>();
                    left_set = false;
                } else {
                    left = event.get<double>();
                    left_set = true;
                }
                break;
            default:
                std::cout << "Unhandled event type: " << event.event_type() << " " << "\n";
                break;
        }
		cursor.next();
	}

	is.close();
}

size_t check_contained_in(VV_f& test_set, vector<Condition>& cs, size_t index) {
    for(size_t i = 0; i != cs.size(); i++) {
        if(test_set[index] > cs[i])
            return i;
    }
    return -1;
}
#include "validation.h"
#include <boost/dynamic_bitset.hpp>

using namespace boost;

class Condition {
    public:
        Condition();

        Condition(short dim);

        Condition(short dim, short n_hypers);

        Condition(HyperRectangle hyp);

		u16string get_id() const;

		void set_id(u16string id);

        void set_n_feats(unsigned short dim);

        void insert_atomic_condition(unsigned short feature, const Interval& interval);

        void insert_atomic_condition(unsigned short feature, const Interval& interval, u16string id);

        void insert_atomic_condition_adding_no_complexity(unsigned short feature, const Interval& interval);

        size_t get_n_constr_feats() const;

        size_t get_n_constr_feats(vector<size_t> numerical_indexes, vector<vector<size_t>> categorical_indexes) const;

        size_t get_constr_feats_size() const;

        size_t get_hypers_excluded_count() const;

        unsigned short get_complexity() const;

        const dynamic_bitset<>& get_constrained_features_ref() const;

        void insert_in_hypers_excluded(unsigned short id);

        const Interval& operator()(size_t f) const;

        bool is_constrained(unsigned short f) const;

        void print() const;

        void print_with_hypers_excluded() const;

        const pair<unsigned short, Interval>& get_last_feature_constr() const;

        //bool has_same_prefix(const Condition& right) const;

        friend bool operator <(const Condition& left, const Condition& right);
        friend bool operator ==(const Condition& left, const Condition& right);
        friend bool operator >(const vector<float>& v, const Condition& right);
        friend bool meet(const Condition& left, const Condition& right, Condition& meet_res, const vector<size_t>& numerical_values_columns_index, const vector<vector<size_t>>& categorical_values_columns_index);
        friend bool subsumes(const Condition& left, const Condition& right, bool same_prefix, vector<size_t>& numerical_features_indexes, vector<vector<size_t>>& categorical_features_indexes, dynamic_bitset<>& numerical_feats_bitmask);
        friend bool compare(const Condition& cond1, const Condition& cond2);

    private:
        dynamic_bitset<> constrained_features;
        dynamic_bitset<> hypers_excluded;
        map<unsigned short, Interval> atomic_conditions;
        u16string id;
        unsigned short hypers_excluded_count;
        unsigned short complexity;
        pair<unsigned short, Interval> last_feature_constr;
};

namespace std
{
    template<> struct hash<Condition>
    {
        size_t operator() (const Condition& p) const {
            size_t hash_output = 0;
            const dynamic_bitset<>& constr_feats = p.get_constrained_features_ref();
            for(size_t index = constr_feats.find_first(); index != dynamic_bitset<>::npos; index = constr_feats.find_next(index)) {
                hash_output ^= index ^ std::hash<Interval>()(p(index));
            }
            return hash_output;
        }
    };
}

Condition::Condition() {
	constrained_features = dynamic_bitset<>(0);
	atomic_conditions = map<unsigned short, Interval>();
	hypers_excluded = dynamic_bitset<>(0);
    complexity = 0;
}

Condition::Condition(short dim) {
	constrained_features = dynamic_bitset<>(dim);
	atomic_conditions = map<unsigned short, Interval>();
	hypers_excluded = dynamic_bitset<>(0);
    complexity = 0;
}

Condition::Condition(short dim, short n_hypers) {
	constrained_features = dynamic_bitset<>(dim);
	hypers_excluded = dynamic_bitset<>(n_hypers);
	atomic_conditions = map<unsigned short, Interval>();
    complexity = 0;
}

Condition::Condition(HyperRectangle hyp) {
	constrained_features = dynamic_bitset<>(hyp.get_dim());
	atomic_conditions = map<unsigned short, Interval>();
	hypers_excluded = dynamic_bitset<>(0);
	hypers_excluded_count = 0;
    complexity = 0;

	const vector<Interval>& intvs_ref = hyp.get_intvs_ref();
	for(size_t i = 0; i < intvs_ref.size(); i++) {
		//Inserisci solo la coppia feature-intervallo se l'intv è differente da (-inf, +inf)
		if(intvs_ref[i].get_left() != -INF || intvs_ref[i].get_right() != INF) {
			constrained_features[i] = 1;
			atomic_conditions[i] = intvs_ref[i];
		}
	}
}

void Condition::insert_atomic_condition(unsigned short feature, const Interval& interval) {
	if(constrained_features[feature]) {
        if(atomic_conditions[feature].get_left() != 0)
            complexity--;
        if(atomic_conditions[feature].get_right() != 1)
            complexity--;
    }
	atomic_conditions[feature] = atomic_conditions[feature].intersection(interval);
	if(interval.get_left() != 0)
		complexity++;
	if(interval.get_right() != 1)
		complexity++;
	constrained_features[feature] = 1;
	last_feature_constr = pair<unsigned short, Interval>(feature, atomic_conditions[feature]); 
}

void Condition::insert_atomic_condition(unsigned short feature, const Interval& interval, u16string id) {
    if(constrained_features[feature]) {
        if(atomic_conditions[feature].get_left() != 0)
            complexity--;
        if(atomic_conditions[feature].get_right() != 1)
            complexity--;
    }
    atomic_conditions[feature] = atomic_conditions[feature].intersection(interval);
    if(atomic_conditions[feature].get_left() != 0)
        complexity++;
    if(atomic_conditions[feature].get_right() != 1)
        complexity++;
	constrained_features[feature] = 1;
	last_feature_constr = pair<unsigned short, Interval>(feature, atomic_conditions[feature]); 
	this->id += id;
}

void Condition::insert_atomic_condition_adding_no_complexity(unsigned short feature, const Interval& interval){
	atomic_conditions[feature] = atomic_conditions[feature].intersection(interval);
	constrained_features[feature] = 1;
	last_feature_constr = pair<unsigned short, Interval>(feature, atomic_conditions[feature]); 
}

u16string Condition::get_id() const {
	return id;
}

void Condition::set_id(u16string id) {
	this->id = id;
}

void Condition::set_n_feats(unsigned short dim) {
	constrained_features.resize(dim);
}

size_t Condition::get_n_constr_feats() const {
	return constrained_features.count();
}

size_t Condition::get_n_constr_feats(vector<size_t> numerical_indexes, vector<vector<size_t>> categorical_indexes) const {
    size_t counter = 0;
    size_t index = constrained_features.find_first();
    while(index != dynamic_bitset<>::npos) {
        if(find(numerical_indexes.begin(), numerical_indexes.end(), index) != numerical_indexes.end()) {
            counter++;
            index = constrained_features.find_next(index);
        } else {
            size_t i = 0;
            bool found = false;
            while(!found && i < categorical_indexes.size()) {
                if(find(categorical_indexes[i].begin(), categorical_indexes[i].end(), index) != categorical_indexes[i].end())
                    found = true;
                else i++;
            }
            size_t sub_counter=1;
            bool found_1 = false;
            index = constrained_features.find_next(index);
            while(index != dynamic_bitset<>::npos && find(categorical_indexes[i].begin(), categorical_indexes[i].end(), index) != categorical_indexes[i].end()){
                sub_counter++;
                if(atomic_conditions.at(index).get_right() == 1)
                    found_1 = true;
                index = constrained_features.find_next(index);
            }
            if(found_1)
                counter++;
            else counter += sub_counter;
        }
    }
	return counter;
}

size_t Condition::get_constr_feats_size() const {
	return constrained_features.size();
}

size_t Condition::get_hypers_excluded_count() const {
	return hypers_excluded.count();
}

const dynamic_bitset<>& Condition::get_constrained_features_ref() const {
	return constrained_features;
}

unsigned short Condition::get_complexity() const {
    return complexity;
}

void Condition::insert_in_hypers_excluded(unsigned short id) {
	hypers_excluded[id] = 1;
	hypers_excluded_count++;
}

const Interval& Condition::operator()(size_t f) const{
	return atomic_conditions.at(f);
}

bool Condition::is_constrained(unsigned short f) const {
	return constrained_features[f];
}

void Condition::print() const {
	for(size_t index = constrained_features.find_first(); index != dynamic_bitset<>::npos; index = constrained_features.find_next(index)) {
		cout << "FEATURE " << index << ": ";
		atomic_conditions.at(index).print();
	}
}

void Condition::print_with_hypers_excluded() const {
	print();

    size_t index = hypers_excluded.find_first();

    cout << "HYPERS EXCLUDED: [";
    while(index != dynamic_bitset<>::npos) {
        cout << index << ", ";
        index = hypers_excluded.find_next(index);
    }
    cout << " ]" << endl;
}

/*bool Condition::has_same_prefix(const Condition& right, const vector<size_t>& numerical_features_indexes, const vector<vector<size_t>>& categorical_features_indexes) const {
    
}*/

const pair<unsigned short, Interval>& Condition::get_last_feature_constr() const {
	return last_feature_constr;
}

bool compare(const Condition& cond1, const Condition& cond2) {
    size_t n_constr_feats1 = cond1.constrained_features.count();
    size_t n_constr_feats2 = cond2.constrained_features.count();

    if(n_constr_feats1 == n_constr_feats2) {
        if(cond1.constrained_features == cond2.constrained_features) {
            size_t index = cond1.constrained_features.find_first();

            while(index != dynamic_bitset<>::npos) {
                if(cond1(index).get_left() == 0 && cond2(index).get_left() == 0 && cond1(index).get_right() != cond2(index).get_right())
                    return cond1(index).get_right() > cond2(index).get_right();
                else if(cond1(index).get_right() == 0 && cond2(index).get_right() == 0 && cond1(index).get_left() != cond2(index).get_left())
                    return cond1(index).get_left() < cond2(index).get_left();
                else if (cond1(index).get_left() == cond2(index).get_left() && cond1(index).get_right() != cond2(index).get_right())  
                    return cond1(index).get_right() > cond2(index).get_right();
                else if(cond1(index).get_right() == cond2(index).get_right() && cond1(index).get_left() != cond2(index).get_left())
                    return cond1(index).get_left() < cond2(index).get_left();
                else if(cond1(index).get_left() != cond2(index).get_left() && cond1(index).get_right() != cond2(index).get_right()) {
                    if(cond1(index).get_left() == 0)
                        return true;
                    else if (cond2(index).get_left() == 0)
                        return false;
                    else return cond1(index).get_right() > cond2(index).get_right();
                }  
                else index = cond1.constrained_features.find_next(index);
            }
            return false;
        } else return cond1.constrained_features < cond2.constrained_features;
    } else return n_constr_feats1 < n_constr_feats2;
}

//funziona solo per le condizioni atomiche (ed è utilizzato solo lì) TODO: convertire
bool operator <(const Condition& left, const Condition& right){
    /*if(left.constrained_features == right.constrained_features) {
        size_t index = left.constrained_features.find_first();

        while(index != dynamic_bitset<>::npos) {
            if(left(index).get_left() == 0 && right(index).get_left() == 0)
                return left(index).get_right() > right(index).get_right();
            else return left(index).get_left() < right(index).get_left();
        }
        return false;

    } else return left.constrained_features < right.constrained_features;*/
    return compare(left, right);
}

bool operator ==(const Condition& left, const Condition& right) {
    return left.complexity == right.complexity && left.constrained_features == right.constrained_features && left.atomic_conditions == right.atomic_conditions;
}

bool operator >(const vector<float>& v, const Condition& right) {
    if(v.size() != right.constrained_features.size())
        return false;
    
    bool contained = true;
    size_t index = right.constrained_features.find_first();

    while(contained && index != dynamic_bitset<>::npos) {
        if(!(v[index] > right(index)))
            contained = false;
        index = right.constrained_features.find_next(index);
    }
    return contained;
}

bool subsumes(const Condition& left, const Condition& right, bool same_prefix, vector<size_t>& numerical_features_indexes, vector<vector<size_t>>& categorical_features_indexes, dynamic_bitset<>& numerical_feats_bitmask) {

    if(same_prefix) {
        const pair<unsigned short, Interval>& left_last_feature_constr = left.get_last_feature_constr();
        const pair<unsigned short, Interval>& right_last_feature_constr = right.get_last_feature_constr();
        if(left_last_feature_constr.first == right_last_feature_constr.first) {
            if(left_last_feature_constr.second.intersection(right_last_feature_constr.second) == right_last_feature_constr.second)
                return true;
        }
        return false;
    } else { //sfrutto il pre-proc dei dataset per cui gli indici di var cat sono dopo gli indici delle var numeriche
        const dynamic_bitset<>& constr_feats_left = left.get_constrained_features_ref();
        bool subsumes = true;
        size_t index = constr_feats_left.find_first();
        size_t last_numerical_features_index = numerical_features_indexes[numerical_features_indexes.size()-1];
        while(subsumes && index != dynamic_bitset<>::npos && index <= last_numerical_features_index){
            subsumes = right.is_constrained(index) && ((right.atomic_conditions).at(index).intersection(left.atomic_conditions.at(index)) == right.atomic_conditions.at(index));
            index = constr_feats_left.find_next(index);
        }   
        auto iter = categorical_features_indexes.begin();
        while(subsumes && index != dynamic_bitset<>::npos) {
            while(iter != categorical_features_indexes.end() && find(iter->begin(), iter->end(), index) == iter->end())
                iter++; 
            size_t last_cat_feats_index_iter = iter->operator[](iter->size()-1);
            bool same_cat_feats_constr = false;
            while(index != dynamic_bitset<>::npos && index <= last_cat_feats_index_iter) {
                if(!same_cat_feats_constr)
                    same_cat_feats_constr = right.is_constrained(index);
                index = constr_feats_left.find_next(index);
            }
            subsumes = same_cat_feats_constr;
        }         
        return subsumes;
    }
}

//si suppone che left e right abbiano lo stesso prefisso e che left non subsumi right
bool meet(const Condition& left, const Condition& right, Condition& meet_res, const vector<size_t>& numerical_feature_indexes, const vector<vector<size_t>>& categorical_feature_indexes) {
    
    const pair<unsigned short, Interval>& right_last_feature_constr = right.get_last_feature_constr();
    const pair<unsigned short, Interval>& left_last_feature_constr = left.get_last_feature_constr();
    if(left_last_feature_constr.first == right_last_feature_constr.first && left_last_feature_constr.second.intersection(right_last_feature_constr.second).is_empty())
        return false;
    #if HYPERS_CHECK==1
        dynamic_bitset<> meet_res_hypers_excluded = left.hypers_excluded | right.hypers_excluded;
        if(meet_res_hypers_excluded.count() == left.hypers_excluded.count() || meet_res_hypers_excluded.count() == right.hypers_excluded.count())
            return false;
    #endif

    bool defined = true; 

    if(find(numerical_feature_indexes.begin(), numerical_feature_indexes.end(), right_last_feature_constr.first) != numerical_feature_indexes.end()) {
        if(left.is_constrained(right_last_feature_constr.first))
            defined = !(right_last_feature_constr.second.intersection(left.atomic_conditions.at(right_last_feature_constr.first)).is_empty());
    } else {
        unsigned short i = 0;
        bool found = false;
        while(!found && i < categorical_feature_indexes.size())
            if(find(categorical_feature_indexes[i].begin(), categorical_feature_indexes[i].end(), right_last_feature_constr.first) != categorical_feature_indexes[i].end())
                found = true;
            else i++;
        unsigned short j = 0;
        //cout << "INDEX: " << right_last_feature_constr.first << ", GROUP: " << i << endl;
        while(defined && j < categorical_feature_indexes[i].size()) {
            //cout << "INDEX " << categorical_feature_indexes[i][j] << " OF SAME GROUP IS CONSTRAINED?" << endl;
            if(left.is_constrained(categorical_feature_indexes[i][j])) {
                //cout << "YES, THE LEFT EXTREME OF THE INTV IS: " << left(categorical_feature_indexes[i][j]).get_left() << endl;
                defined = !(left(categorical_feature_indexes[i][j]).get_left() == 0.5);
            }
            j++;
        }
    }
    if(defined) {
        meet_res = left;
        meet_res.insert_atomic_condition(right_last_feature_constr.first, right_last_feature_constr.second, right.id.substr(right.id.length()-1, 1));
        #if HYPERS_CHECK==1
            meet_res.hypers_excluded = meet_res_hypers_excluded;
        #else
            meet_res.hypers_excluded = left.hypers_excluded;
            meet_res.hypers_excluded |= right.hypers_excluded;
        #endif
            meet_res.hypers_excluded_count = meet_res.hypers_excluded.count();
        //meet_res.print();
    }
    return defined;
}
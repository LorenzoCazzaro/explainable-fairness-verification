#include <cstdlib>
#include <iterator>
#include <queue>
#include <algorithm>

Ensemble ens;

struct priority_pairs_cut_ind {
  bool operator() (const SymbolicAttack& a, const SymbolicAttack& b) {
    bool less = a < b;
    if (!less && !(b < a)) //for removing sym atk duplicates
        return a < b;
    else {
        size_t a_split_counter = a.get_split_counter(); 
        size_t b_split_counter = b.get_split_counter();
        if(a_split_counter != b_split_counter) 
            return a_split_counter < b_split_counter;
        else {
            std::vector<Tree*> uT1;
            validate_s2(a, ens, &uT1);
            std::vector<Tree*> uT2;
            validate_s2(b, ens, &uT2);
            size_t uTsizea = uT1.size();
            size_t uTsizeb = uT2.size();
            if (uTsizea != uTsizeb)
                return uTsizea < uTsizeb;
            else return less;
        }
    }
  }
};


/* analyze tree (same as in the article) */
//TODO: eventualmente ritornare a std::set<SymbolicAttack>& U
void analyze_tree(std::set<SymbolicAttack, priority_pairs_cut_ind>& U, Tree &t, Attacker& A) {
    t.annotate(A);
    std::vector<Leaf*> leaves;
    t.get_leaves(leaves);

    float ext_lb, int_lb;

    for (Leaf* ext_lf : leaves) {
        ext_lb = ext_lf->get_label();
        for (const SymbolicAttack& ext_s : *ext_lf->get_sym()) {
            if (ext_s.get_paid() == 0) {
                for (Leaf* int_lf : leaves) {
                    int_lb = int_lf->get_label();
                    if (ext_lb != int_lb) {
                        for (const SymbolicAttack& int_s : *int_lf->get_sym()) {
                            if (int_s.get_paid() > 0) {
                                HyperRectangle inter_s_pre = ext_s.get_pre().intersection(int_s.get_pre());
                                if (!inter_s_pre.is_empty()) {
                                    HyperRectangle inter_s_post = int_s.get_post().intersection(inter_s_pre + A.get_atk()); // OLD
                                    //HyperRectangle inter_s_post = int_s.get_post().intersection(int_s.get_atk_f() + inter_s_pre);
                                    U.emplace(SymbolicAttack(inter_s_pre, inter_s_post, int_s.get_paid(), int_s.get_atk_f()));
                                }
                            }
                        }
                    }
                }
            }  
        }
    }
}

void analyze_tree(std::set<SymbolicAttack>& U, Tree &t, Attacker& A) {
    t.annotate(A);
    std::vector<Leaf*> leaves;
    t.get_leaves(leaves);

    float ext_lb, int_lb;

    for (Leaf* ext_lf : leaves) {
        ext_lb = ext_lf->get_label();
        for (const SymbolicAttack& ext_s : *ext_lf->get_sym()) {
            if (ext_s.get_paid() == 0) {
                for (Leaf* int_lf : leaves) {
                    int_lb = int_lf->get_label();
                    if (ext_lb != int_lb) {
                        for (const SymbolicAttack& int_s : *int_lf->get_sym()) {
                            if (int_s.get_paid() > 0) {
                                HyperRectangle inter_s_pre = ext_s.get_pre().intersection(int_s.get_pre());
                                if (!inter_s_pre.is_empty()) {
                                    //HyperRectangle inter_s_post = int_s.get_post().intersection(inter_s_pre + A.get_atk()); // OLD
                                    HyperRectangle inter_s_post = int_s.get_post().intersection(int_s.get_atk_f() + inter_s_pre);
                                    U.emplace(SymbolicAttack(inter_s_pre, inter_s_post, int_s.get_paid(), int_s.get_atk_f()));
                                }
                            }
                        }
                    }
                }
            }  
        }
    }
}

/* mode indica se per l'SA 'u' è stato usato tutto il budget disponibile o no */
void sa_from_intvs(const SymbolicAttack& u, std::set<SymbolicAttack>& U, Attacker& A, std::vector<Interval>& intvs, size_t f, bool mode) {
    const std::vector<Interval>& u_pre = u.get_pre().get_intvs_ref();
    const std::vector<Interval>& u_post = u.get_post().get_intvs_ref();

    SymbolicAttack s_u;
    for (Interval& i : intvs) {
        s_u = u;
        s_u.get_pre_noconst()[f] = i;

        if (mode) {
            if (u_pre[f] == u_post[f])
                s_u.get_post_noconst()[f] = i;
            else
                //s_u.get_post_noconst()[f] = u_post[f].intersection(u.get_atk_f()[f] + i);
                s_u.get_post_noconst()[f] = u_post[f].intersection(i + A[f]); // OLD
        } else
            s_u.get_post_noconst()[f] = u_post[f].intersection(i + A[f]);

        s_u.split_increment();

        if (!s_u.is_empty())
            U.insert(s_u);
    }
}

void split_interval_in_4(const Interval& i, std::vector<Interval>& intvs, float v, const Interval atk) {
    float atk_left = atk.get_left();
    float atk_right = atk.get_right();

    if (v > i) {
        if (v + atk_left > i) {
            intvs.emplace_back(i.get_left(), v + atk_left, i.get_left_type(), 0);
            intvs.emplace_back(v + atk_left, v, 1, 0);
        } else
            intvs.emplace_back(i.get_left(), v, i.get_left_type(), 0);

        if (v + atk_right > i) {
            intvs.emplace_back(v, v + atk_right, 1, 0);
            intvs.emplace_back(v + atk_right, i.get_right(), 1, i.get_right_type());    
        } else
            intvs.emplace_back(v, i.get_right(), 1, i.get_right_type());
    } else {
        cout << "SBRANZO" << endl; //lasciamolo intanto, è un alert simpatico
        if (v + atk_left > i) {
            intvs.emplace_back(v + atk_left, i.get_right(), 0, i.get_right_type());
            intvs.emplace_back(i.get_left(), v + atk_left, i.get_left_type(), 1);
        }

        if (v + atk_right > i) {
            intvs.emplace_back(i.get_left(), v + atk_right, i.get_left_type(), 0);
            intvs.emplace_back(v + atk_right, i.get_right(), 1, i.get_right_type());    
        }
    }
}

/* MOLTO PIÙ VELOCE U più piccolo */
/* restituisce la threshold con l'occorrenza più alta tra gli laberi indecisi */
VM_fst get_ths_occ(std::vector<Tree*>& T, size_t dim) {
    VM_fst ths_occ = VM_fst(dim);
   
    for (Tree* t : T) {
        VS_f ths_p = VS_f(dim);
        t->fill_ths(ths_p);

        for (size_t f = 0; f < dim; ++f) {
            for (float v : ths_p[f]) {
                if (ths_occ[f].find(v) == ths_occ[f].end())
                    ths_occ[f][v] = 1;
                else
                    ths_occ[f][v]++;
            }   
        }
    }

    return ths_occ;
}

/* MOLTO PIÙ VELOCE U più piccolo */
/* restituisce la threshold con l'occorrenza più alta tra gli laberi indecisi */
VM_fst get_leaf_ths_occ(std::vector<Tree*>& T, size_t dim) {
    VM_fst leaf_ths_occ = VM_fst(dim);
            
    for (Tree* t : T) {
        VS_f ths_p = VS_f(dim);
        t->fill_leaf_ths(ths_p);

        for (size_t f = 0; f < dim; ++f) {
            for (float v : ths_p[f]) {
                if (leaf_ths_occ[f].find(v) == leaf_ths_occ[f].end())
                    leaf_ths_occ[f][v] = 1;
                else
                    leaf_ths_occ[f][v]++;
            }   
        }
    }

    return leaf_ths_occ;
}

/* seleziona la threshold con l'occorrenza più alta tra gli laberi indecisi */
bool split_max_occ(const SymbolicAttack& u, std::vector<Tree*>& T, Attacker& A, size_t& out_f, float& out_v, double p = 0.7) {
    size_t dim = A.get_dim();
    const std::vector<Interval>& u_pre = u.get_pre().get_intvs_ref();

    /* cerca le threshold degli alberi indecisi */
    VM_fst ths = get_ths_occ(T, dim);

    size_t f_max = 0;
    float v_max = 0;
    size_t max_occ = 0;
    bool flag = false;
    for (size_t f = 0; f < dim; ++f) {
        Interval I = Interval(u_pre[f].get_left(), u_pre[f].get_right(), 1, 1);
        for (P_fst p : ths[f]) {
            if (p.second > max_occ) {
                float v = p.first;
                /* qui senza  v > I || più veloce */
                if ( v > I  ) {//|| (v + A[f].get_left()) > I || (v + A[f].get_right()) > I ) {
                    f_max = f;
                    v_max = p.first;
                    max_occ = p.second;
                    flag = true;
                }
            }
        }
    }

    out_f = f_max;
    out_v = v_max;
    return flag;
}

/* seleziona la threshold con l'occorrenza più alta tra gli laberi indecisi */
bool split_max_leaf_occ(const SymbolicAttack& u, std::vector<Tree*>& T, Attacker& A, size_t& out_f, float& out_v, double p = 0.7) {
    size_t dim = A.get_dim();
    const std::vector<Interval>& u_pre = u.get_pre().get_intvs_ref();

    /* cerca le threshold degli alberi indecisi */
    VM_fst ths = get_leaf_ths_occ(T, dim);

    size_t f_max = 0;
    float v_max = 0;
    size_t max_occ = 0;
    bool flag = false;
    for (size_t f = 0; f < dim; ++f) {
        Interval I = Interval(u_pre[f].get_left(), u_pre[f].get_right(), 1, 1);
        for (P_fst p : ths[f]) {
            if (p.second > max_occ) {
                float v = p.first;
                /* qui senza  v > I || più veloce */
		if ((v + A[f].get_left()) > I || (v + A[f].get_right()) > I) {
			cout << "SBRANZO2" << endl;
		}
                if ( v > I  || (v + A[f].get_left()) > I || (v + A[f].get_right()) > I ) {
                    f_max = f;
                    v_max = p.first;
                    max_occ = p.second;
                    flag = true;
                }
            }
        }
    }

    out_f = f_max;
    out_v = v_max;
    return flag;
}

/* seleziona la feature con l'occorrenza più alta tra gli laberi indecisi */
bool split_max_leaf_occ_f(const SymbolicAttack& u, std::vector<Tree*>& T, Attacker& A, size_t& out_f, float& out_v) {
    size_t dim = A.get_dim();
    const std::vector<Interval>& u_pre = u.get_pre().get_intvs_ref();

    std::vector<int> f_occ(dim, 0);
    VM_fst f_ths = VM_fst(dim);

    bool flag = false;
    for (Tree* t : T) {
        bool tree_flag = true;
        VS_f ths_p = VS_f(dim);
        t->fill_leaf_ths(ths_p);

        for (size_t f = 0; f < dim; ++f) {
            Interval I = Interval(u_pre[f].get_left(), u_pre[f].get_right(), 1, 1);

            for (float v : ths_p[f]) {
                if ( v > I || (v + A[f].get_left()) > I || (v + A[f].get_right()) > I ) {
                    /* conta sola un'occorrenza della feature per albero */
                    if (tree_flag) {
                        f_occ[f]++;
                        tree_flag = false;
                    }

                    if (f_ths[f].find(v) == f_ths[f].end())
                        f_ths[f][v] = 1;
                    else
                        f_ths[f][v]++;

                    flag = true;
                }
            }                
        }
    }

    size_t max = 0;
    size_t arg_max = 0;
    size_t i = 0;

    for (size_t occ : f_occ) {
        if (occ > max) {
            max = occ;
            arg_max = i;
        }
        ++i;
    }

    out_f = arg_max;
    if (f_ths[arg_max].size()==0)
        return false;

    max_map(f_ths[arg_max], &out_v);
    return true;
}



/* seleziona la threshold con l'occorrenza più alta tra gli laberi indecisi */
bool split_min_occ(const SymbolicAttack& u, std::vector<Tree*>& T, Attacker& A, size_t& out_f, float& out_v, double p = 0.7) {
    size_t dim = A.get_dim();
    const std::vector<Interval>& u_pre = u.get_pre().get_intvs_ref();

    /* cerca le threshold degli alberi indecisi */
    VM_fst ths = get_ths_occ(T, dim);

    size_t f_min = 0;
    float v_min = 0;
    size_t min_occ = (size_t)-1;
    bool flag = false;
    for (size_t f = 0; f < dim; ++f) {
        Interval I = Interval(u_pre[f].get_left(), u_pre[f].get_right(), 1, 1);
        for (P_fst p : ths[f]) {
            if (p.second < min_occ) {
                float v = p.first;
                /* qui senza  v > I || più veloce */
                if ( v > I || (v + A[f].get_left()) > I || (v + A[f].get_right()) > I ) {
                    f_min = f;
                    v_min = p.first;
                    min_occ = p.second;
                    flag = true;
                }
            }
        }
    }

    out_f = f_min;
    out_v = v_min;
    return flag;
}

/* funzione utile se cambiamo split (da rimuovere a lavoro finito) */
void split_sa(const SymbolicAttack& u, std::set<SymbolicAttack>& U, std::vector<Tree*>& T, Attacker& A) {
    size_t f;
    float v;
    bool res = split_max_occ(u, T, A, f, v, 0);
    //bool res = split_max_leaf_occ_f(u, T, A, f, v);
    // bool res = split_min_occ(u, T, A, f, v, 0);
    if (res) {
        const std::vector<Interval>& u_pre = u.get_pre().get_intvs_ref();
        std::vector<Interval> intvs;
        split_interval_in_4(u_pre[f], intvs, v, A[f]);
        sa_from_intvs(u, U, A, intvs, f, u.get_paid() == A.get_budget());
    }
}
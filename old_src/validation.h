#include "certification.h"

#include <vector>
#include <set>
#include <set>
#include <map>
#include <algorithm>
#include <random>

#include <cstddef>


/* verifica che la maggioranze della foresta concorda su s*/
bool validate_s(const SymbolicAttack& s, Ensemble& T, std::vector<Tree*>* uT=std::nullptr_t(), bool early_stopping=false) {
    std::map<float, size_t> counter;
    Counter count_pre, count_post;
    size_t count_err = 0;
    size_t e_size = T.size();

    for (float l : T.get_labels())
        counter[l] = 0;

    for (Tree* t : T.get_trees()) {
        bool robust = false;
        count_post = (*t)(s.get_post(), T.get_n_labels());
        if (count_post.size() == 1) {
            count_pre = (*t)(s.get_pre(), T.get_n_labels());
            /* TODO rimettere count_pre == count_post invece di confrontare i primi valori (FORSE PIÙ EFFICIENTE COSÌ) */
            auto first_a = count_pre.get_set().begin();
            auto first_b = count_post.get_set().begin();
            robust = (count_pre.size() == 1 && (*first_a) == (*first_b));
        }
        if (robust) {
            auto first = count_pre.get_set().begin();
            counter[*first]++;
        } else {
            count_err++;
            if (uT!=std::nullptr_t())
                uT->push_back(t);
        }

        // TODO se non ci interessano gli alberi indecisi allora va bene terminare prima
        if (early_stopping && count_err > int(e_size / 2))
            return false;
    }
    return (max_map(counter) > int(e_size / 2));
}

/* verifica che la maggioranze della foresta concorda su s*/
/* ritorna 0 se safe, altrimenti un codice del tipo */
size_t validate_s23(const SymbolicAttack& s, Ensemble& T, std::vector<Tree*>* uT=std::nullptr_t(), float* l=nullptr) {
    float pre_label = 0;
    float post_label = 0;

    bool pre_check = T.predict_ensemble(s.get_pre(), pre_label);
    bool post_check = T.predict_ensemble(s.get_post(), post_label);

    size_t type = 0;

    /* pre e post hanno una sola label ed è la stessa (SAFE) */
    if (pre_check==true && post_check==true && pre_label == post_label)
        return type;

    /* TODO nel caso di multi labels (> 2) occorre controllare se una delle pre è quella della post ANCHE QUI */
    /* pre e post hanno una sola label ma diversa. Questo SA non può essere raffinato ulteriormente (ESAUSTO) */
    if (pre_check==true && post_check==true && pre_label != post_label)
        type = 1;

    /* pre con più label, post con una sola TODO nel caso di multi labels (> 2) occorre controllare se una delle pre è quella della post */
    if (pre_check==false && post_check==true) {
        type = 2;
        if (l!=nullptr)
            *l = post_label;
    }

    /* pre con una label, post con più label */
    if (pre_check==true && post_check==false)
        type = 3;

    /* pre e post più label */
    if (pre_check==false && post_check==false)
        type = 4;
    
    /* cerca gli alberi indecisi TODO soluzine implementativamente veloce ma non efficinete */
    validate_s(s, T, uT);
    return type;
}


/* verifica che la maggioranze della foresta concorda su u*/
bool validate_u(const SymbolicAttack& u, Ensemble& T, std::vector<Tree*>* uT=std::nullptr_t(), bool early_stopping=false) {
    std::map<float, size_t> counter;
    Counter count_pre, count_post;
    size_t count_err = 0;
    size_t e_size = T.size();

    /* initializes a counter map with labels as key and 0 as vlaue. */
    for (float l : T.get_labels())
        counter[l] = 0;

    /* for each tree inside the ensemble computes its agreement wrt pre and post images. */
    for (Tree* t : T.get_trees()) {
        bool robust = false;
        count_post = (*t)(u.get_post(), T.get_n_labels());

        if (count_post.size() == 1) {
            count_pre = (*t)(u.get_pre(), T.get_n_labels());
            robust = (count_pre.size() == 1 && (count_pre == count_post));
        }

        /* if the tree is robust wrt u, it increments the map counter. */
        if (robust) {
            auto first = count_pre.get_set().begin();
            counter[*first]++;

            /* if the majority of the forest agree on the label return true (robust). */
            bool res = (max_map(counter) > int(e_size / 2));
            if (res)
                return true;
            
        } else {
            count_err++;
            /* adds undecided tree */
            if (uT!=std::nullptr_t())
                uT->push_back(t);
        }

        /* if the majority of the forest are undecided return false (unsafe). */
        if (early_stopping && count_err > int(e_size / 2))
            return false;
    }

    /* if the majority of the forest agree on the label return true (robust), false otherwise. */
    return (max_map(counter) > int(e_size / 2));
}

/* checks type of a SA `u` passed as argument. Returns throught `uT` a vector of undecided tree of `u`. */
size_t validate_s2(const SymbolicAttack& u, Ensemble& T, std::vector<Tree*>* uT=std::nullptr_t(), float* l=nullptr) {
    float pre_label = 0;
    float post_label = 0;

    bool pre_check = T.predict_ensemble(u.get_pre(), pre_label);
    bool post_check = T.predict_ensemble(u.get_post(), post_label);
    size_t type = 3;

    /* pre e post hanno una sola label ed è la stessa (SAFE) */
    if (pre_check==true && post_check==true)
        /* TODO nel caso di multi labels (> 2) occorre controllare se una delle pre è quella della post ANCHE QUI */
        return (pre_label == post_label) ? 0 : 1; /* pre e post hanno una sola label ma diversa. Questo SA non può essere raffinato ulteriormente (ESAUSTO) */

    /* pre con più label, post con una sola TODO nel caso di multi labels (> 2) occorre controllare se una delle pre è quella della post */
    if (pre_check==false && post_check==true) {
        if (l!=nullptr)
            *l = post_label;
        type = 2;
    }
    
    /* TODO soluzine implementativamente veloce ma non efficinete */
    /* finds undecided */
    validate_u(u, T, uT);
    return type;
}

/* TRUE ROBUSTNESS */
bool atk_rec(Ensemble& T, Attacker& A, VS_f& ths, V_f& x, float y, size_t f, size_t b) {
    if (b <= 0 || f >= A.get_dim()) {
        return (T.predict(x) == y);
    }
    
    bool out = atk_rec(T, A, ths, x, y, f + 1, b);
    if (out == false)
        return false;

    float orig = x[f];
    float c = A.get_cost()[f];
    /* bisogna aprire l'estremo perché se orig + atk = v allora non deve saltare */
    Interval intv(A[f].get_left(), A[f].get_right(), 0, 1);
    for (float v : ths[f]) {
        if (v > (orig + intv)) {
            float val = (orig <= v)? std::nextafterf(v, INF) : v;
            x[f] = val;
            out = atk_rec(T, A, ths, x, y, f + 1, b - c);
            x[f] = orig;
            if (!out)
                return false;
        }
    }

    x[f] = orig;
    return true;
}

float atk(Ensemble& T, Attacker& A, VS_f& ths, VV_f& X, V_f& Y, size_t b) {
    float count = 0;
    for (size_t i = 0; i < Y.size(); ++i) {
        if(T.predict(X[i]) == Y[i] && atk_rec(T, A, ths, X[i], Y[i], 0, b))
            count++;
    }
    return count / X.size();
}

float atk(Ensemble& T, Attacker& A, VS_f& ths, VV_f& X, V_f& Y, size_t b, vector<size_t>& rob_instances_index) {
    float count = 0;
    for (size_t i = 0; i < Y.size(); ++i) {
        if(T.predict(X[i]) == Y[i] && atk_rec(T, A, ths, X[i], Y[i], 0, b)) {
            count++;
            rob_instances_index.push_back(i);
        }
    }
    return count / X.size();
}

float true_robustness(Ensemble& T, Attacker& A, VV_f& X, V_f& Y, size_t b) {
    VS_f ths = T.get_ths();
    return atk(T, A, ths, X, Y, b);
}

float true_robustness(Ensemble& T, Attacker& A, VV_f& X, V_f& Y, size_t b, vector<size_t>& rob_instances_index) {
    VS_f ths = T.get_ths();
    return atk(T, A, ths, X, Y, b, rob_instances_index);
}











/* verifica la correttezza della predict ensemble. Per ogni SA campiona randomicamente un'istanza nella pre e una nella posta e verifica che la predizione sia la stessa.
Se cos' non fosse c'è un essore o nella predizione dell'ensemble o nella verifica di safety dell'SA. */
void validate_S(std::set<SymbolicAttack>& S, Ensemble& T, size_t iter=1000) {
    for (const SymbolicAttack& s : S) {
        const std::vector<Interval>& lpre = s.get_pre().get_intvs_ref();
        const std::vector<Interval>& lpost = s.get_post().get_intvs_ref();
        size_t dim = s.get_pre().get_dim();

        for (int j = 0; j < iter; ++j) {
            V_f ist_pre;
            V_f ist_post;

            /* per ogni dimensione... */
            for (int i = 0; i < dim; ++i) {
                /* campiono dalla pre image un valore dall'intervallo della dimensione i */
                /* TODO SISTEMARE TUTTI GLI NEXTAFTER PERCHE SONO ROTTI (mettere + 0.0...1 minimo float perchè after prec è buggato o buggato codice. Da indagare )*/
                float left_pre = (lpre[i].get_left_type()==OPEN) ? std::nextafterf(lpre[i].get_left(), INF) + 0.0001 : lpre[i].get_left();
                float right_pre = (lpre[i].get_right_type()==OPEN) ? std::nextafterf(lpre[i].get_right(), -INF) - 0.0001 : lpre[i].get_right();
                float val_pre = left_pre + static_cast <float> (rand()) / ( static_cast <float> (RAND_MAX / (right_pre - left_pre)));
                /* assegno il valore all'isanza finale che identifica la pre */
                ist_pre.push_back(val_pre);

                /* campiono dalla post image un valore dall'intervallo della dimensione i */
                float left_post = (lpost[i].get_left_type()==OPEN) ? std::nextafterf(lpost[i].get_left(), INF) + 0.0001 : lpost[i].get_left();
                float right_post = (lpost[i].get_right_type()==OPEN) ? std::nextafterf(lpost[i].get_right(), -INF) - 0.0001 : lpost[i].get_right();
                float val_post = left_post + static_cast <float> (rand()) / ( static_cast <float> (RAND_MAX / (right_post - left_post)));
                /* assegno il valore all'isanza finale che identifica la post */
                ist_post.push_back(val_post);
            }

            /* eseguo la predizone */
            float pred_pre = T.predict(ist_pre);
            float pred_post = T.predict(ist_post);
            /* se le due predizioni non coincidono c'+ un problema */
            if (pred_pre != pred_post) {
                /* stampe di debug dell'errore */
                std::cout << pred_pre << " " << pred_post << std::endl;
                s.print();
                
                float pre_label, post_label;
                bool pre_check = T.score_ensemble(s.get_pre(), pre_label);
                bool post_check = T.score_ensemble(s.get_post(), post_label);

                std::cout << pre_check << " " << post_check << " " << pre_label << " " << post_label << std::endl;

                print_ist(ist_pre);
                print_ist(ist_post);
                break;
            }
        }
    }
}

/* verifica che tutto cio che è fuori u pre non sia mai attaccabile */
/* in breve, cerca di campionare un istanza x che sta al di fuori delle pre image di U e poi verifica se attaccando (in modo randomico) è possibile generare un attacco.
Se trova un attacco vorrebbe dire che l'insimeme U non copre tutti i possibili attacchi. Dato che x non appartiene alle pre image di U vuol dire che c'è un errore. */
bool validate_U(std::set<SymbolicAttack>& U, Ensemble& T, Attacker& A, size_t iter=1000) {
    size_t dim = A.get_dim();
    VS_f ths = T.get_ths();
    size_t ist_creates = 0;

    /* vettore delle features, da 0 a dim - 1 */
    std::vector<size_t> rand_pos;
    for (size_t i = 0; i < dim; ++i)
        rand_pos.push_back(i);

    /* crea al massimo iter istanze */
    for (size_t i = 0; i < iter; ++i) {
        bool robust = false;
        size_t try_b = 0;
        V_f ist_org;
        V_f ist_atk(dim);

        /* campionare al di fuori delle pre image di U non è semplice.
        L'algoritmo ci prova per almeno iter volte, se non riesci rinuncia */
        while(!robust && try_b < 100) {
            V_f ist_org_app(dim);

            /* per ogni feature seleziona a caso un valore.
            Dato che ci sono intervalli a -inf, inf campioniamo le threshold degli alberi.
            Questo perché se dobbiamo poi generare un attacco non a senso farlo se siamo troppo distanti dalla soglia (non sarebbe mai possibile generare un attacco efficacie).
            */
            for (size_t f = 0; f < ths.size(); ++f) {
                S_f v_list = ths[f];
                if (v_list.size() != 0) {
                    /* campiona la soglia */
                    size_t pos =  static_cast <size_t> (rand()) / ( static_cast <size_t> (RAND_MAX / v_list.size()));
                    float v = *std::next(v_list.begin(), pos);
                    /* a destra o a sinistra della soglia */
                    size_t p = static_cast <size_t> (rand()) / ( static_cast <size_t> (RAND_MAX / 2));
                    //v = (p) ? std::nextafterf(v, INF) : std::nextafterf(v, -INF);
                    v = (p) ? v + 0.001 : v - 0.001;
                    /* genera l'istanza originale */
                    ist_org_app[f] = v;
                } else
                    /* se la feature non è usata nel ensemble mette 0 */
                    ist_org_app[f] = 0;
            }

            /* verifica se sta dentro una pre image di U. Se non è sentro a nessuna pre image allora è robusta */
            robust = !check_in(ist_org_app, U);
            if (robust)
                ist_org = std::move(ist_org_app);
            try_b++;
        }

        if (robust) {
            ist_creates++;
            /* ordina le feature in modo randomico per generare gli attacchi */
            std::shuffle(rand_pos.begin(), rand_pos.end(), std::default_random_engine(rand()));
            size_t b = A.get_budget();

            /* seleziona una feature */
            for (size_t f : rand_pos) {
                /* se l'attaccante ha ancora sufficiente budget */
                if (ths[f].size() && b - A.get_cost()[f] >= 0) {
                    /* paga il costo dell'atttacco */
                    b = b - A.get_cost()[f];
                    /* campiona randomicamente l'intensità dell'attacco */      
                    float atk_left = (A[f].get_left_type()==OPEN) ? std::nextafterf(A[f].get_left(), INF) : A[f].get_left();
                    float atk_right = (A[f].get_right_type()==OPEN) ? std::nextafterf(A[f].get_right(), -INF) : A[f].get_right();
                    float atk_delta = atk_left + static_cast <float> (rand()) / ( static_cast <float> (RAND_MAX / (atk_right - atk_left)));
                    /* genera l'istanza attaccata partendo da quella originale */
                    ist_atk[f] = ist_org[f] + atk_delta;
                } else
                    /* se la feature non è usata nel ensemble lascia l'originale */
                    ist_atk[f] = ist_org[f];
            }

            /* esegue le predizioni */
            float pred_pre = T.predict(ist_org);
            float pred_post = T.predict(ist_atk);
            /* se le predizioni sono diverse vuol dire che è possibile generare un attacco partendo da un istanza al di fuori di U, il che vìola la teoria del paper */
            if (pred_pre != pred_post) {
                if (robust) {
                    std::cout << "ERRORE" << std::endl;
                    std::cout << pred_pre << " " << pred_post << std::endl;
                    print_ist(ist_org);
                    print_ist(ist_atk);
                    return false;
                }
            }
        }
    }

    std::cout << "\tchecked istances " << ist_creates << std::endl;
    return true;
}
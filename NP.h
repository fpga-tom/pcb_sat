//
// Created by tomas on 1/1/19.
//

#ifndef PCB_SAT_NP_H
#define PCB_SAT_NP_H

#include <vector>
#include <random>
#include <map>
#include <cryptominisat5/solvertypesmini.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>

class NP {
    const std::vector<std::vector<CMSat::Lit>>& clauses;

    typedef Eigen::SparseMatrix<int> mat_t;
    typedef Eigen::Array<bool, Eigen::Dynamic, 1> vec_t;
    typedef Eigen::Triplet<int> T;

    mat_t kuramoto;
    uint32_t terminal;

    std::map<uint32_t , uint32_t > lit_terminal_map;
    std::map<std::pair<uint32_t , uint8_t >, std::pair<uint32_t, bool> > clause_terminal_map;
    std::map<uint32_t , CMSat::Lit > terminal_lit_map;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution;

public:
    NP(const std::vector<std::vector<CMSat::Lit>>& clauses) : clauses(clauses), terminal(0), distribution(0.0,1.0) {



        std::map<uint32_t , std::vector<uint32_t> > var_terminal;


        for (uint32_t c = 0; c < clauses.size(); c++) {
            const std::vector<CMSat::Lit> &clause = clauses[c];
            if (clause.size() == 1) {
                std::abort();
            } else if (clause.size() == 2) {
                std::abort();
            } else if (clause.size() == 3) {
                for (uint64_t i = 0; i < clause.size(); i++) {
                    const CMSat::Lit &lit = clause[i];
                    clause_terminal_map.insert(std::make_pair(std::make_pair(c, i), std::make_pair(terminal, lit.sign())));
                    terminal_lit_map.insert(std::make_pair(terminal, lit));
                    lit_terminal_map.insert(std::make_pair(lit.toInt(), terminal));
                    if(var_terminal.count(lit.var()) == 0){
                        var_terminal.insert(std::make_pair(lit.var(), std::vector<uint32_t >()));
                    }
                    var_terminal[lit.var()].emplace_back(terminal);
                    terminal++;
                }
            }
        }

        kuramoto = mat_t(terminal, terminal);
        /*
        for(int i = 0 ; i < terminal; i++) {
            for(int j = 0; j < terminal; j++) {
                kuramoto.insert(i,j) = 0;
            }
        }
         */

        std::vector<T> triplets;
        for (uint32_t c = 0; c < clauses.size(); c++) {
            const std::vector<CMSat::Lit> &clause = clauses[c];
            if (clause.size() == 1) {
                std::abort();
            } else if (clause.size() == 2) {
                std::abort();
            } else if (clause.size() == 3) {
                for (uint64_t i = 0; i < clause.size(); i++) {
                    auto& a = clause_terminal_map[std::make_pair(c, i)];
                    auto& b = clause_terminal_map[std::make_pair(c, (i + 1) % 3)];

                    CMSat::Lit& lit = terminal_lit_map[b.first];
                    CMSat::Lit& lit1 = terminal_lit_map[a.first];

                    if(lit.sign() || lit1.sign())
//                        kuramoto.insert(a.first, b.first) = -1;
                        triplets.emplace_back(T(a.first, b.first, -1));
                    else
//                        kuramoto.insert(a.first, b.first) = 1;
                        triplets.emplace_back(T(a.first, b.first, 1));


                    if(lit1.sign() || lit.sign())
//                        kuramoto.insert(b.first, a.first) = -1;
                        triplets.emplace_back(T(b.first, a.first, -1));
                    else
//                        kuramoto.insert(b.first, a.first) = 1;
                        triplets.emplace_back(T(b.first, a.first,1));

                    std::vector<uint32_t > l = var_terminal[terminal_lit_map[a.first].var()];
                    for(int q = 0; q < l.size(); q++) {
                        int t = l[q];
                        CMSat::Lit& lit = terminal_lit_map[t];
                        if(a.second == false && lit.sign() == true) {
                            triplets.emplace_back(T(a.first, t,-1));
                            triplets.emplace_back(T(t, a.first,-1));
//                            kuramoto.insert(a.first, t) = -1;
//                            kuramoto.insert(t, a.first) = -1;
                        }
                    }
                }
            }
        }
        kuramoto.setFromTriplets(triplets.begin(), triplets.end());

    }

    void solve() {

        std::string label[] = {"a", "b", "c", "~b", "~c", "~d", "~a", "c", "d", "a", "~b", "~d"};
        vec_t s(terminal);
        for(int i = 0; i < terminal; i++) {

            if(terminal_lit_map[i].sign())
                s[i] = false;
            else
                s[i] = true;
        }
        for(int i = 0; i < 158125000; i++) {
            s = mul(s, kuramoto);
//            auto r = (s[0] or s[1] or s[2]) and (not s[3] or not s[4] or not s[5]) and (not s[6] or s[7] or s[8]) and (s[9] or not s[10] or not s[11]);
            auto r = check(s);
//            for(int j = 0; j < terminal;j++) {
//                std::cout << "[" << j << "," << s[j] << "] ";
//            }
            std::cout << i << std::endl;



//            auto _a = s[0];
//            auto _b = s[1];
//            auto _c = s[2];
//            auto _d = s[8];
//            r = (_a or _b or _c) and (not _b or not _c or not _d) and (not _a or _c or _d) and (_a or not _b or not _d);
//            r = (s[0] or s[1] or s[2]) and ( s[3] or  s[4] or  s[5]) and ( s[6] or s[7] or s[8]) and (s[9] or s[10] or  s[11]);
//            std::cout << s[0] << "," << s[1] << "," << s[2] << "," << s[8] << "," << r << std::endl;
            if(r) {

                for(int j = 0; j < terminal;j++) {
                std::cout << "[" << j << "," << s[j] << "] ";
            }
            std::cout << r << std::endl;
                std::cout << "got it " << i << std::endl;
                break;
            }
        }



    }

    bool check(const vec_t& s) {
        bool ret = true;
        for(int c = 0; c < clauses.size(); c++) {
            const std::vector<CMSat::Lit>& clause = clauses[c];
            bool cnf = false;
            for(int i = 0; i < clause.size(); i++) {
                const CMSat::Lit& l = clause[i];
                auto t = clause_terminal_map[std::make_pair(c, i)];

                if(l.sign())
                    cnf |= s[t.first];
                else
                    cnf |= s[t.first];

#if 1
                for(int c1 = c; c1 < clauses.size(); c1++) {
                    const std::vector<CMSat::Lit>& clause1 = clauses[c1];
                    for(int i1 = 0; i1 < clause1.size(); i1++) {
                        const CMSat::Lit& l1 = clause1[i1];
                        if(l.var() == l1.var()) {
                            auto t1 = clause_terminal_map[std::make_pair(c1, i1)];
                            if(l.sign() != l1.sign()) {
                                if(s[t.first] == s[t1.first])
                                    return false;
                            }
                            if(l.sign() == l1.sign()) {
                                if(s[t.first] != s[t1.first])
                                    return false;
                            }
                        }
                    }
                }
#endif


            }
            ret &= cnf;
            if(ret == false)
                return false;
        }
        return ret;
    }

    vec_t mul(const vec_t& s, const mat_t& m) {
        vec_t ret(terminal);
        for(int i = 0; i < m.outerSize(); i++) {

            ret[i] = s[i];

            int maj = 0;
            for(mat_t::InnerIterator j(m, i); j ; ++j) {
                if(j.value() == 1) {
                    if(distribution(generator) < .97)
//                        ret[i] |= s[i] ^ s[j];
                        maj += (s[j.col()] == s[j.row()] ? 1 : -1);
//                        ret[i] ^= not (s[j]);
                } else if(j.value() == -1) {
                    if(distribution(generator) < .97)
//                        ret[i] |= s[i] ^ (!s[j]);
                        maj += (s[j.col()] != s[j.row()] ? 1 : -1);
//                        ret[i] ^= (s[j]);
                }
            }
            if(maj < 0) {
                ret[i] = not s[i];
            }
        }

        return ret;
    }



};


#endif //PCB_SAT_NP_H

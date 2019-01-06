#include "Dynamic.h"
#include <iostream>
#include <cryptominisat5/cryptominisat.h>
#include <assert.h>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include "Expression.h"
#include "Solg_solver.h"
#include "Circuit.h"
#include <boost/algorithm/string.hpp>
#include "NP.h"

using std::vector;
using namespace CMSat;

using namespace expression::boolean;

typedef std::map<uint64_t , std::map<uint64_t , expr_t>> mat_t;

lbool find(uint64_t num_pages);

int main(int argc, char** argv) {
    for (uint64_t np = 2; np < 45; np++) {
        lbool ret = find(np);
        std::cout << np << " " << ret << std::endl;
        if (ret == l_True) {
            break;
        }
    }
}


lbool find(uint64_t num_pages) {
#if 1
//    uint64_t num_pages = 10;
#if 0
    std::vector<uint64_t> V = {1,2,3,4,5};
    std::vector<std::pair<uint64_t ,uint64_t >> E = {
            std::make_pair(1,2),
            std::make_pair(2,3),
            std::make_pair(3,4),
            std::make_pair(4,5),
            std::make_pair(5,1),
            std::make_pair(1,3),
            std::make_pair(1,4),
            std::make_pair(2,4),
            std::make_pair(2,5),
            std::make_pair(3,5)
    };


#else
    std::vector<uint64_t > V;
    std::vector<std::pair<int ,int >> E;

    srand(3);
    int V_count = 9;
    for(int i = 0; i < V_count; i++) {
        V.emplace_back(i);
    }
    int pairs = 9;
    int max_pins = 5;
    for(int i = 0; i < pairs; i++) {
        auto r1 = rand() % V_count;
        auto pins = rand() % max_pins + 2;
        for( int j = 0; j < pins; ++j) {
            auto r2 = rand() % V_count;
            while (r1 == r2) {
                r2 = rand() % V_count;
            }
            auto p = std::make_pair(r1, r2);
            auto p1 = std::make_pair(r2, r1);
            auto it = std::find(E.begin(), E.end(), p);
            if (it == E.end()) {
                auto it1 = std::find(E.begin(), E.end(), p1);
                if (it1 == E.end()) {
                    if(E.size() < pairs)
                        E.emplace_back(p);
                }
            }
        }
        if(E.size() > pairs) {
            break;
        }
    }

    std::cout << E.size() << std::endl;

#endif

    std::cout << "generating tree..." << std::endl;
    mat_t left_of;


    std::cout << "generating directions..." << std::endl;
    expr_t direction;
#if 0
    for(uint64_t i = 0; i < V.size(); ++i) {
        for(uint64_t j = 0; j < V.size(); j++) {
//            if (i == j) continue;
            if(left_of.count(V[i]) == 0) {
                left_of.insert(std::pair<uint64_t , std::map<uint64_t ,expr_t>>(V[i], std::map<uint64_t, expr_t>()));
            }
            auto &var = left_of[V[i]];

            if (V[i] > V[j]) {
                var.insert(std::pair<uint64_t , expr_t>(V[j], !left_of[V[j]][V[i]]) );
            } else {
                expr_t v(new expr_s("dir_" + std::to_string(V[i]) + "_" + std::to_string(V[j])));
                if( i == j )
                    var.insert(std::pair<uint64_t , expr_t>(V[j], !v));
                else
                    var.insert(std::pair<uint64_t , expr_t>(V[j], v));
                if(direction.get() == nullptr)
                    direction = v;
                else
                    direction = direction && v;
            }
        }
    }
#else

    for(uint64_t i = 0; i < V.size(); ++i) {
        for(uint64_t j = i; j < V.size(); j++) {
            if (i == j) continue;
            if(left_of.count(V[i]) == 0) {
                left_of.insert(std::pair<uint64_t , std::map<uint64_t ,expr_t>>(V[i], std::map<uint64_t, expr_t>()));
            }
            if(left_of.count(V[j]) == 0) {
                left_of.insert(std::pair<uint64_t , std::map<uint64_t ,expr_t>>(V[j], std::map<uint64_t, expr_t>()));
            }
            auto &var = left_of[V[i]];

            expr_t v(new expr_s("dir_" + std::to_string(V[i]) + "_" + std::to_string(V[j])));
            if(i != j) {
                var.insert(std::pair<uint64_t, expr_t>(V[j], v));
                left_of[V[j]].insert(std::pair<uint64_t, expr_t>(V[i], !left_of[V[i]][V[j]]));
                if (direction.get() == nullptr)
                    direction = left_of[V[i]][V[j]] || left_of[V[j]][V[i]];
                else
                    direction = direction && (left_of[V[i]][V[j]] || left_of[V[j]][V[i]]);
            } else {
                var.insert(std::pair<uint64_t, expr_t>(V[j], !v));
                if (direction.get() == nullptr)
                    direction = left_of[V[i]][V[j]];
                else
                    direction = direction && left_of[V[i]][V[j]];
            }
        }
    }
#endif

    std::cout << "generating pages..." << std::endl;
    mat_t page_assign;

    expr_t all_pages;
    for(int e = 0; e < E.size(); ++e) {
        expr_t edge;
        for (int p = 0; p < num_pages; ++p) {
            if(page_assign.count(p) == 0) {
                page_assign.insert(std::pair<uint64_t, std::map<uint64_t, expr_t>>(p, std::map<uint64_t, expr_t>()));
            }
            auto &var = page_assign[p];

            expr_t v(new expr_s("page_" + std::to_string(p) + "_" + std::to_string(e)));
            var.insert(std::pair<uint64_t , expr_t>(e, v));
            if (edge.get() == nullptr)
                edge = v;
            else
                edge = edge || v;

        }
        if( all_pages.get() == nullptr )
            all_pages = edge;
        else
            all_pages = all_pages && edge;
    }


    std::cout << "generating same..." << std::endl;
    mat_t same_page_rule;

    for (int i = 0;i < E.size(); ++i) {
        for (int j = i+1; j < E.size(); ++j) {
            if (i == j) continue;
            expr_t rule;
            for(int p = 0; p < num_pages; ++p) {
                if ( rule.get() == nullptr )
                    rule = page_assign[p][i] && page_assign[p][j];
                else
                    rule = (rule || (page_assign[p][i] && page_assign[p][j]));
            }

            if(same_page_rule.count(i) == 0) {
                same_page_rule.insert(std::pair<uint64_t , std::map<uint64_t , expr_t>>(i, std::map<uint64_t, expr_t>()));
            }
            same_page_rule[i].insert(std::pair<uint64_t, expr_t>(j, rule));
        }
    }

    std::cout << "generating cross..." << std::endl;
    expr_t cross;

    int count = 0;
    uint64_t delta_sum = 0;
    for(int i = 0 ; i < E.size(); ++i) {
        for(int j = i+1; j < E.size(); ++j) {
            if (i == j) continue;
            uint64_t vi = E[i].first;
            uint64_t vj = E[i].second;
            uint64_t vk = E[j].first;
            uint64_t vl = E[j].second;

            if(vi == vk || vi == vl || vj == vk  ||  vj == vl) continue;

            expr_t c = (
                    (same_page_rule[i][j]) >>
                    (
                            !(left_of[vi][vk] && left_of[vk][vj] && left_of[vj][vl]) && !(left_of[vi][vl] && left_of[vl][vj] && left_of[vj][vk]) &&
                            !(left_of[vj][vk] && left_of[vk][vi] && left_of[vi][vl]) && !(left_of[vj][vl] && left_of[vl][vi] && left_of[vi][vk]) &&
                            !(left_of[vk][vi] && left_of[vi][vl] && left_of[vl][vj]) && !(left_of[vk][vj] && left_of[vj][vl] && left_of[vl][vi]) &&
                            !(left_of[vl][vi] && left_of[vi][vk] && left_of[vk][vj]) && !(left_of[vl][vj] && left_of[vj][vk] && left_of[vk][vi])
                    ));

            if(cross.get() == nullptr)
                cross = c;
            else {
                cross = cross && c;
            }

        }
    }

    expr_t target = all_pages && direction && cross;


    std::cout << "generating solver.." << std::endl;
    SATSolver solver;
    solver.set_num_threads(1);


    std::vector<uint64_t > linear;
    std::cout << "linearizing..." << std::endl;
    linearize lz(linear);
    lz(target);
    std::cout << "tseitin..." << std::endl;
    solver.new_vars(linear.size());
    std::vector<std::vector<Lit>> clauses;
    tseitin tseitin(clauses, linear);
    tseitin(target);

    std::cout << "vars " << linear.size() << std::endl;
    std::cout << "clauses " << clauses.size() << std::endl;


    auto v1 = linear.size();
    auto v2 = v1 + 1;
    solver.new_var();
    solver.new_var();

    clauses.emplace_back(std::vector<Lit>{Lit(tseitin.indexof(target), false), Lit(v1, false), Lit(v2, false)});
    clauses.emplace_back(std::vector<Lit>{Lit(tseitin.indexof(target), false), Lit(v1, false), Lit(v2, true)});
    clauses.emplace_back(std::vector<Lit>{Lit(tseitin.indexof(target), false), Lit(v1, true), Lit(v2, false)});
    clauses.emplace_back(std::vector<Lit>{Lit(tseitin.indexof(target), false), Lit(v1, true), Lit(v2, true)});

    for(auto&c : clauses) {
        if(c.size() == 1) {
            solver.new_var();
            solver.new_var();
        }
        if(c.size() == 2) {
            solver.new_var();
        }
    }

    std::vector<std::vector<Lit>> clauses_to_add;
    int new_lit = 0;
    for(auto&c : clauses) {
        if(c.size() == 1) {
            c.emplace_back(Lit(linear.size() + new_lit, false));
            c.emplace_back(Lit(linear.size() + new_lit+1, false));

            clauses_to_add.emplace_back(std::vector<Lit>{Lit(c[0]), Lit(new_lit, false), Lit(new_lit+1, true)});
            clauses_to_add.emplace_back(std::vector<Lit>{Lit(c[0]), Lit(new_lit, true), Lit(new_lit+1, false)});
            clauses_to_add.emplace_back(std::vector<Lit>{Lit(c[0]), Lit(new_lit, true), Lit(new_lit+1, true)});
            new_lit+=2;
        }
        if(c.size() == 2) {
            c.emplace_back(Lit(linear.size() + new_lit, false));
            clauses_to_add.emplace_back(std::vector<Lit>{Lit(c[0]), Lit(c[1]), Lit(new_lit, true)});
            new_lit++;
        }
    }

    for(const auto& c : clauses) {
        solver.add_clause(c);
    }


    solver.set_verbosity(0);
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    lbool ret = solver.solve();
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    uint64_t delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
    std::cout << "Cryptominisat done " << (delta_us/1e6) << std::endl;

#if 1
    if(ret == l_True) {
        std::vector<lbool> model = solver.get_model();

//        for(int i = 0 ; i < E.size(); ++i) {
//            for (int j = i + 1; j < E.size(); ++j) {
//                if (i == j) continue;
//                uint64_t vi = E[i].first;
//                uint64_t vj = E[i].second;
//                uint64_t vk = E[j].first;
//                uint64_t vl = E[j].second;
//
//                if (vi == vk || vi == vl || vj == vk || vj == vl) continue;
//
//                uint64_t r = tseitin.indexof(same_page_rule[i][j]);
//                std::cout << model[r] << " " << r << " " << i << " " << j << std::endl;
//            }
//        }

        uint64_t r = tseitin.indexof(all_pages);
        std::cout << model[r] << " " << r  << std::endl;

//        for(int p = 0; p < num_pages; p++) {
//            for(int e = 0; e < E.size(); e++) {
//                uint64_t r = tseitin.indexof(page_assign[p][e]);
//                if(model[r] == l_True) {
//                    std::cout << model[r] << " " << p << " " << e << std::endl;
//                }
//            }
//        }

//        for(int i = 0; i < V.size(); ++i) {
//            for(int j = 0; j < V.size(); ++j) {
//                if(i == j) continue;
//                uint64_t a = V[i];
//                uint64_t b = V[j];
//                uint64_t r = tseitin.indexof(left_of[a][b]);
//                if(model[r] == l_True) {
//                    std::cout << model[r] << " left_of " << a << " " << b << std::endl;
//                }
//            }
//        }

        std::vector<uint64_t > s(V.size());
        std:iota(s.begin(), s.end(), 0);
        std::sort(s.begin(), s.end(), [&tseitin, &model,&left_of](uint64_t a, uint64_t b) -> bool {
            uint64_t r = tseitin.indexof(left_of[a][b]);
            return model[r] == l_True;
        });

//        std::copy(s.begin(), s.end(), std::ostream_iterator<uint64_t >(std::cout));
        for(const auto& a : s) {
            std::cout << a << ",";
        }

        std::cout << std::endl;

    }
#endif
#endif

//    Solg_solver solg(Solg_solver::gate_t::_or);
//    solg.rk4();
#if 0
    std::cout << "Creating circuit" << std::endl;
    solg::circuit::Circuit circ(linear, clauses);
//    std::cout << "Circuit nodes " << circ.nodes.size() << std::endl;
//    std::cout << "Circuit edges " << circ.edges.size() << std::endl;
    std::cout << "Solving circuit" << std::endl;


    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    circ.solve();
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
    std::cout << "Solving done " << (delta_us/1e6) << std::endl;

#endif

#if 1
    std::vector<std::vector<Lit>> clauses1;

    std::vector<Lit> clause;
    clause.emplace_back(Lit(0, false));
    clause.emplace_back(Lit(1, false));
    clause.emplace_back(Lit(2, false));
    clauses1.emplace_back(clause);
    clause.clear();

    clause.emplace_back(Lit(1, false));
    clause.emplace_back(Lit(2, true));
    clause.emplace_back(Lit(3, true));
    clauses1.emplace_back(clause);
    clause.clear();

    clause.emplace_back(Lit(0, true));
    clause.emplace_back(Lit(2, false));
    clause.emplace_back(Lit(3, false));
    clauses1.emplace_back(clause);
    clause.clear();

    clause.emplace_back(Lit(0, false));
    clause.emplace_back(Lit(1, true));
    clause.emplace_back(Lit(3, true));
    clauses1.emplace_back(clause);
    clause.clear();

    clause.emplace_back(Lit(0, false));
    clause.emplace_back(Lit(1, true));
    clause.emplace_back(Lit(4, true));
    clauses1.emplace_back(clause);
    clause.clear();

    clause.emplace_back(Lit(0, false));
    clause.emplace_back(Lit(2, true));
    clause.emplace_back(Lit(4, true));
    clauses1.emplace_back(clause);
    clause.clear();

    clause.emplace_back(Lit(1, false));
    clause.emplace_back(Lit(3, false));
    clause.emplace_back(Lit(4, false));
    clauses1.emplace_back(clause);
    clause.clear();

    clause.emplace_back(Lit(0, false));
    clause.emplace_back(Lit(3, true));
    clause.emplace_back(Lit(4, false));
    clauses1.emplace_back(clause);
    clause.clear();

    clause.emplace_back(Lit(5, false));
    clause.emplace_back(Lit(3, true));
    clause.emplace_back(Lit(4, false));
    clauses1.emplace_back(clause);
    clause.clear();

    clause.emplace_back(Lit(5, false));
    clause.emplace_back(Lit(2, true));
    clause.emplace_back(Lit(1, false));
    clauses1.emplace_back(clause);
    clause.clear();

    clause.emplace_back(Lit(5, false));
    clause.emplace_back(Lit(2, false));
    clause.emplace_back(Lit(3, true));
    clauses1.emplace_back(clause);
    clause.clear();

    NP np(clauses);
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    np.solve();
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
    std::cout << "Solving done " << (delta_us/1e6) << std::endl;

//    Dynamic d(7, 1, 1050);
//    for(int i = 0; i < 25000; i++)
//        d.step(i);

#endif
    return ret;
}

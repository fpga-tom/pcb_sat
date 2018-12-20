#include <iostream>
#include <cryptominisat5/cryptominisat.h>
#include <boost/algorithm/string.hpp>
#include <assert.h>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include "Expression.h"

using std::vector;
using namespace CMSat;

typedef std::map<uint64_t , std::map<uint64_t , expr_t>> mat_t;

lbool find(uint64_t num_pages);

int main(int argc, char** argv) {
    for (uint64_t np = 2; np < 15; np++) {
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

    srand(2);
    int V_count = 200;
    for(int i = 0; i < V_count; i++) {
        V.emplace_back(i);
    }
    int pairs = 200;
    for(int i = 0; i < pairs; i++) {
        auto r1 = rand() % V_count;
        auto r2 = rand() % V_count;
        while(r1 == r2) {
            r2 = rand() % V_count;
        }
        auto p = std::make_pair(r1, r2);
        auto it = std::find(E.begin(), E.end(), p);
        if(it == E.end()) {
            E.emplace_back(p);
        }
    }

#endif

    std::cout << "generating tree..." << std::endl;
    mat_t left_of;


    std::cout << "generating directions..." << std::endl;
    expr_t direction;
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

    std::cout << "generating cross..." << E.size() << std::endl;
    expr_t cross;

    int count = 0;
    uint64_t delta_sum = 0;
    for(int i = 0 ; i < E.size(); ++i) {
        for(int j = i + 1; j < E.size(); ++j) {
            if (i == j) continue;
            uint64_t vi = E[i].first;
            uint64_t vj = E[i].second;
            uint64_t vk = E[j].first;
            uint64_t vl = E[j].second;

//            if(vi == vk || vi == vl || vj == vk  ||  vj == vl) continue;

            expr_t c = (
                    same_page_rule[i][j] >>
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
    solver.set_num_threads(12);


    std::vector<uint64_t > linear;
    std::cout << "linearizing..." << std::endl;
    linearize lz(linear);
    lz(target);
    std::cout << "tseitin..." << std::endl;
    tseitin tseitin(solver, linear);
    tseitin(target);

    std::cout << linear.size() << std::endl;
    solver.add_clause(std::vector<Lit>{Lit(tseitin.indexof(target), false)});

    solver.set_verbosity(0);
    lbool ret = solver.solve();
#if 0
    if(ret == l_True) {
        std::vector<lbool> model = solver.get_model();

        for(int i = 0 ; i < E.size(); ++i) {
            for (int j = i + 1; j < E.size(); ++j) {
                if (i == j) continue;
                uint64_t vi = E[i].first;
                uint64_t vj = E[i].second;
                uint64_t vk = E[j].first;
                uint64_t vl = E[j].second;

//                if (vi == vk || vi == vl || vj == vk || vj == vl) continue;

                uint64_t r = tseitin.indexof(same_page_rule[i][j]);
                std::cout << model[r] << " " << r << " " << i << " " << j << std::endl;
            }
        }

        uint64_t r = tseitin.indexof(all_pages);
        std::cout << model[r] << " " << r  << std::endl;

        for(int p = 0; p < num_pages; p++) {
            for(int e = 0; e < E.size(); e++) {
                uint64_t r = tseitin.indexof(page_assign[p][e]);
                if(model[r] == l_True) {
                    std::cout << model[r] << " " << p << " " << e << std::endl;
                }
            }
        }
    }
#endif
#endif
    return ret;
}
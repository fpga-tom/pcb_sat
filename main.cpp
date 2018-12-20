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
    for (uint64_t np = 5; np < 15; np++) {
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
    int V_count = 70;
    for(int i = 0; i < V_count; i++) {
        V.emplace_back(i);
    }
    int pairs = 70;
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

#if 0
    for(uint64_t j = 1; j < V.size(); ++j) {
        if(left_of.count(V[0]) == 0) {
            left_of.insert(std::pair<uint64_t , std::map<uint64_t ,expr_t>>(V[0], std::map<uint64_t, expr_t>()));
        }
        auto &var = left_of[V[0]];

        expr_t v{"dir_" + std::to_string(V[0]) + "_" + std::to_string(V[j])};
        vars.emplace_back(v);
        var.insert(std::pair<uint64_t , expr_t>(V[j], v));
        directions.emplace_back(v);
    }

    if(left_of.count(V[1]) == 0) {
        left_of.insert(std::pair<uint64_t , std::map<uint64_t ,expr_t>>(V[2], std::map<uint64_t, expr_t>()));
    }
    auto &var = left_of[V[1]];

    expr_t v{"dir_" + std::to_string(V[1]) + "_" + std::to_string(V[2])};
    vars.emplace_back(v);
    var.insert(std::pair<uint64_t , expr_t>(V[2], v));
    directions.emplace_back(v);

#else
    std::cout << "generating directions..." << std::endl;
    expr_t direction;
    for(uint64_t i = 0; i < V.size(); ++i) {
        for(uint64_t j = 0; j < V.size(); j++) {
            if (i == j) continue;
            if(left_of.count(V[i]) == 0) {
                left_of.insert(std::pair<uint64_t , std::map<uint64_t ,expr_t>>(V[i], std::map<uint64_t, expr_t>()));
            }
            auto &var = left_of[V[i]];

            if (V[i] > V[j]) {
                var.insert(std::pair<uint64_t , expr_t>(V[j], !left_of[V[j]][V[i]]) );
            } else {
                expr_t v(new expr_s("dir_" + std::to_string(V[i]) + "_" + std::to_string(V[j])));
                var.insert(std::pair<uint64_t , expr_t>(V[j], v));
                if(direction.get() == nullptr)
                    direction = v;
                else
                    direction = direction && v;
            }
        }
    }
#endif


    std::cout << "generating pages..." << std::endl;
    mat_t page_assign;

//    expr_t pages_formula;
    expr_t all_pages;
    for(int e = 0; e < E.size(); ++e) {
        expr_t page;
        for (int p = 0; p < num_pages; ++p) {
            if(page_assign.count(p) == 0) {
                page_assign.insert(std::pair<uint64_t, std::map<uint64_t, expr_t>>(p, std::map<uint64_t, expr_t>()));
            }
            auto &var = page_assign[p];

            expr_t v(new expr_s("page_" + std::to_string(p) + "_" + std::to_string(e)));
            var.insert(std::pair<uint64_t , expr_t>(e, v));
            if (page.get() == nullptr)
                page = v;
            else
                page = page || v;

        }
        if( e == 0)
            all_pages = page;
        else
            all_pages = all_pages && page;
    }


    std::cout << "generating same..." << std::endl;
    mat_t same_page_rule;

    for (int i = 0;i < E.size(); ++i) {
        for (int j = i + 1; j < E.size(); ++j) {
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
        for(int j = i+1; j < E.size(); ++j) {
            if (i == j) continue;
            uint64_t vi = E[i].first;
            uint64_t vj = E[i].second;
            uint64_t vk = E[j].first;
            uint64_t vl = E[j].second;

            if(vi == vk || vi == vl || vj == vk || vj == vl) continue;

            expr_t c = (
                    same_page_rule[i][j] >>
                    (
                            !(left_of[vi][vk] && left_of[vk][vj] && left_of[vj][vl]) &&
                            !(left_of[vi][vl] && left_of[vl][vj] && left_of[vj][vk]) &&
                            !(left_of[vj][vk] && left_of[vk][vi] && left_of[vi][vl]) &&
                            !(left_of[vj][vl] && left_of[vl][vi] && left_of[vi][vk]) &&
                            !(left_of[vk][vi] && left_of[vi][vl] && left_of[vl][vj]) &&
                            !(left_of[vk][vj] && left_of[vj][vl] && left_of[vl][vi]) &&
                            !(left_of[vl][vi] && left_of[vi][vk] && left_of[vk][vj]) &&
                            !(left_of[vl][vj] && left_of[vj][vk] && left_of[vk][vi])
                    ));

            if(cross.get() == nullptr)
                cross = c;
            else {
                cross = cross && c;
            }

        }
    }



    expr_t target = all_pages && direction && cross;

#if 0
    std::stringstream ss;

    ss << "import sympy as sp" << std::endl;
    ss << "from sympy.abc import *" << std::endl;

    for(const expr_t& vs : vars) {
        ss << vs << ", ";
    }

    ss << " = symbols('";

    for(const expr_t& vs : vars) {
        ss << vs << " ";
    }

    ss << "')" << std::endl;

    ss << "a = sp.to_cnf(";
    ss << target;
    ss << ", False)" << std::endl;
    ss << "with open('/tmp/cnf.txt', 'w') as f:" << std::endl;
    ss << "\tfor i in a.args:\n"
          " \t   print >> f, str(i).";
    for(int i = 0; i < vars.size(); i++) {
        ss << "replace('" << vars[i] << "','" << std::to_string(i) << "')";
        if (i < vars.size() - 1)
            ss << ".";
    }

    ss << ".replace(' ', '')" << std::endl;


    Py_SetProgramName(argv[0]);
    Py_Initialize();
    PyRun_SimpleString(ss.str().c_str());
    Py_Finalize();
#endif

    std::cout << "generating solver.." << std::endl;
    SATSolver solver;
    solver.set_num_threads(10);
#if 0
    solver.new_vars(vars.size());

    std::ifstream infile("/tmp/cnf.txt");
    std::string line;
    std::vector<Lit> clause;
    while(std::getline(infile, line)) {
        clause.clear();
        std::vector<std::string> l;
        boost::split(l, line, boost::is_any_of("|"));
        for(auto c : l) {
            if(c.front() == '~') {
                clause.emplace_back(Lit(std::stoi(c.substr(1)), true));
            } else {
                clause.emplace_back(Lit(std::stoi(c), false));
            }
        }
        solver.add_clause(clause);
    }
#endif


    std::vector<uint64_t > linear;
    std::cout << "linearizing..." << std::endl;
    linearize lz(linear);
    lz(target);
    std::cout << "tseitin..." << std::endl;
    tseitin tseitin(solver, linear);
    tseitin(target);

    std::cout << linear.size() << std::endl;
    solver.add_clause(std::vector<Lit>{Lit(0, false)});


    solver.set_verbosity(0);
    lbool ret = solver.solve();
//    solver.print_stats();
//    std::cout << ret << std::endl;

#endif
    return ret;
}
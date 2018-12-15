#include <iostream>
#include <cryptominisat5/cryptominisat.h>
#include <boost/algorithm/string.hpp>
#include <assert.h>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include "Expression.h"
#include <Python.h>

using std::vector;
using namespace CMSat;

typedef std::map<uint64_t , std::map<uint64_t , expr_t>> mat_t;

extern void demorgan();

int main(int argc, char** argv) {

    demorgan();
    uint64_t num_pages = 2;
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

    ExpressionFactory f;
    mat_t left_of;
    std::vector<expr_t> directions;
    std::vector<expr_t> vars;

    for(uint64_t i = 0; i < V.size(); ++i) {
        for(uint64_t j = 0; j < V.size(); j++) {
            if (i == j) continue;
            if(left_of.count(V[i]) == 0) {
                left_of.insert(std::pair<uint64_t , std::map<uint64_t ,expr_t>>(V[i], std::map<uint64_t, expr_t>()));
            }
            auto &var = left_of[V[i]];

            if (V[i] > V[j]) {
                expr_t v{"dir_" + std::to_string(V[i]) + "_" + std::to_string(V[j])};
                auto op = f._not(v);
                vars.emplace_back(v);
                var.insert(std::pair<uint64_t , expr_t>(V[j], op));
                directions.emplace_back(op);
            } else {
                expr_t v{"dir_" + std::to_string(V[i]) + "_" + std::to_string(V[j])};
                vars.emplace_back(v);
                var.insert(std::pair<uint64_t , expr_t>(V[j], v));
                directions.emplace_back(v);
            }
        }
    }

    expr_t direction = f._and(directions);

    mat_t page_assign;
    std::vector<expr_t> pages_formula;

    for(int e = 0; e < E.size(); ++e) {
        std::vector<expr_t> pages;
        for (int p = 0; p < num_pages; ++p) {
            if(page_assign.count(p) == 0) {
                page_assign.insert(std::pair<uint64_t, std::map<uint64_t, expr_t>>(p, std::map<uint64_t, expr_t>()));
            }
            auto &var = page_assign[p];

            std::string v("page_" + std::to_string(p) + "_" + std::to_string(e));
            vars.emplace_back(v);
            var.insert(std::pair<uint64_t , expr_t>(e, v));
            pages.emplace_back(v);

        }
        expr_t page = f._or(pages);
        pages_formula.emplace_back(page);
    }

    expr_t all_pages = f._and(pages_formula);

    mat_t same_page_rule;

    for (int i = 0;i < E.size(); ++i) {
        for (int j = i + 1; j < E.size(); ++j) {
            if (i == j) continue;
            std::vector<expr_t> rule;
            for(int p = 0; p < num_pages; ++p) {
                rule.emplace_back(f._and({page_assign[p][i], page_assign[p][j]}));
            }

            if(same_page_rule.count(i) == 0) {
                same_page_rule.insert(std::pair<uint64_t , std::map<uint64_t , expr_t>>(i, std::map<uint64_t, expr_t>()));
            }
            same_page_rule[i].insert(std::pair<uint64_t, expr_t>(j, f._or(rule)));
        }
    }

    std::vector<expr_t> cross;

    for(int i = 0 ; i < E.size(); ++i) {
        for(int j = i+1; j < E.size(); ++j) {
            if (i == j) continue;
            uint64_t vi = E[i].first;
            uint64_t vj = E[i].second;
            uint64_t vk = E[j].first;
            uint64_t vl = E[j].second;

            if(vi == vk || vi == vl || vj == vk || vj == vl) continue;


            cross.emplace_back(f._implication(
                    (same_page_rule[i][j]),
                    (f._and({
                                    f._not(f._and({
                                                          left_of[vi][vk],
                                                          left_of[vk][vj],
                                                          left_of[vj][vl]})),
                                    f._not(f._and({
                                                          left_of[vi][vl],
                                                          left_of[vl][vj],
                                                          left_of[vj][vk]})),
                                    f._not(f._and({
                                                          left_of[vj][vk],
                                                          left_of[vk][vi],
                                                          left_of[vi][vl]})),
                                    f._not(f._and({
                                                          left_of[vj][vl],
                                                          left_of[vl][vi],
                                                          left_of[vi][vk]})),
                                    f._not(f._and({
                                                          left_of[vk][vi],
                                                          left_of[vi][vl],
                                                          left_of[vl][vj]})),
                                    f._not(f._and({
                                                          left_of[vk][vj],
                                                          left_of[vj][vl],
                                                          left_of[vl][vi]})),
                                    f._not(f._and({
                                                          left_of[vl][vi],
                                                          left_of[vi][vk],
                                                          left_of[vk][vj]})),
                                    f._not(f._and({
                                                          left_of[vl][vj],
                                                          left_of[vj][vk],
                                                          left_of[vk][vi]}))
                            }))));
        }
    }


    expr_t all_cross = f._or(cross);
    expr_t target = f._and({all_pages, direction, all_cross});

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

    SATSolver solver;
    solver.set_num_threads(4);
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

    lbool ret = solver.solve();
    std::cout << ret << std::endl;

    return 0;
}
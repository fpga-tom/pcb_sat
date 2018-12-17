//
// Created by tomas on 12/15/18.
//

#include<vector>
#include <map>
#include"z3++.h"

using namespace z3;

typedef std::map<uint64_t , std::map<uint64_t , expr>> mat_t;

void demorgan() {
    std::cout << "de-Morgan example\n";

    context c;

    expr x = c.bool_const("x");
    expr y = c.bool_const("y");
    expr conjecture = (!(x && y)) == (!x || !y);

    solver s(c);
    // adding the negation of the conjecture as a constraint.
    s.add(!conjecture);
    std::cout << s << "\n";
    std::cout << s.to_smt2() << "\n";
    switch (s.check()) {
        case unsat:   std::cout << "de-Morgan is valid\n"; break;
        case sat:     std::cout << "de-Morgan is not valid\n"; break;
        case unknown: std::cout << "unknown\n"; break;
    }
}

class Z3ExpressionFactory {
public:
    expr _and_head(expr exp1, std::vector<expr> exp2) {
        if(exp2.empty()) {
            return exp1;
        }
        expr exp = _and_head(exp2.front(), std::vector<expr>(exp2.begin() + 1, exp2.end()));
        return (exp1 && exp);
    }
    expr _and(std::vector<expr> exp) { return _and_head(exp.front(), std::vector<expr>(exp.begin() + 1, exp.end())); }
    expr _and(std::initializer_list<expr> exp) { return _and(std::vector<expr>(exp)); }

    expr _or_head(expr exp1, std::vector<expr> exp2) {
        if(exp2.empty()) {
            return exp1;
        }
        expr exp = _or_head(exp2.front(), std::vector<expr>(exp2.begin() + 1, exp2.end()));
        return (exp1 || exp);
    }
    expr _or(std::vector<expr> exp) { return _or_head(exp.front(), std::vector<expr>(exp.begin() + 1, exp.end())); }
    expr _or(std::initializer_list<expr> exp) { return _or(std::vector<expr>(exp)); }

    expr _implication(expr premise, expr conclusion) { return _or({_not(premise), conclusion}); }
    expr _not(expr exp) { return (!exp); }
};

void book_embedding() {
    context c;

    uint64_t num_pages = 4;
#if 1
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
    std::vector<std::pair<uint64_t ,uint64_t >> E;

    srand(1);
    int V_count = 50;
    for(int i = 0; i < V_count; i++) {
        V.emplace_back(i);
    }
    int pairs = 100;
    for(int i = 0; i < pairs; i++) {
        auto r1 = rand() % V_count;
        auto r2 = rand() % V_count;
        while(r1 == r2) {
            r2 = rand() % V_count;
        }
        E.emplace_back(std::make_pair(r1, r2));
    }
#endif

    Z3ExpressionFactory f;
    mat_t left_of;
    std::vector<expr> directions;
    std::vector<expr> vars;

    for(uint64_t i = 0; i < V.size(); ++i) {
        for(uint64_t j = 0; j < V.size(); j++) {
            if (i == j) continue;
            if(left_of.count(V[i]) == 0) {
                left_of.insert(std::pair<uint64_t , std::map<uint64_t ,expr>>(V[i], std::map<uint64_t, expr>()));
            }
            auto &var = left_of[V[i]];

            if (V[i] > V[j]) {
                expr v = !c.bool_const(("dir_" + std::to_string(V[i]) + "_" + std::to_string(V[j])).c_str());
                vars.emplace_back(v);
                var.insert(std::pair<uint64_t , expr>(V[j], v));
                directions.emplace_back(v);
            } else {
                expr v = c.bool_const(("dir_" + std::to_string(V[i]) + "_" + std::to_string(V[j])).c_str());
                vars.emplace_back(v);
                var.insert(std::pair<uint64_t , expr>(V[j], v));
                directions.emplace_back(v);
            }
        }
    }

    expr direction = f._and(directions);

    mat_t page_assign;
    std::vector<expr> pages_formula;

    for(int e = 0; e < E.size(); ++e) {
        std::vector<expr> pages;
        for (int p = 0; p < num_pages; ++p) {
            if(page_assign.count(p) == 0) {
                page_assign.insert(std::pair<uint64_t, std::map<uint64_t, expr>>(p, std::map<uint64_t, expr>()));
            }
            auto &var = page_assign[p];

            expr v = c.bool_const(("page_" + std::to_string(p) + "_" + std::to_string(e)).c_str());
            vars.emplace_back(v);
            var.insert(std::pair<uint64_t , expr>(e, v));
            pages.emplace_back(v);

        }
        expr page = f._or(pages);
        pages_formula.emplace_back(page);
    }

    expr all_pages = f._and(pages_formula);

    mat_t same_page_rule;

    for (int i = 0;i < E.size(); ++i) {
        for (int j = i + 1; j < E.size(); ++j) {
            if (i == j) continue;
            std::vector<expr> rule;
            for(int p = 0; p < num_pages; ++p) {
                rule.emplace_back(page_assign.at(p).at(i) && page_assign.at(p).at(j));
            }

            if(same_page_rule.count(i) == 0) {
                same_page_rule.insert(std::pair<uint64_t , std::map<uint64_t , expr>>(i, std::map<uint64_t, expr>()));
            }
            same_page_rule[i].insert(std::pair<uint64_t, expr>(j, f._or(rule)));
        }
    }

    std::vector<expr> cross;

    for(int i = 0 ; i < E.size(); ++i) {
        for(int j = i+1; j < E.size(); ++j) {
            if (i == j) continue;
            uint64_t vi = E[i].first;
            uint64_t vj = E[i].second;
            uint64_t vk = E[j].first;
            uint64_t vl = E[j].second;

            if(vi == vk || vi == vl || vj == vk || vj == vl) continue;

            cross.emplace_back(f._implication(
                    (same_page_rule.at(i).at(j)),
                                    !(left_of.at(vi).at(vk) && left_of.at(vk).at(vj) && left_of.at(vj).at(vl)) &&
                                    !(left_of.at(vi).at(vl) && left_of.at(vl).at(vj) && left_of.at(vj).at(vk)) &&
                                    !(left_of.at(vj).at(vk) && left_of.at(vk).at(vi) && left_of.at(vi).at(vl)) &&
                                    !(left_of.at(vj).at(vl) && left_of.at(vl).at(vi) && left_of.at(vi).at(vk)) &&
                                    !(left_of.at(vk).at(vi) && left_of.at(vi).at(vl) && left_of.at(vl).at(vj)) &&
                                    !(left_of.at(vk).at(vj) && left_of.at(vj).at(vl) && left_of.at(vl).at(vi)) &&
                                    !(left_of.at(vl).at(vi) && left_of.at(vi).at(vk) && left_of.at(vk).at(vj)) &&
                                    !(left_of.at(vl).at(vj) && left_of.at(vj).at(vk) && left_of.at(vk).at(vi))
                            ));
        }
    }


    expr all_cross = f._and(cross);
    expr target = f._and({all_pages, direction, all_cross});
    std::cout << target.to_string() << std::endl;

    solver s(c);
    s.add(target);
    switch (s.check()) {
        case unsat:   std::cout << "false\n"; break;
        case sat:     std::cout << "true\n"; break;
        case unknown: std::cout << "unknown\n"; break;
    }
}
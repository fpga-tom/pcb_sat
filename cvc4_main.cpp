/*********************                                                        */
/*! \file simple_vc_cxx.cpp
 ** \verbatim
 ** Top contributors (to current version):
 **   Morgan Deters, Dejan Jovanovic
 ** This file is part of the CVC4 project.
 ** Copyright (c) 2009-2018 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file COPYING in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief A simple demonstration of the C++ interface
 **
 ** A simple demonstration of the C++ interface.  Compare to the Java
 ** interface in SimpleVC.java; they are virtually line-by-line
 ** identical.
 **/

#include <iostream>

#include <cvc4/cvc4.h> // use this after CVC4 is properly installed
//#include "smt/smt_engine.h"

using namespace std;
using namespace CVC4;

int cvc4_main() {
    ExprManager em;
    SmtEngine smt(&em);

    // Prove that for integers x and y:
    //   x > 0 AND y > 0  =>  2x + y >= 3

    Type integer = em.integerType();

    Expr x = em.mkVar("x", integer);
    Expr y = em.mkVar("y", integer);
    Expr zero = em.mkConst(Rational(0));

    Expr x_positive = em.mkExpr(kind::GT, x, zero);
    Expr y_positive = em.mkExpr(kind::GT, y, zero);

    Expr two = em.mkConst(Rational(2));
    Expr twox = em.mkExpr(kind::MULT, two, x);
    Expr twox_plus_y = em.mkExpr(kind::PLUS, twox, y);

    Expr three = em.mkConst(Rational(3));
    Expr twox_plus_y_geq_3 = em.mkExpr(kind::GEQ, twox_plus_y, three);

    Expr formula =
            em.mkExpr(kind::AND, x_positive, y_positive).
                    impExpr(twox_plus_y_geq_3);

    cout << "Checking validity of formula " << formula << " with CVC4." << endl;
    cout << "CVC4 should report VALID." << endl;
    cout << "Result from CVC4 is: " << smt.query(formula) << endl;

    return 0;
}

typedef std::map<uint64_t , std::map<uint64_t , Expr>> mat_t;
static ExprManager em;

class CVC4ExpressionFactory {
public:
    Expr _and_head(Expr exp1, std::vector<Expr> exp2) {
        if(exp2.empty()) {
            return exp1;
        }
        Expr exp = _and_head(exp2.front(), std::vector<Expr>(exp2.begin() + 1, exp2.end()));
        return em.mkExpr(kind::AND, exp1, exp);
    }
    Expr _and(std::vector<Expr> exp) { return _and_head(exp.front(), std::vector<Expr>(exp.begin() + 1, exp.end())); }
    Expr _and(std::initializer_list<Expr> exp) { return _and(std::vector<Expr>(exp)); }

    Expr _or_head(Expr exp1, std::vector<Expr> exp2) {
        if(exp2.empty()) {
            return exp1;
        }
        Expr exp = _or_head(exp2.front(), std::vector<Expr>(exp2.begin() + 1, exp2.end()));
        return em.mkExpr(kind::OR, exp1, exp);
    }
    Expr _or(std::vector<Expr> exp) { return _or_head(exp.front(), std::vector<Expr>(exp.begin() + 1, exp.end())); }
    Expr _or(std::initializer_list<Expr> exp) { return _or(std::vector<Expr>(exp)); }

    Expr _implication(Expr premise, Expr conclusion) { return _or({_not(premise), conclusion}); }
    Expr _not(Expr exp) { return em.mkExpr(kind::NOT, exp); }
};

void cvc4_book_embedding() {

    uint64_t num_pages = 30;
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
    std::vector<std::pair<uint64_t ,uint64_t >> E;

    srand(2);
    int V_count = 100;
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

    CVC4ExpressionFactory f;
    mat_t left_of;
    Type _bool = em.booleanType();
    std::vector<Expr> directions;
    std::vector<Expr> vars;

    for(uint64_t i = 0; i < V.size(); ++i) {
        for(uint64_t j = 0; j < V.size(); j++) {
            if (i == j) continue;
            if(left_of.count(V[i]) == 0) {
                left_of.insert(std::pair<uint64_t , std::map<uint64_t ,Expr>>(V[i], std::map<uint64_t, Expr>()));
            }
            auto &var = left_of[V[i]];

            if (V[i] > V[j]) {
                Expr v = f._not(em.mkVar("dir_" + std::to_string(V[i]) + "_" + std::to_string(V[j]), _bool));
                vars.emplace_back(v);
                var.insert(std::pair<uint64_t , Expr>(V[j], v));
                directions.emplace_back(v);
            } else {
                Expr v = em.mkVar("dir_" + std::to_string(V[i]) + "_" + std::to_string(V[j]), _bool);
                vars.emplace_back(v);
                var.insert(std::pair<uint64_t , Expr>(V[j], v));
                directions.emplace_back(v);
            }
        }
    }

    Expr direction = f._and(directions);

    mat_t page_assign;
    std::vector<Expr> pages_formula;

    for(int e = 0; e < E.size(); ++e) {
        std::vector<Expr> pages;
        for (int p = 0; p < num_pages; ++p) {
            if(page_assign.count(p) == 0) {
                page_assign.insert(std::pair<uint64_t, std::map<uint64_t, Expr>>(p, std::map<uint64_t, Expr>()));
            }
            auto &var = page_assign[p];

            Expr v = em.mkVar("page_" + std::to_string(p) + "_" + std::to_string(e), _bool);
            vars.emplace_back(v);
            var.insert(std::pair<uint64_t , Expr>(e, v));
            pages.emplace_back(v);

        }
        Expr page = f._or(pages);
        pages_formula.emplace_back(page);
    }

    Expr all_pages = f._and(pages_formula);

    mat_t same_page_rule;

    for (int i = 0;i < E.size(); ++i) {
        for (int j = i + 1; j < E.size(); ++j) {
            if (i == j) continue;
            std::vector<Expr> rule;
            for(int p = 0; p < num_pages; ++p) {
                rule.emplace_back(f._and({page_assign[p][i], page_assign[p][j]}));
            }

            if(same_page_rule.count(i) == 0) {
                same_page_rule.insert(std::pair<uint64_t , std::map<uint64_t , Expr>>(i, std::map<uint64_t, Expr>()));
            }
            same_page_rule[i].insert(std::pair<uint64_t, Expr>(j, f._or(rule)));
        }
    }

    std::vector<Expr> cross;

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


    Expr all_cross = f._and(cross);
    Expr target = f._and({all_pages, direction, all_cross});

    SmtEngine smt(&em);

    std::cout << Configuration::isBuiltWithCryptominisat() << std::endl;
    std::cout << smt.checkSat(target) << std::endl;
}

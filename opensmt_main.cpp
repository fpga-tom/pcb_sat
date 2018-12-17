//
// Created by tomas on 12/16/18.
//

#include <opensmt/opensmt2.h>


typedef std::map<uint64_t , std::map<uint64_t , PTRef>> mat_t;

class SMTPTRefessionFactory {
    Logic& logic;
public:
    SMTPTRefessionFactory(Logic& logic) : logic(logic) {}
    PTRef _and_head(PTRef exp1, std::vector<PTRef> exp2) {
        if(exp2.empty()) {
            return exp1;
        }
        PTRef exp = _and_head(exp2.front(), std::vector<PTRef>(exp2.begin() + 1, exp2.end()));
        return logic.mkAnd(exp1, exp);
    }
    PTRef _and(std::vector<PTRef> exp) { return _and_head(exp.front(), std::vector<PTRef>(exp.begin() + 1, exp.end())); }
    PTRef _and(std::initializer_list<PTRef> exp) { return _and(std::vector<PTRef>(exp)); }

    PTRef _or_head(PTRef exp1, std::vector<PTRef> exp2) {
        if(exp2.empty()) {
            return exp1;
        }
        PTRef exp = _or_head(exp2.front(), std::vector<PTRef>(exp2.begin() + 1, exp2.end()));
        return logic.mkOr(exp1, exp);
    }
    PTRef _or(std::vector<PTRef> exp) { return _or_head(exp.front(), std::vector<PTRef>(exp.begin() + 1, exp.end())); }
    PTRef _or(std::initializer_list<PTRef> exp) { return _or(std::vector<PTRef>(exp)); }

    PTRef _implication(PTRef premise, PTRef conclusion) { return _or({_not(premise), conclusion}); }
    PTRef _not(PTRef exp) { return logic.mkNot(exp); }
};

void opensmt2_book_embedding() {
    const char *msg;
    SMTConfig c;
//    c.sat_split_threads(6);
//    c.setOption(SMTConfig::o_dump_query, SMTOption(1), msg);
//    c.setOption(SMTConfig::o_sat_dump_learnts, SMTOption(1), msg);
    UFTheory uftheory(c);
    THandler thandler(c, uftheory);
    SimpSMTSolver solver(c, thandler);
    MainSolver mainSolver(thandler, c, &solver, "main_solver");

    Logic& logic = thandler.getLogic();
    SMTPTRefessionFactory f(logic);


    uint64_t num_pages = 13;
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


    mat_t left_of;
    std::vector<PTRef> directions;
    std::vector<PTRef> vars;

    for(uint64_t i = 0; i < V.size(); ++i) {
        for(uint64_t j = 0; j < V.size(); j++) {
            if (i == j) continue;
            if(left_of.count(V[i]) == 0) {
                left_of.insert(std::pair<uint64_t , std::map<uint64_t ,PTRef>>(V[i], std::map<uint64_t, PTRef>()));
            }
            auto &var = left_of[V[i]];

            if (V[i] > V[j]) {
                PTRef v = f._not(logic.mkBoolVar(("dir_" + std::to_string(V[i]) + "_" + std::to_string(V[j])).c_str()));
                vars.emplace_back(v);
                var.insert(std::pair<uint64_t , PTRef>(V[j], v));
                directions.emplace_back(v);
            } else {
                PTRef v = logic.mkBoolVar(("dir_" + std::to_string(V[i]) + "_" + std::to_string(V[j])).c_str());
                vars.emplace_back(v);
                var.insert(std::pair<uint64_t , PTRef>(V[j], v));
                directions.emplace_back(v);
            }
        }
    }

    PTRef direction = f._and(directions);

    mat_t page_assign;
    std::vector<PTRef> pages_formula;

    for(int e = 0; e < E.size(); ++e) {
        std::vector<PTRef> pages;
        for (int p = 0; p < num_pages; ++p) {
            if(page_assign.count(p) == 0) {
                page_assign.insert(std::pair<uint64_t, std::map<uint64_t, PTRef>>(p, std::map<uint64_t, PTRef>()));
            }
            auto &var = page_assign[p];

            PTRef v = logic.mkBoolVar(("page_" + std::to_string(p) + "_" + std::to_string(e)).c_str());
            vars.emplace_back(v);
            var.insert(std::pair<uint64_t , PTRef>(e, v));
            pages.emplace_back(v);

        }
        PTRef page = f._or(pages);
        pages_formula.emplace_back(page);
    }

    PTRef all_pages = f._and(pages_formula);

    mat_t same_page_rule;

    for (int i = 0;i < E.size(); ++i) {
        for (int j = i + 1; j < E.size(); ++j) {
            if (i == j) continue;
            std::vector<PTRef> rule;
            for(int p = 0; p < num_pages; ++p) {
                rule.emplace_back(f._and({page_assign[p][i], page_assign[p][j]}));
            }

            if(same_page_rule.count(i) == 0) {
                same_page_rule.insert(std::pair<uint64_t , std::map<uint64_t , PTRef>>(i, std::map<uint64_t, PTRef>()));
            }
            same_page_rule[i].insert(std::pair<uint64_t, PTRef>(j, f._or(rule)));
        }
    }

    std::vector<PTRef> cross;

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


    PTRef all_cross = f._and(cross);
    PTRef target = f._and({all_pages, direction, all_cross});

    sstat in = mainSolver.push(target);
    if (in == l_False)    std::cout << "false\n";
    if (in == l_True)     std::cout << "true\n";
    if (in == l_Undef) std::cout << "unknown\n";

    char **msg1;
//    mainSolver.simplifyFormulas();
    sstat r = mainSolver.check();
    mainSolver.writeSolverState("/tmp/solver1.cnf", msg1);
    mainSolver.writeCnfState("/tmp/solver.cnf", msg1);



    if (r == l_False)    std::cout << "false\n";
    if (r == l_True)     std::cout << "true\n";
    if (r == l_Undef) std::cout << "unknown\n";

}

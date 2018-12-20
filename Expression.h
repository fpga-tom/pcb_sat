//
// Created by tomas on 12/15/18.
//

#ifndef PCB_SAT_expr_tESSION_H
#define PCB_SAT_expr_tESSION_H

#include <string>
#include <vector>
#include <boost/variant.hpp>
#include <boost/spirit/home/x3/support/ast/variant.hpp>
#include <boost/fusion/adapted/struct/adapt_struct.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <cryptominisat5/cryptominisat.h>
#include <map>
#include <numeric>


struct expr_s;
typedef std::shared_ptr<expr_s> expr_t;
enum class expr_type_e {_var, _op_not, _op_and, _op_or};
struct expr_s {

    expr_s(const expr_t oper1, expr_type_e type) : oper1(oper1), type(type) {}
    explicit expr_s(const std::string& var) : var(var), type(expr_type_e::_var) {}
    expr_s(const expr_t oper1, const expr_t oper2, expr_type_e type) : oper1(oper1), oper2(oper2), type(type) {}

    expr_t oper1;
    expr_t oper2;
    const std::string var;
    expr_type_e type;

};

struct linearize
{
    linearize(std::vector<uint64_t >& os) : _os(os) {}
    std::vector<uint64_t >& _os;

    void operator()(expr_t& ex) const {
        if(ex->type == expr_type_e::_var) {
            operator_var(ex);
        } else if(ex->type == expr_type_e::_op_and) {
            operator_and(ex);
        } else if(ex->type == expr_type_e::_op_or) {
            operator_or(ex);
        } else if(ex->type == expr_type_e::_op_not) {
            operator_not(ex);
        }
    }

    void operator_var(expr_t& v) const { _os.emplace_back((uint64_t )v.get()); }

    void operator_and(expr_t& b) const {
        _os.emplace_back((uint64_t )b.get());
        this->operator()(b->oper1);
        this->operator()(b->oper2);
    }
    void operator_or(expr_t& b) const {
        _os.emplace_back((uint64_t )b.get());
        this->operator()(b->oper1);
        this->operator()(b->oper2);
    }

    void operator_not(expr_t& u) const
    {
        _os.emplace_back((uint64_t )u.get());
        this->operator()(u->oper1);
    }
};
struct tseitin
{
    std::vector<uint64_t >& _lz;
    std::map<uint64_t , uint64_t > _index;
    CMSat::SATSolver& _os;

    tseitin(CMSat::SATSolver& os, std::vector<uint64_t >& lz) : _os(os), _lz(lz) {
        _os.new_vars(lz.size());
        for(int i = 0; i < lz.size(); ++i) {
            _index.insert(std::make_pair(lz[i], i));
        }

    }


    uint64_t indexof(const expr_t&  ex) const {
        return _index.at((uint64_t )ex.get());
    }

    void operator_var(expr_t& v) const {  }

    void operator()(expr_t& ex) const {
        if(ex->type == expr_type_e::_var) {
            operator_var(ex);
        } else if(ex->type == expr_type_e::_op_and) {
            operator_and(ex);
        } else if(ex->type == expr_type_e::_op_or) {
            operator_or(ex);
        } else if(ex->type == expr_type_e::_op_not) {
            operator_not(ex);
        }
    }

    void operator_and(expr_t& b) const {
        std::vector<CMSat::Lit> clause;

        uint64_t _a = indexof(b);
        uint64_t _b = indexof(b->oper1);
        uint64_t _c = indexof(b->oper2);

        this->operator()(b->oper1);
        this->operator()(b->oper2);


        clause.emplace_back(CMSat::Lit(_a, true));
        clause.emplace_back(CMSat::Lit(_b, false));
        _os.add_clause(clause);
        clause.clear();

        clause.emplace_back(CMSat::Lit(_a, true));
        clause.emplace_back(CMSat::Lit(_c, false));
        _os.add_clause(clause);
        clause.clear();

        clause.emplace_back(CMSat::Lit(_a, false));
        clause.emplace_back(CMSat::Lit(_b, true));
        clause.emplace_back(CMSat::Lit(_c, true));
        _os.add_clause(clause);
        clause.clear();
    }
    void operator_or(expr_t& b) const {
        std::vector<CMSat::Lit> clause;
        uint64_t _a = indexof(b);
        uint64_t _b = indexof(b->oper1);
        uint64_t _c = indexof(b->oper2);

        this->operator()(b->oper1);
        this->operator()(b->oper2);


        clause.emplace_back(CMSat::Lit(_a, false));
        clause.emplace_back(CMSat::Lit(_b, true));
        _os.add_clause(clause);
        clause.clear();

        clause.emplace_back(CMSat::Lit(_a, false));
        clause.emplace_back(CMSat::Lit(_c, true));
        _os.add_clause(clause);
        clause.clear();

        clause.emplace_back(CMSat::Lit(_a, true));
        clause.emplace_back(CMSat::Lit(_b, false));
        clause.emplace_back(CMSat::Lit(_c, false));
        _os.add_clause(clause);
        clause.clear();
    }

    void operator_not(expr_t& u) const
    {
        std::vector<CMSat::Lit> clause;

        uint64_t _a = indexof(u);
        uint64_t _b = indexof(u->oper1);

        this->operator()(u->oper1);

        clause.emplace_back(CMSat::Lit(_a, true));
        clause.emplace_back(CMSat::Lit(_b, true));
        _os.add_clause(clause);
        clause.clear();

        clause.emplace_back(CMSat::Lit(_a, false));
        clause.emplace_back(CMSat::Lit(_b, false));
        _os.add_clause(clause);
        clause.clear();
    }
};



expr_t operator&&(const expr_t& a, const expr_t& b) {
    assert(a.get() != nullptr);
    assert(b.get() != nullptr);
    return expr_t(new expr_s(a, b, expr_type_e::_op_and));
}

expr_t operator||(const expr_t& a, const expr_t& b) {
    assert(a.get() != nullptr);
    assert(b.get() != nullptr);
    return expr_t(new expr_s(a, b, expr_type_e::_op_or));
}

expr_t operator!(const expr_t& a) {
    assert(a.get() != nullptr);
    return expr_t(new expr_s(a, expr_type_e::_op_not));
}

expr_t operator>>(const expr_t& a,const expr_t& b) {
    assert(a.get() != nullptr);
    assert(b.get() != nullptr);
    return !a || b;
}






#endif //PCB_SAT_expr_tESSION_H

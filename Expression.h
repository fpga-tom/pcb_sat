//
// Created by tomas on 12/15/18.
//

#ifndef PCB_SAT_expr_tESSION_H
#define PCB_SAT_expr_tESSION_H

#include <string>
#include <vector>
#include <boost/variant.hpp>
#include <cryptominisat5/cryptominisat.h>
#include <map>


struct op_or  { }; // tag
struct op_and { }; // tag
struct op_not { }; // tag

typedef std::string var;
template <typename tag> struct binop;
template <typename tag> struct unop;

typedef boost::variant<var,
boost::recursive_wrapper<unop <op_not> >,
boost::recursive_wrapper<binop<op_and> >,
boost::recursive_wrapper<binop<op_or> >
> expr_t;



template <typename tag> struct binop
{
    explicit binop(const expr_t& l, const expr_t& r, uint32_t id) : oper1(l), oper2(r), id(id) { }
    expr_t oper1, oper2;
    uint32_t id;

    bool operator==(const binop& other) const {
        return this->id == other.id;
    }

    bool operator<(const binop& other) const {
        return this->id < other.id;
    }
};

template <typename tag> struct unop
{
    explicit unop(const expr_t& o, uint32_t id) : oper1(o), id(id) { }
    expr_t oper1;
    uint32_t id;

    bool operator==(const unop& other) const {
        return this->id == other.id;
    }

    bool operator<(const unop& other) const {
        return this->id < other.id;
    }
};

struct printer : boost::static_visitor<void>
{
    printer(std::ostream& os) : _os(os) {}
    std::ostream& _os;

    //
    void operator()(const var& v) const { _os << v; }

    void operator()(const binop<op_and>& b) const { print(" & ", b.oper1, b.oper2); }
    void operator()(const binop<op_or >& b) const { print(" | ", b.oper1, b.oper2); }

    void print(const std::string& op, const expr_t& l, const expr_t& r) const
    {
        _os << "(";
        boost::apply_visitor(*this, l);
        _os << op;
        boost::apply_visitor(*this, r);
        _os << ")";
    }

    void operator()(const unop<op_not>& u) const
    {
        _os << "(";
        _os << "~";
        boost::apply_visitor(*this, u.oper1);
        _os << ")";
    }
};

std::ostream& operator<<(std::ostream& os, const expr_t& e);

struct linearize : boost::static_visitor<void >
{
    linearize(std::vector<expr_t>& os) : _os(os) {}
    std::vector<expr_t>& _os;

    //
    void operator()(const var& v) const { _os.emplace_back(v); }

    void operator()(const binop<op_and>& b) const {
        _os.emplace_back(b);
        boost::apply_visitor(*this, b.oper1);
        boost::apply_visitor(*this, b.oper2);
    }
    void operator()(const binop<op_or >& b) const {
        _os.emplace_back(b);
        boost::apply_visitor(*this, b.oper1);
        boost::apply_visitor(*this, b.oper2);
    }

    void operator()(const unop<op_not>& u) const
    {
        _os.emplace_back(u);
        boost::apply_visitor(*this, u.oper1);
    }
};
struct tseitin : boost::static_visitor<void >
{
    std::vector<expr_t>& _lz;
    std::map<expr_t, uint32_t > _index;
    CMSat::SATSolver& _os;

    tseitin(CMSat::SATSolver& os, std::vector<expr_t>& lz) : _os(os), _lz(lz) {
        _os.new_vars(lz.size());
        for(int i = 0; i < lz.size(); ++i) {
            _index.insert(std::make_pair(lz[i], i));
        }

    }


    uint32_t indexof(const expr_t& ex) const {
        return _index.at(ex);
//        auto it = std::find(_lz.begin(), _lz.end(), ex);
//        if(it == _lz.end()) {
//            std::cerr << "not found" << std::endl;
//        }
//        return std::distance(_lz.begin(), it);
    }
    //
    void operator()(const var& v) const {  }

    void operator()(const binop<op_and>& b) const {
        std::vector<CMSat::Lit> clause;

        uint32_t _a = indexof(b);
        uint32_t _b = indexof(b.oper1);
        uint32_t _c = indexof(b.oper2);

        boost::apply_visitor(*this, b.oper1);
        boost::apply_visitor(*this, b.oper2);


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
    void operator()(const binop<op_or >& b) const {
        std::vector<CMSat::Lit> clause;
        uint32_t _a = indexof(b);
        uint32_t _b = indexof(b.oper1);
        uint32_t _c = indexof(b.oper2);

        boost::apply_visitor(*this, b.oper1);
        boost::apply_visitor(*this, b.oper2);


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

    void operator()(const unop<op_not>& u) const
    {
        std::vector<CMSat::Lit> clause;

        uint32_t _a = indexof(u);
        uint32_t _b = indexof(u.oper1);

        boost::apply_visitor(*this, u.oper1);

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


class ExpressionFactory {
    int count;
public:
    ExpressionFactory() : count(0) {}
    expr_t _and_head(expr_t exp1, std::vector<expr_t> exp2) {
        if(exp2.empty()) {
            return exp1;
        }
        expr_t exp = _and_head(exp2.front(), std::vector<expr_t>(exp2.begin() + 1, exp2.end()));
        return binop<op_and>(exp1, exp, count++);
    }
    expr_t _and(std::vector<expr_t> exp) { return _and_head(exp.front(), std::vector<expr_t>(exp.begin() + 1, exp.end())); }
    expr_t _and(std::initializer_list<expr_t> exp) { return _and(std::vector<expr_t>(exp)); }

    expr_t _or_head(expr_t exp1, std::vector<expr_t> exp2) {
        if(exp2.empty()) {
            return exp1;
        }
        expr_t exp = _or_head(exp2.front(), std::vector<expr_t>(exp2.begin() + 1, exp2.end()));
        return binop<op_or>(exp1, exp, count++);
    }
    expr_t _or(std::vector<expr_t> exp) { return _or_head(exp.front(), std::vector<expr_t>(exp.begin() + 1, exp.end())); }
    expr_t _or(std::initializer_list<expr_t> exp) { return _or(std::vector<expr_t>(exp)); }

    expr_t _implication(expr_t premise, expr_t conclusion) { return _or({_not(premise), conclusion}); }
    expr_t _not(expr_t exp) { return unop<op_not>(exp, count++); }
};


#endif //PCB_SAT_expr_tESSION_H

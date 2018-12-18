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
#include <numeric>


struct op_or  { }; // tag
struct op_and { }; // tag
struct op_not { }; // tag

typedef struct _var {
    _var() {}
    _var(const std::string& str, uint32_t id) : str(str), id(id) {}
    std::string str;
    uint32_t id;

    bool operator==(const _var& other) const {
        return this->id == other.id;
    }

    bool operator<(const _var& other) const {
        return this->id < other.id;
    }

    operator int32_t() {
        return id;
    }
} var;
template <typename tag> struct binop;
template <typename tag> struct unop;

typedef boost::variant<var,
boost::recursive_wrapper<unop <op_not> >,
boost::recursive_wrapper<binop<op_and> >,
boost::recursive_wrapper<binop<op_or> >
> expr_t;



template <typename tag> struct binop
{
    binop(const expr_t&& l, const expr_t&& r, uint32_t id) : oper1(l), oper2(r), id(id) { }

//    explicit binop(expr_t&& l, expr_t&& r, uint32_t id) : oper1(l), oper2(r), id(id) { }
    /*
    binop(const binop&& other) : oper1(std::move(other.oper1)), oper2(std::move(other.oper2)), id(std::move(other.id)) {}
    binop(const binop& other) : oper1(other.oper1), oper2(other.oper2), id(other.id) {}
    binop<tag>& operator=(binop<tag>&& other)  {
        oper1 = std::move(other.oper1);
        oper2 = std::move(other.oper2);
        id = std::move(other.id);
        return *this;
    }
     */
    expr_t oper1, oper2;
    uint32_t id;

    bool operator==(const binop& other) const {
        return this->id == other.id;
    }

    bool operator<(const binop& other) const {
        return this->id < other.id;
    }

    operator int32_t() {
        return id;
    }
};

template <typename tag> struct unop
{
    unop(const expr_t&& o, uint32_t id) : oper1(o), id(id) { }
//    explicit unop(expr_t&& o, uint32_t id) : oper1(o), id(id) { }
    /*
    unop(const unop&& other) : oper1(std::move(other.oper1)), id(std::move(other.id)) {}
    unop(const unop& other) : oper1(std::move(other.oper1)), id(std::move(other.id)) {}
     unop<tag>& operator=(unop<tag>&& other)  {
        oper1 = std::move(other.oper1);
        id = std::move(other.id);
        return *this;
    }
     */
    expr_t oper1;
    uint32_t id;



    bool operator==(const unop& other) const {
        return this->id == other.id;
    }

    bool operator<(const unop& other) const {
        return this->id < other.id;
    }

    operator int32_t() {
        return id;
    }
};

/*
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
 */

std::ostream& operator<<(std::ostream& os, const expr_t& e);

struct linearize : boost::static_visitor<void >
{
    linearize(std::vector<uint32_t >& os) : _os(os) {}
    std::vector<uint32_t >& _os;

    //
    void operator()(const var& v) const { _os.emplace_back(v.id); }

    void operator()(const binop<op_and>& b) const {
        _os.emplace_back(b.id);
        boost::apply_visitor(*this, b.oper1);
        boost::apply_visitor(*this, b.oper2);
    }
    void operator()(const binop<op_or >& b) const {
        _os.emplace_back(b.id);
        boost::apply_visitor(*this, b.oper1);
        boost::apply_visitor(*this, b.oper2);
    }

    void operator()(const unop<op_not>& u) const
    {
        _os.emplace_back(u.id);
        boost::apply_visitor(*this, u.oper1);
    }
};
struct tseitin : boost::static_visitor<void >
{
    std::vector<uint32_t >& _lz;
    std::map<uint32_t , uint32_t > _index;
    CMSat::SATSolver& _os;

    tseitin(CMSat::SATSolver& os, std::vector<uint32_t >& lz) : _os(os), _lz(lz) {
        _os.new_vars(lz.size());
        for(int i = 0; i < lz.size(); ++i) {
            _index.insert(std::make_pair(lz[i], i));
        }

    }


    uint32_t indexof(const expr_t&  ex) const {
        return _index.at(boost::apply_visitor(find_id(), ex));
//        auto it = std::find(_lz.begin(), _lz.end(), ex);
//        if(it == _lz.end()) {
//            std::cerr << "not found" << std::endl;
//        }
//        return std::distance(_lz.begin(), it);
    }
    //
    void operator()(const var& v) const {  }

    struct find_id : public boost::static_visitor<uint32_t > {

        uint32_t operator()(const var& v) const { return v.id; }
        uint32_t operator()(const binop<op_and>& b) const { return b.id; }
        uint32_t operator()(const binop<op_or >& b) const { return b.id; }
        uint32_t operator()(const unop<op_not>& u) const { return u.id; }
    };

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

public:
    int count;
    ExpressionFactory() : count(0) {}
    expr_t _and_head(const expr_t& exp1, const std::vector<expr_t>::const_iterator exp2, const std::vector<expr_t>::const_iterator exp2_end) {
        if(exp2 == exp2_end) {
            return exp1;
        }
        const expr_t exp = _and_head(*exp2, std::next(exp2), exp2_end);
        return binop<op_and>(std::move(exp1), std::move(exp), count++);
    }
    expr_t _and(const std::vector<expr_t>& exp) { return _and_head(exp.front(), std::next(exp.begin()), exp.end()); }
    expr_t _and(const std::initializer_list<expr_t>& exp) { return _and(std::vector<expr_t>(exp)); }

    expr_t _or_head(const expr_t& exp1, const std::vector<expr_t>::const_iterator exp2, const std::vector<expr_t>::const_iterator exp2_end) {
        if(exp2 == exp2_end) {
            return exp1;
        }
        const expr_t exp = _or_head(*exp2, std::next(exp2), exp2_end);
        return binop<op_or>(std::move(exp1), std::move(exp), count++);
    }
    expr_t _or(const std::vector<expr_t>& exp) { return _or_head(exp.front(), std::next(exp.begin()), exp.end()); }
    expr_t _or(const std::initializer_list<expr_t>& exp) { return _or(std::vector<expr_t>(exp)); }
//    expr_t _and(std::vector<expr_t>&& exp) { return std::accumulate(std::next(exp.begin()), exp.end(), exp[0], [this](expr_t& a, expr_t& b) {
//        return binop<op_and>(std::move(a), std::move(b), count++);
//    }); }
//    expr_t _and(const std::vector<expr_t>& exp) { return std::accumulate(std::next(exp.begin()), exp.end(), exp[0], [this](const expr_t& a, const expr_t& b) {
//            return binop<op_and>(std::move(a),std::move(b), count++);
//        }); }
//    expr_t _and(const std::initializer_list<expr_t>& exp) { return std::accumulate(std::next(exp.begin()), exp.end(), *exp.begin(), [this](const expr_t& a, const expr_t& b) {
//            return binop<op_and>(std::move(a), std::move(b), count++);
//        }); }

//    expr_t _or_head(const expr_t& exp1, const std::vector<expr_t>::iterator& exp2, const std::vector<expr_t>::iterator& exp2_end) {
//        if(exp2 == exp2_end) {
//            return exp1;
//        }
//        expr_t exp = _or_head(*exp2, std::next(exp2), exp2_end);
//        return binop<op_or>(exp1, exp, count++);
//    }
//    expr_t _or(std::vector<expr_t> exp) { return _or_head(exp.front(), std::next(exp.begin()), exp.end()); }
//    expr_t _or(const std::vector<expr_t>& exp) { return std::accumulate(std::next(exp.begin()), exp.end(), exp[0], [this](const expr_t& a, const expr_t& b) {
//            return binop<op_or>(std::move(a), std::move(b), count++);
//        }); }
//    expr_t _or(const std::initializer_list<expr_t>& exp) { return std::accumulate(std::next(exp.begin()), exp.end(), *exp.begin(), [this](const expr_t& a, const expr_t& b) {
//            return binop<op_or>(std::move(a), std::move(b), count++);
//        }); }

    expr_t _implication(const expr_t& premise, const expr_t& conclusion) { return _or({_not(premise), conclusion}); }
    expr_t _not(const expr_t& exp) { return unop<op_not>(std::move(exp), count++); }
    expr_t _var(const std::string& v) { return var(v, count++); }
};


#endif //PCB_SAT_expr_tESSION_H

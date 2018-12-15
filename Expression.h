//
// Created by tomas on 12/15/18.
//

#ifndef PCB_SAT_expr_tESSION_H
#define PCB_SAT_expr_tESSION_H

#include <string>
#include <vector>
#include <boost/variant.hpp>

struct op_or  {}; // tag
struct op_and {}; // tag
struct op_xor {}; // tag
struct op_not {}; // tag

typedef std::string var;
template <typename tag> struct binop;
template <typename tag> struct unop;

typedef boost::variant<var,
boost::recursive_wrapper<unop <op_not> >,
boost::recursive_wrapper<binop<op_and> >,
boost::recursive_wrapper<binop<op_xor> >,
boost::recursive_wrapper<binop<op_or> >
> expr_t;


template <typename tag> struct binop
{
    explicit binop(const expr_t& l, const expr_t& r) : oper1(l), oper2(r) { }
    expr_t oper1, oper2;
};

template <typename tag> struct unop
{
    explicit unop(const expr_t& o) : oper1(o) { }
    expr_t oper1;
};

struct printer : boost::static_visitor<void>
{
    printer(std::ostream& os) : _os(os) {}
    std::ostream& _os;

    //
    void operator()(const var& v) const { _os << v; }

    void operator()(const binop<op_and>& b) const { print(" & ", b.oper1, b.oper2); }
    void operator()(const binop<op_or >& b) const { print(" | ", b.oper1, b.oper2); }
    void operator()(const binop<op_xor>& b) const { print(" ^ ", b.oper1, b.oper2); }

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


class ExpressionFactory {
public:
    expr_t _and_head(expr_t exp1, std::vector<expr_t> exp2) {
        if(exp2.empty()) {
            return exp1;
        }
        expr_t exp = _and_head(exp2.front(), std::vector<expr_t>(exp2.begin() + 1, exp2.end()));
        return binop<op_and>(exp1, exp);
    }
    expr_t _and(std::vector<expr_t> exp) { return _and_head(exp.front(), std::vector<expr_t>(exp.begin() + 1, exp.end())); }
    expr_t _and(std::initializer_list<expr_t> exp) { return _and(std::vector<expr_t>(exp)); }

    expr_t _or_head(expr_t exp1, std::vector<expr_t> exp2) {
        if(exp2.empty()) {
            return exp1;
        }
        expr_t exp = _or_head(exp2.front(), std::vector<expr_t>(exp2.begin() + 1, exp2.end()));
        return binop<op_or>(exp1, exp);
    }
    expr_t _or(std::vector<expr_t> exp) { return _or_head(exp.front(), std::vector<expr_t>(exp.begin() + 1, exp.end())); }
    expr_t _or(std::initializer_list<expr_t> exp) { return _or(std::vector<expr_t>(exp)); }

    expr_t _implication(expr_t premise, expr_t conclusion) { return _or({_not(premise), conclusion}); }
    expr_t _not(expr_t exp) { return unop<op_not>(exp); }
};


#endif //PCB_SAT_expr_tESSION_H

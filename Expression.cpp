//
// Created by tomas on 12/15/18.
//

#include "Expression.h"


void expression::boolean::linearize::operator()(expression::boolean::expr_t &ex) const {
    if (ex->type == expr_type_e::_var) {
        operator_var(ex);
    } else if (ex->type == expr_type_e::_op_and) {
        operator_and(ex);
    } else if (ex->type == expr_type_e::_op_or) {
        operator_or(ex);
    } else if (ex->type == expr_type_e::_op_not) {
        operator_not(ex);
    }
}

void expression::boolean::linearize::operator_var(expression::boolean::expr_t &v) const { _os.emplace_back((uint64_t) v.get()); }

void expression::boolean::linearize::operator_and(expression::boolean::expr_t &b) const {
    _os.emplace_back((uint64_t) b.get());
    this->operator()(b->oper1);
    this->operator()(b->oper2);
}

void expression::boolean::linearize::operator_or(expression::boolean::expr_t &b) const {
    _os.emplace_back((uint64_t) b.get());
    this->operator()(b->oper1);
    this->operator()(b->oper2);
}

void expression::boolean::linearize::operator_not(expression::boolean::expr_t &u) const {
    _os.emplace_back((uint64_t) u.get());
    this->operator()(u->oper1);
}

expression::boolean::expr_s::expr_s(const expression::boolean::expr_t oper1, expression::boolean::expr_type_e type) : oper1(oper1), type(type) {}

expression::boolean::expr_s::expr_s(const std::string &var) : var(var), type(expr_type_e::_var) {}

expression::boolean::expr_s::expr_s(const expression::boolean::expr_t oper1, const expression::boolean::expr_t oper2,
                                    expression::boolean::expr_type_e type) : oper1(oper1), oper2(oper2), type(type) {}

expression::boolean::tseitin::tseitin(std::vector<std::vector<CMSat::Lit>> &os, std::vector<uint64_t> &lz) : _os(os), _lz(lz) {
//        _os.new_vars(lz.size());
    for (int i = 0; i < lz.size(); ++i) {
        _index.insert(std::make_pair(lz[i], i));
        _index_reverse.insert(std::make_pair(i, lz[i]));
    }

}

uint64_t expression::boolean::tseitin::indexof(const expression::boolean::expr_t &ex) const {
    return _index.at((uint64_t) ex.get());
}

expression::boolean::expr_t expression::boolean::tseitin::indexof(const uint64_t node) const {
    return expr_t(*((expr_t *) _index_reverse.at(node)));
}

void expression::boolean::tseitin::operator_var(expression::boolean::expr_t &v) const {}

void expression::boolean::tseitin::operator()(expression::boolean::expr_t &ex) const {
    if (ex->type == expr_type_e::_var) {
        operator_var(ex);
    } else if (ex->type == expr_type_e::_op_and) {
        operator_and(ex);
    } else if (ex->type == expr_type_e::_op_or) {
        operator_or(ex);
    } else if (ex->type == expr_type_e::_op_not) {
        operator_not(ex);
    }
}

void expression::boolean::tseitin::operator_and(expression::boolean::expr_t &b) const {
    std::vector<CMSat::Lit> clause;

    uint64_t _a = indexof(b);
    uint64_t _b = indexof(b->oper1);
    uint64_t _c = indexof(b->oper2);

    this->operator()(b->oper1);
    this->operator()(b->oper2);


    clause.emplace_back(CMSat::Lit(_a, true));
    clause.emplace_back(CMSat::Lit(_b, false));
//        _os.add_clause(clause);
    _os.emplace_back(clause);
    clause.clear();

    clause.emplace_back(CMSat::Lit(_a, true));
    clause.emplace_back(CMSat::Lit(_c, false));
//        _os.add_clause(clause);
    _os.emplace_back(clause);
    clause.clear();

    clause.emplace_back(CMSat::Lit(_a, false));
    clause.emplace_back(CMSat::Lit(_b, true));
    clause.emplace_back(CMSat::Lit(_c, true));
//        _os.add_clause(clause);
    _os.emplace_back(clause);
    clause.clear();
}

void expression::boolean::tseitin::operator_or(expression::boolean::expr_t &b) const {
    std::vector<CMSat::Lit> clause;
    uint64_t _a = indexof(b);
    uint64_t _b = indexof(b->oper1);
    uint64_t _c = indexof(b->oper2);

    this->operator()(b->oper1);
    this->operator()(b->oper2);


    clause.emplace_back(CMSat::Lit(_a, false));
    clause.emplace_back(CMSat::Lit(_b, true));
//        _os.add_clause(clause);
    _os.emplace_back(clause);
    clause.clear();

    clause.emplace_back(CMSat::Lit(_a, false));
    clause.emplace_back(CMSat::Lit(_c, true));
//        _os.add_clause(clause);
    _os.emplace_back(clause);
    clause.clear();

    clause.emplace_back(CMSat::Lit(_a, true));
    clause.emplace_back(CMSat::Lit(_b, false));
    clause.emplace_back(CMSat::Lit(_c, false));
//        _os.add_clause(clause);
    _os.emplace_back(clause);
    clause.clear();
}

void expression::boolean::tseitin::operator_not(expression::boolean::expr_t &u) const {
    std::vector<CMSat::Lit> clause;

    uint64_t _a = indexof(u);
    uint64_t _b = indexof(u->oper1);

    this->operator()(u->oper1);

    clause.emplace_back(CMSat::Lit(_a, true));
    clause.emplace_back(CMSat::Lit(_b, true));
//        _os.add_clause(clause);
    _os.emplace_back(clause);
    clause.clear();

    clause.emplace_back(CMSat::Lit(_a, false));
    clause.emplace_back(CMSat::Lit(_b, false));
//        _os.add_clause(clause);
    _os.emplace_back(clause);
    clause.clear();
}

expression::real::expr_s::expr_s(const expression::real::expr_t oper1, expression::real::expr_type_e type) : oper1(oper1), type(type), constant(0) {}

expression::real::expr_s::expr_s(const std::string &var) : var(var), type(expr_type_e::_var),  constant(0) {}

expression::real::expr_s::expr_s(const std::string &var, const uint64_t constant) : var(var), type(expr_type_e::_var),  constant(constant) {}

expression::real::expr_s::expr_s(const uint64_t constant) : type(expr_type_e::_op_constant),  constant(constant) {}

expression::real::expr_s::expr_s(const expression::real::expr_t oper1, const expression::real::expr_t oper2,
                                 expression::real::expr_type_e type) : oper1(oper1), oper2(oper2), type(type),  constant(0) {}

expression::real::generator::generator(const std::vector<uint32_t> &stack, uint32_t max_vars, uint32_t max_constant) : stack(stack), max_vars(max_vars),
                                                                                                                       max_constant(max_constant) {}

expression::real::expr_t expression::real::generator::operator()() {
    uint32_t run = 0;
    do {
        for (uint32_t i = 0; i < stack.size(); i++) {
            int e = stack[i];
            switch (e % (int) expr_type_e::SIZE) {
                case (int) expr_type_e::_var:
                    if(run == 0)
                        operator_var(stack[(++i) % stack.size()] % max_vars);
                    break;
                case (int) expr_type_e::_op_constant:
                    if(run == 0)
                        operator_constant(stack[(++i) % stack.size()] % max_constant);
                    break;
                case (int) expr_type_e::_op_mul:
                    operator_mul();
                    break;
                case (int) expr_type_e::_op_pow:
                    operator_pow();
                    break;
                case (int) expr_type_e::_op_plus:
                    operator_plus();
                    break;
                case (int) expr_type_e::_op_neg:
                    operator_neg();
                    break;
                case (int) expr_type_e::_op_sin:
                    operator_sin();
                    break;
                case (int) expr_type_e::_op_cos:
                    operator_cos();
                    break;
                case (int) expr_type_e::_op_exp:
                    operator_exp();
                    break;
                default:
                    std::abort();
            }
        }
    } while(expr_queue.size() != 1 && run++ < 1);

    if(expr_queue.size() == 1)
        return expr_queue.back();
    return nullptr;
}

void expression::real::generator::operator_var(const uint32_t i) {
    expr_queue.emplace_back(expr_t(new expr_s{"x" + std::to_string(i), i}));
}

void expression::real::generator::operator_constant(const uint64_t i) {
    expr_queue.emplace_back(expr_t(new expr_s{i}));
}

void expression::real::generator::operator_mul() {
    if(expr_queue.size() >= 2) {
        expr_t op1 = expr_queue.back();
        expr_queue.pop_back();
        expr_t op2 = expr_queue.back();
        expr_queue.pop_back();
        expr_queue.emplace_back(op1 * op2);
    }
}

void expression::real::generator::operator_pow() {
    if(expr_queue.size() >= 2) {
        expr_t op1 = expr_queue.back();
        expr_queue.pop_back();
        expr_t op2 = expr_queue.back();
        expr_queue.pop_back();
        expr_queue.emplace_back(op1 ^ op2);
    }
}

void expression::real::generator::operator_plus() {
    if(expr_queue.size() >= 2) {
        expr_t op1 = expr_queue.back();
        expr_queue.pop_back();
        expr_t op2 = expr_queue.back();
        expr_queue.pop_back();
        expr_queue.emplace_back(op1 + op2);
    }
}

void expression::real::generator::operator_neg() {
    if(expr_queue.size() >= 1) {
        expr_t op1 = expr_queue.back();
        expr_queue.pop_back();
        expr_queue.emplace_back(-op1);
    }
}

void expression::real::generator::operator_cos() {
    if(expr_queue.size() >= 1) {
        expr_t op1 = expr_queue.back();
        expr_queue.pop_back();
        expr_queue.emplace_back(expr_t(new expr_s(op1, expr_type_e::_op_cos)));
    }
}

void expression::real::generator::operator_sin() {
    if(expr_queue.size() >= 1) {
        expr_t op1 = expr_queue.back();
        expr_queue.pop_back();
        expr_queue.emplace_back(expr_t(new expr_s(op1, expr_type_e::_op_sin)));
    }
}

void expression::real::generator::operator_exp() {
    if(expr_queue.size() >= 1) {
        expr_t op1 = expr_queue.back();
        expr_queue.pop_back();
        expr_queue.emplace_back(expr_t(new expr_s(op1, expr_type_e::_op_exp)));
    }
}

expression::real::printer::printer() { }

void expression::real::printer::operator()(const expression::real::expr_t &e) {
    switch ((int)e->type) {
        case (int) expr_type_e::_var:
            operator_var(e);
            break;
        case (int) expr_type_e::_op_constant:
            operator_constant(e);
            break;
        case (int) expr_type_e::_op_mul:
            operator_mul(e);
            break;
        case (int) expr_type_e::_op_pow:
            operator_pow(e);
            break;
        case (int) expr_type_e::_op_plus:
            operator_plus(e);
            break;
        case (int) expr_type_e::_op_neg:
            operator_neg(e);
            break;
        case (int) expr_type_e::_op_sin:
            operator_sin(e);
            break;
        case (int) expr_type_e::_op_cos:
            operator_cos(e);
            break;
        case (int) expr_type_e::_op_exp:
            operator_exp(e);
            break;
        default:
            std::abort();
    }
}

void expression::real::printer::operator_var(const expression::real::expr_t &e) {
    std::cout << e->var;
}

void expression::real::printer::operator_constant(const expression::real::expr_t &e) {
    std::cout << e->constant;
}

void expression::real::printer::operator_mul(const expression::real::expr_t &e) {
    std::cout << "(";
    this->operator()(e->oper1);
    std::cout << ")";
    std::cout << "*";
    std::cout << "(";
    this->operator()(e->oper2);
    std::cout << ")";
}

void expression::real::printer::operator_pow(const expression::real::expr_t &e) {
    std::cout << "pow(";
    this->operator()(e->oper1);
    std::cout << ",";
    this->operator()(e->oper2);
    std::cout << ")";
}

void expression::real::printer::operator_plus(const expression::real::expr_t &e) {
    std::cout << "(";
    this->operator()(e->oper1);
    std::cout << ")";
    std::cout << "+";
    std::cout << "(";
    this->operator()(e->oper2);
    std::cout << ")";
}

void expression::real::printer::operator_neg(const expression::real::expr_t &e) {
    std::cout << "-(";
    this->operator()(e->oper1);
    std::cout << ")";
}

void expression::real::printer::operator_sin(const expression::real::expr_t &e) {
    std::cout << "sin(";
    this->operator()(e->oper1);
    std::cout << ")";
}

void expression::real::printer::operator_cos(const expression::real::expr_t &e) {
    std::cout << "cos(";
    this->operator()(e->oper1);
    std::cout << ")";
}

void expression::real::printer::operator_exp(const expression::real::expr_t &e) {
    std::cout << "exp(";
    this->operator()(e->oper1);
    std::cout << ")";
}

expression::real::eval::eval(const std::vector<double> vars) : vars(vars) { }

double expression::real::eval::operator()(const expression::real::expr_t &e) {
    switch ((int)e->type) {
        case (int) expr_type_e::_var:
            return operator_var(e);
        case (int) expr_type_e::_op_constant:
            return operator_constant(e);
        case (int) expr_type_e::_op_mul:
            return operator_mul(e);
        case (int) expr_type_e::_op_pow:
            return operator_pow(e);
        case (int) expr_type_e::_op_plus:
            return operator_plus(e);
        case (int) expr_type_e::_op_neg:
            return operator_neg(e);
        case (int) expr_type_e::_op_sin:
            return operator_sin(e);
        case (int) expr_type_e::_op_cos:
            return operator_cos(e);
        case (int) expr_type_e::_op_exp:
            return operator_exp(e);
        default:
            std::abort();
    }
}

double expression::real::eval::operator_var(const expression::real::expr_t &e) {
    return vars[e->constant];
}

double expression::real::eval::operator_constant(const expression::real::expr_t &e) {
    return e->constant;
}

double expression::real::eval::operator_mul(const expression::real::expr_t &e) {
    return
            this->operator()(e->oper1) *
            this->operator()(e->oper2);
}

double expression::real::eval::operator_pow(const expression::real::expr_t &e) {
    return std::pow(
            std::abs(this->operator()(e->oper1)),
            this->operator()(e->oper2));
}

double expression::real::eval::operator_plus(const expression::real::expr_t &e) {
    return
            this->operator()(e->oper1) +
            this->operator()(e->oper2);
}

double expression::real::eval::operator_neg(const expression::real::expr_t &e) {
    return -this->operator()(e->oper1);
}

double expression::real::eval::operator_sin(const expression::real::expr_t &e) {
    return std::sin(
            this->operator()(e->oper1));
}

double expression::real::eval::operator_cos(const expression::real::expr_t &e) {
    return std::cos(
            this->operator()(e->oper1));
}

double expression::real::eval::operator_exp(const expression::real::expr_t &e) {
    return std::exp(
            this->operator()(e->oper1));
}


expression::boolean::expr_t expression::boolean::operator&&(const expression::boolean::expr_t &a, const expression::boolean::expr_t &b) {
    assert(a.get() != nullptr);
    assert(b.get() != nullptr);
    return expression::boolean::expr_t(new expression::boolean::expr_s(a, b, expression::boolean::expr_type_e::_op_and));
}

expression::boolean::expr_t
expression::boolean::operator||(const expression::boolean::expr_t &a, const expression::boolean::expr_t &b) {
    assert(a.get() != nullptr);
    assert(b.get() != nullptr);
    return expr_t(new expr_s(a, b, expr_type_e::_op_or));
}

expression::boolean::expr_t expression::boolean::operator!(const expression::boolean::expr_t &a) {
    assert(a.get() != nullptr);
    return expr_t(new expr_s(a, expr_type_e::_op_not));
}

expression::boolean::expr_t
expression::boolean::operator>>(const expression::boolean::expr_t &a, const expression::boolean::expr_t &b) {
    assert(a.get() != nullptr);
    assert(b.get() != nullptr);
    return !a || b;
}

expression::boolean::expr_t
expression::boolean::operator*(const expression::boolean::expr_t &a, const expression::boolean::expr_t &b) {
    assert(a.get() != nullptr);
    assert(b.get() != nullptr);
    return (a >> b) && (b >> a);
}

expression::real::expr_t expression::real::operator*(const expression::real::expr_t &a, const expression::real::expr_t &b) {
    assert(a.get() != nullptr);
    assert(b.get() != nullptr);
    return expression::real::expr_t(new expression::real::expr_s(a, b, expression::real::expr_type_e::_op_mul));
}

expression::real::expr_t
expression::real::operator+(const expression::real::expr_t &a, const expression::real::expr_t &b) {
    assert(a.get() != nullptr);
    assert(b.get() != nullptr);
    return expr_t(new expr_s(a, b, expr_type_e::_op_plus));
}

expression::real::expr_t expression::real::operator-(const expression::real::expr_t &a) {
    assert(a.get() != nullptr);
    return expr_t(new expr_s(a, expr_type_e::_op_neg));
}

expression::real::expr_t
expression::real::operator^(const expression::real::expr_t &a, const expression::real::expr_t &b) {
    assert(a.get() != nullptr);
    assert(b.get() != nullptr);
    return expr_t(new expr_s(a, b, expr_type_e::_op_pow));
}

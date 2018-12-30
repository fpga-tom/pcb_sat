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
#include <cmath>

namespace expression {

    namespace boolean {

        struct expr_s;
        typedef std::shared_ptr<expr_s> expr_t;
        enum class expr_type_e {
            _var, _op_not, _op_and, _op_or
        };

        struct expr_s {

            expr_s(const expr_t oper1, expr_type_e type);

            explicit expr_s(const std::string &var);

            expr_s(const expr_t oper1, const expr_t oper2, expr_type_e type);

            expr_t oper1;
            expr_t oper2;
            const std::string var;
            expr_type_e type;

        };


        struct linearize {
            linearize(std::vector<uint64_t> &os) : _os(os) {}

            std::vector<uint64_t> &_os;

            void operator()(expr_t &ex) const;

            void operator_var(expr_t &v) const;

            void operator_and(expr_t &b) const;

            void operator_or(expr_t &b) const;

            void operator_not(expr_t &u) const;
        };

        struct tseitin {
            std::vector<uint64_t> &_lz;
            std::map<uint64_t, uint64_t> _index;
            std::map<uint64_t, uint64_t> _index_reverse;
            std::vector<std::vector<CMSat::Lit>> &_os;
//    CMSat::SATSolver& _os;

            tseitin(std::vector<std::vector<CMSat::Lit>> &os, std::vector<uint64_t> &lz);


            uint64_t indexof(const expr_t &ex) const;

            expr_t indexof(const uint64_t node) const;

            void operator_var(expr_t &v) const;

            void operator()(expr_t &ex) const;

            void operator_and(expr_t &b) const;

            void operator_or(expr_t &b) const;

            void operator_not(expr_t &u) const;
        };


        expr_t operator&&(const expr_t &a, const expr_t &b);

        expr_t operator||(const expr_t &a, const expr_t &b);

        expr_t operator!(const expr_t &a);

        expr_t operator>>(const expr_t &a, const expr_t &b);

        expr_t operator*(const expr_t &a, const expr_t &b);
    }

    namespace real {

        struct expr_s;
        typedef std::shared_ptr<expr_s> expr_t;
        enum class expr_type_e {
            _var, _op_neg, _op_plus, _op_mul, _op_pow, _op_sin, _op_cos, _op_exp, _op_constant, SIZE
        };

        struct expr_s {

            expr_s(const expr_t oper1, expr_type_e type);

            explicit expr_s(const std::string &var);
            explicit expr_s(const std::string &var, const uint64_t constant);
            explicit expr_s(const uint64_t constant);

            expr_s(const expr_t oper1, const expr_t oper2, expr_type_e type);

            expr_t oper1;
            expr_t oper2;
            const std::string var;
            const uint64_t constant;
            expr_type_e type;

        };

        expr_t operator*(const expr_t &a, const expr_t &b);

        expr_t operator+(const expr_t &a, const expr_t &b);

        expr_t operator-(const expr_t &a);

        expr_t operator^(const expr_t &a, const expr_t &b);

        struct generator {
            generator(const std::vector<uint32_t>& stack, uint32_t max_vars, uint32_t max_constant);

            std::vector<uint32_t > stack;
            std::vector<expr_t> expr_queue;
            const uint32_t max_vars;
            const uint32_t max_constant;

            expr_t operator()();

            void operator_var(const uint32_t i);

            void operator_constant(const uint64_t i);

            void operator_mul();

            void operator_pow();

            void operator_plus();

            void operator_neg();

            void operator_cos();
            void operator_sin();

            void operator_exp();
        };

        struct printer {
            printer();

            void operator()(const expr_t& e);

            void operator_var(const expr_t& e);

            void operator_constant(const expr_t& e);

            void operator_mul(const expr_t& e);

            void operator_pow(const expr_t& e);

            void operator_plus(const expr_t& e);

            void operator_neg(const expr_t& e);

            void operator_sin(const expr_t& e);

            void operator_cos(const expr_t& e);

            void operator_exp(const expr_t& e);
        };


        struct eval {
            eval(const std::vector<double> vars);

            const std::vector<double> vars;

            double operator()(const expr_t& e);

            double operator_var(const expr_t& e);

            double operator_constant(const expr_t& e);

            double operator_mul(const expr_t& e);

            double operator_pow(const expr_t& e);

            double operator_plus(const expr_t& e);

            double operator_neg(const expr_t& e);

            double operator_sin(const expr_t& e);

            double operator_cos(const expr_t& e);

            double operator_exp(const expr_t& e);
        };

    }
}



#endif //PCB_SAT_expr_tESSION_H

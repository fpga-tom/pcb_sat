//
// Created by tomas on 12/29/18.
//

#ifndef PCB_SAT_DYNAMIC_H
#define PCB_SAT_DYNAMIC_H

#include <boost/numeric/odeint.hpp>
#include <vector>
#include <cstdint>
#include <cstdlib>
#include <random>
#include "Expression.h"

typedef std::vector< double > state_type;

typedef boost::numeric::odeint::runge_kutta_cash_karp54< state_type > error_stepper_type;
typedef boost::numeric::odeint::controlled_runge_kutta< error_stepper_type > controlled_stepper_type;


class Dynamic {

    std::vector<std::vector<uint32_t >> ge;
    std::vector<uint32_t > current;
    double Temp;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution;

    uint64_t system;
    uint64_t gen;

    controlled_stepper_type controlled_stepper;
    boost::numeric::odeint::runge_kutta4< state_type > stepper;

    std::vector<expression::real::expr_t> generate_exp(std::vector<uint32_t >&);
    double err_old;
    double err_best;

    struct push_back_state_and_time
    {
        std::vector< state_type >& m_states;
        std::vector< double >& m_times;

        push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
                : m_states( states ) , m_times( times ) { }

        void operator()( const state_type &x , double t )
        {
            m_states.push_back( x );
            m_times.push_back( t );
        }
    };

    double error(std::vector< state_type >& m_states,
    std::vector< double >& m_times);

public:

    Dynamic(uint64_t gen, uint64_t system, double Temp) : system(system), gen(gen), Temp(Temp), distribution(0,1), err_old(0), current(gen), err_best(1e10) {
        for(int j = 0; j < system;j++) {
            std::vector<uint32_t> v;
            ge.emplace_back(v);
            for (int i = 0; i < gen; i++) {
                ge[j].emplace_back(rand());
            }
        }
    }

    void step(uint64_t i);

    struct osc {



        const std::vector<expression::real::expr_t> expr;
        osc(const std::vector<expression::real::expr_t>& expr);

        void operator()(const state_type &x, state_type &dxdt, const double /* t */);
    };

};


#endif //PCB_SAT_DYNAMIC_H

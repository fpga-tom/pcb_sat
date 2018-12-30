//
// Created by tomas on 12/29/18.
//

#include "Dynamic.h"


void Dynamic::step(uint64_t i) {

    std::copy(ge[0].begin(), ge[0].end(), current.begin());
    current[rand() % gen] = rand();
    expression::real::generator generat(current, system, 10);
    expression::real::expr_t e = generat();
    while(e == nullptr) {
        std::copy(ge[0].begin(), ge[0].end(), current.begin());
        current[rand() % gen] = rand();
        expression::real::generator generat(current, system, 10);
        e = generat();
    }

    osc osc1(generate_exp(current));
    std::vector<double> x({0});
    std::vector<state_type> x_vec;
    std::vector<double> times;
    integrate_adaptive( stepper , osc1, x , 0.0 , 12.0 , 0.01, push_back_state_and_time( x_vec , times ) );
    double err = error(x_vec, times);
    err_old = err;
    Temp = Temp*0.9995;
    std::cout << Temp << "/" << err_best << "/" << err << "/" << 1./(1+std::exp(-(err_best - err)/(Temp))) << std::endl;
    if(err < err_best || 1./(1+std::exp(-(err_best - err)/(Temp))) > distribution(generator)) {
        std::copy(current.begin(), current.end(), ge[0].begin());
        err_best = err;
        expression::real::printer p;
        p(e);
        std::cout << std::endl;
    }


}

std::vector<expression::real::expr_t> Dynamic::generate_exp(std::vector<uint32_t >& v) {
    std::vector<expression::real::expr_t> ret;
    for(int i = 0; i < system; i++) {
        expression::real::generator generator(v, system, 10);
        expression::real::expr_t e = generator();

        ret.emplace_back(e);
    }
    return ret;
}

double Dynamic::error(std::vector<state_type> &m_states, std::vector<double> &m_times) {
    double ret = 0;
    double ret1 = 0;
    for(int i = 0; i < m_times.size(); i++) {
        ret += m_states[i][0] * std::cos(m_times[i]);
        ret1 += std::abs(m_states[i][0] - std::cos(m_times[i]));
    }
    return  + ret1/m_times.size();
}

Dynamic::osc::osc(const std::vector<expression::real::expr_t> &expr) : expr(expr) {}

void Dynamic::osc::operator()(const state_type &x, state_type &dxdt, const double) {
    expression::real::eval ev(x);
    for(int i = 0; i < expr.size(); i++) {
        double v = ev(expr[i]);
        if(v > 1e10) {
            v = 1e10;
        }
        if(v < -1e10) {
            v = -1e10;
        }
        if(isnan(v)) {
            dxdt[i] = 1e10;
        } else {
            dxdt[i] = v;
        }

    }
}

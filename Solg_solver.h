//
// Created by tomas on 12/23/18.
//

#ifndef PCB_SAT_SOLG_SOLVER_H
#define PCB_SAT_SOLG_SOLVER_H

#include <cmath>
#include <cstdint>
#include <iostream>
#include <fstream>

typedef double solver_t;

class Solg_solver {
public:
    const solver_t C = 1e-5;
    const solver_t R = 1;
    const solver_t alpha = 60;
    const solver_t R_on = .01;
    const solver_t R_off = 1;

    enum class voltage_t {v1, v2, v3, V1_R, V2_R, V3_R, V1_M, V2_M, V3_M, SIZE};
    enum class memristor_t {x1, x2, x3, x4, x5, SIZE};
    enum class gate_t {_and, _or};

    const gate_t gate;

    Solg_solver(const gate_t& gate) : gate(gate) {}

    solver_t theta(const solver_t y) {
        if(y > (solver_t )1) return 1;
        if(y < (solver_t )0) return 0;
        solver_t ret = 3*std::pow(y,(solver_t )2) -  2*std::pow(y,(solver_t )3);
        return ret;

    }
    solver_t h(const solver_t x, const solver_t v) {
        return (1-std::exp((solver_t )-2*x)) * theta(v/((solver_t )(2*.1))) + ((1-std::exp(-2*(1-x))) * theta(-v/(solver_t )(2*.1)));
    }
    solver_t g(const solver_t x) {
        return ((solver_t )1.)/((R_off - R_on)*x + R_on);
    }
    solver_t v(const memristor_t j, const solver_t voltage[]) {
        solver_t f = (gate == gate_t ::_and ? (solver_t )1. : (solver_t )-1.);
        switch (j) {
            case memristor_t ::x1:
                return (voltage[(int)voltage_t ::V1_M] - voltage[(int)voltage_t::v1]) * f;
            case memristor_t ::x2:
                return ( voltage[(int)voltage_t ::V2_M] - voltage[(int)voltage_t::v2] ) * f;
            case memristor_t ::x3:
                return ( voltage[(int)voltage_t ::v3] - voltage[(int)voltage_t::V3_M] ) * f;
            case memristor_t ::x4:
                return (voltage[(int)voltage_t ::v1] - voltage[(int)voltage_t::v3]) * f;
            case memristor_t ::x5:
                return (voltage[(int)voltage_t ::v2] - voltage[(int)voltage_t::v3]) * f;
            default:
                std::cerr << "unknown" << std::endl;
                std::abort();
        }
    }
    solver_t dxdt(const memristor_t j, const solver_t voltage[], const solver_t mem_state[]) {
        solver_t vm = v(j, voltage);
        return -alpha*h(mem_state[(int) j], vm) * g(mem_state[(int) j])*vm;

    }

    solver_t dv2dt(const solver_t voltage[], const solver_t mem_state[]) {
        solver_t f = (gate == gate_t ::_and ? (solver_t )1. : (solver_t )-1.);
        return ((solver_t )1./C)*
                                (
                                        ((solver_t)-2) * (v(memristor_t::x5, voltage) * g(mem_state[(int)memristor_t::x5]) +
                                        -v(memristor_t::x2, voltage) * g(mem_state[(int)memristor_t::x2]) +
                                        f*(voltage[(int)voltage_t::v2] - voltage[(int)voltage_t::V2_R])/R) +

                                        v(memristor_t::x4, voltage) * g(mem_state[(int)memristor_t::x4]) +
                                        f*(voltage[(int)voltage_t::V3_R] - voltage[(int)voltage_t::v3])/R +
                                        v(memristor_t::x5, voltage) * g(mem_state[(int)memristor_t::x5]) +
                                        -v(memristor_t::x3, voltage) * g(mem_state[(int)memristor_t::x3])
                                );
    }
    solver_t dv3dt(const solver_t voltage[], const solver_t mem_state[]) {
        solver_t f = (gate == gate_t ::_and ? (solver_t )1. : (solver_t )-1.);
        return ((solver_t ).5/C)*
               (
                       ((solver_t)-3.) * (v(memristor_t::x5, voltage) * g(mem_state[(int)memristor_t::x5]) +
                                  -v(memristor_t::x2, voltage) * g(mem_state[(int)memristor_t::x2]) +
                                  f*(voltage[(int)voltage_t::v2] - voltage[(int)voltage_t::V2_R])/R) +

                                  ((solver_t)2.)     *  (v(memristor_t::x4, voltage) * g(mem_state[(int)memristor_t::x4]) +
                                    f*(voltage[(int)voltage_t::V3_R] - voltage[(int)voltage_t::v3])/R +
                                    v(memristor_t::x5, voltage) * g(mem_state[(int)memristor_t::x5]) +
                                    -v(memristor_t::x3, voltage) * g(mem_state[(int)memristor_t::x3]))
                                    );

    }

    void voltage_update(solver_t voltage[]) {
        solver_t f = (gate == gate_t ::_and ? (solver_t )1. : (solver_t )-1.);
        voltage[(int)voltage_t::V1_M] = 0 * voltage[(int)voltage_t::v1] -1*voltage[(int)voltage_t::v2] + 1*voltage[(int)voltage_t::v3] + 1 * f;
        voltage[(int)voltage_t::V1_R] = 3 * voltage[(int)voltage_t::v1] + 1* voltage[(int)voltage_t::v2] -2*voltage[(int)voltage_t::v3] - 1 * f;
        voltage[(int)voltage_t::V2_M] = -1 * voltage[(int)voltage_t::v1] + 0*voltage[(int)voltage_t::v2] + 1*voltage[(int)voltage_t::v3] +1 *f;
        voltage[(int)voltage_t::V2_R] = 1 * voltage[(int)voltage_t::v1] + 3* voltage[(int)voltage_t::v2] -2*voltage[(int)voltage_t::v3] - 1*f;
        voltage[(int)voltage_t::V3_M] = 2 * voltage[(int)voltage_t::v1] + 2*voltage[(int)voltage_t::v2]  -1*voltage[(int)voltage_t::v3] - 2 *f ;
        voltage[(int)voltage_t::V3_R] = -3 * voltage[(int)voltage_t::v1]-3*voltage[(int)voltage_t::v2] + 5*voltage[(int)voltage_t::v3] + 2 *f;
    }

    void rk4() {

        solver_t voltage[(int)voltage_t::SIZE] = {0,};
        solver_t mem_state[(int)memristor_t ::SIZE] = {0.2,.1,.3,.4,.121};
//        std::fill(mem_state, &mem_state[(int)memristor_t::SIZE], .1);
        std::fill(voltage, &voltage[(int)voltage_t ::SIZE], 0);
        solver_t step = 1e-6;

        voltage[(int)voltage_t::v1] = -1.0;
        voltage[(int)voltage_t::v2] = 1.0;
        voltage[(int)voltage_t::v3] = 1.;
//        mem_state[(int)memristor_t::x3] = .75;
        std::ofstream v1_out("/tmp/v1.txt");
        std::ofstream v2_out("/tmp/v2.txt");
        std::ofstream v3_out("/tmp/v3.txt");

        for(solver_t t = 0; t < .004; t+=step) {
            voltage_update(voltage);
            solver_t mem_state_new_k1[(int)memristor_t::SIZE] = {0,};
            solver_t mem_state_step[(int)memristor_t::SIZE] = {0,};
            solver_t voltage_k1[(int)voltage_t ::SIZE] = {0,};
            solver_t voltage_step[(int)voltage_t ::SIZE] = {0,};

            // rk4 step k1
            for(int i = 0 ; i < (int)memristor_t::SIZE; i++) {
                mem_state_new_k1[i] = dxdt((memristor_t) i, voltage, mem_state) * step;
            }
            for(int i = 0 ; i < (int)memristor_t::SIZE; i++) {
                mem_state_step[i] = mem_state_new_k1[i] / 2 + mem_state[i];

            }
            std::copy(voltage, &voltage[(int)voltage_t::SIZE], voltage_k1);

            voltage_k1[(int)voltage_t::v2] = dv2dt(voltage, mem_state) * step;
            voltage_k1[(int)voltage_t::v3] = dv3dt(voltage, mem_state) * step;
            std::copy(voltage, &voltage[(int)voltage_t::SIZE], voltage_step);

            voltage_step[(int)voltage_t::v2] = voltage_k1[(int)voltage_t::v2] / 2 + voltage[(int)voltage_t::v2];
            voltage_step[(int)voltage_t::v3] = voltage_k1[(int)voltage_t::v3] / 2 + voltage[(int)voltage_t::v3];

            // rk4 step k2
            solver_t mem_state_new_k2[(int)memristor_t::SIZE] = {0,};
            solver_t voltage_k2[(int)voltage_t ::SIZE] = {0,};
            voltage_update(voltage_step);
            for(int i = 0 ; i < (int)memristor_t::SIZE; i++) {
                mem_state_new_k2[i] = dxdt((memristor_t) i, voltage_step, mem_state_step) * step;
            }
            for(int i = 0 ; i < (int)memristor_t::SIZE; i++) {
                mem_state_step[i] = mem_state_new_k2[i] / 2 + mem_state[i];

            }


            std::copy(voltage, &voltage[(int)voltage_t::SIZE], voltage_k2);
            voltage_k2[(int)voltage_t::v2] = dv2dt(voltage_step, mem_state_step) * step;
            voltage_k2[(int)voltage_t::v3] = dv3dt(voltage_step, mem_state_step) * step;
            std::copy(voltage, &voltage[(int)voltage_t::SIZE], voltage_step);

            voltage_step[(int)voltage_t::v2] = voltage_k2[(int)voltage_t::v2] / 2 + voltage[(int)voltage_t::v2];
            voltage_step[(int)voltage_t::v3] = voltage_k2[(int)voltage_t::v3] / 2 + voltage[(int)voltage_t::v3];

            // rk4 step k3
            solver_t mem_state_new_k3[(int)memristor_t::SIZE] = {0,};
            solver_t voltage_k3[(int)voltage_t ::SIZE] = {0,};
            voltage_update(voltage_step);
            for(int i = 0 ; i < (int)memristor_t::SIZE; i++) {
                mem_state_new_k3[i] = dxdt((memristor_t) i, voltage_step, mem_state_step) * step;
            }
            for(int i = 0 ; i < (int)memristor_t::SIZE; i++) {
                mem_state_step[i] = mem_state_new_k3[i] + mem_state[i];

            }

            std::copy(voltage, &voltage[(int)voltage_t::SIZE], voltage_k3);
            voltage_k3[(int)voltage_t::v2] = dv2dt(voltage_step, mem_state_step) * step;
            voltage_k3[(int)voltage_t::v3] = dv3dt(voltage_step, mem_state_step) * step;
            std::copy(voltage, &voltage[(int)voltage_t::SIZE], voltage_step);

            voltage_step[(int)voltage_t::v2] = voltage_k3[(int)voltage_t::v2] + voltage[(int)voltage_t::v2];
            voltage_step[(int)voltage_t::v3] = voltage_k3[(int)voltage_t::v3] + voltage[(int)voltage_t::v3];

            // rk4 step k4
            solver_t mem_state_new_k4[(int)memristor_t::SIZE] = {0,};
            solver_t voltage_k4[(int)voltage_t ::SIZE] = {0,};
            voltage_update(voltage_step);
            for(int i = 0 ; i < (int)memristor_t::SIZE; i++) {
                mem_state_new_k4[i] = dxdt((memristor_t )i, voltage_step, mem_state_step) * step;

            }

            std::copy(voltage, &voltage[(int)voltage_t::SIZE], voltage_k4);
            voltage_k4[(int)voltage_t::v2] = dv2dt(voltage_step, mem_state_step) * step;
            voltage_k4[(int)voltage_t::v3] = dv3dt(voltage_step, mem_state_step) * step;
            std::copy(voltage, &voltage[(int)voltage_t::SIZE], voltage_step);

            for(int i = 0 ; i < (int)memristor_t::SIZE; i++) {
                mem_state[i] += (mem_state_new_k1[i] + 2 * mem_state_new_k2[i] + 2* mem_state_new_k3[i] + mem_state_new_k4[i]) / 6.;
            }

            for(int i = 1 ; i <= (int)voltage_t::v3; i++) {
                voltage[i] += (voltage_k1[i] + 2 * voltage_k2[i] + 2* voltage_k3[i] + voltage_k4[i]) / 6.;
            }



            v1_out << t << "\t" << voltage[(int)voltage_t::v1] << std::endl;
            v2_out << t << "\t" << voltage[(int)voltage_t::v2] << std::endl;
            v3_out << t << "\t" << voltage[(int)voltage_t::v3] << std::endl;

        }

    }

};


#endif //PCB_SAT_SOLG_SOLVER_H

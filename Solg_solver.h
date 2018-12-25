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

//#define V1
//#define V2
#define V3

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

    solver_t theta(const solver_t y);
    solver_t h(const solver_t x, const solver_t v);
    solver_t g(const solver_t x);
    solver_t v(const memristor_t j, const solver_t voltage[]);
    solver_t dxdt(const memristor_t j, const solver_t voltage[], const solver_t mem_state[]);
    solver_t mem_1(const solver_t voltage[], const solver_t mem_state[]);
    solver_t mem_2(const solver_t voltage[], const solver_t mem_state[]);
    solver_t mem_3(const solver_t voltage[], const solver_t mem_state[]);
    solver_t dv1dt(const solver_t voltage[], const solver_t mem_state[]);
    solver_t dv2dt(const solver_t voltage[], const solver_t mem_state[]);
    solver_t dv3dt(const solver_t voltage[], const solver_t mem_state[]);
    void voltage_update(solver_t voltage[]);
    void rk4();

};


#endif //PCB_SAT_SOLG_SOLVER_H

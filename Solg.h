//
// Created by tomas on 12/26/18.
//

#ifndef PCB_SAT_SOLG_H
#define PCB_SAT_SOLG_H

#include <cstdint>
#include <cstring>
#include <random>



//#define V1
//#define V2
#define V3

namespace solg {

    typedef float solg_t;

    namespace gate2 {
        enum class component : uint8_t {
            terminal, memristor, SIZE
        };

        enum class logic : uint8_t { _and, _or, SIZE };

        namespace state {
            enum class terminal : uint8_t {
                v1, v2, v3, V1_R, V2_R, V3_R, V1_M, V2_M, V3_M, i1, i2, i3, SIZE
            };
            enum class memristor : uint8_t {
                x1, x2, x3, x4, x5, SIZE
            };
        }

        typedef struct {
            solg_t term_state[(int)gate2::state::terminal::SIZE];
            solg_t mem_state[(int)gate2::state::memristor::SIZE];
        } state_t;


        class Gate {

            const solg_t step = 1e-7;
            const solg_t C = 1e-5;
            const solg_t R = 1;
            const solg_t alpha = 60;
            const solg_t R_on = .01;
            const solg_t R_off = 1;

            const gate2::logic gate;
        public:
            gate2::state_t state, state_prev;//, k1, k2, k3, k4, k1_step, k2_step, k3_step, k4_step;
            solg_t I[2];

            Gate(const gate2::logic& gate) : gate(gate), I{0,0} {
                memset(&state, 0, sizeof(state));
                memset(&state_prev, 0, sizeof(state_prev));
//                memset(&k1, 0, sizeof(k1));
//                memset(&k2, 0, sizeof(k2));
//                memset(&k3, 0, sizeof(k3));
//                memset(&k4, 0, sizeof(k4));
//
//                memset(&k1_step, 0, sizeof(k1_step));
//                memset(&k2_step, 0, sizeof(k2_step));
//                memset(&k3_step, 0, sizeof(k3_step));
//                memset(&k4_step, 0, sizeof(k4_step));

                std::uniform_real_distribution<double> unif(0,1);
                std::default_random_engine re;
                for(int i = 0;i < (int)gate2::state::memristor::SIZE; i++) {
                    state.mem_state[i] = unif(re);
                }
                state.term_state[(int)gate2::state::terminal::v3] = 1.0;
            }

            solg_t theta(solg_t y);

            solg_t h(solg_t x, solg_t v);

            solg_t g(gate2::state::memristor j, const gate2::state_t& s);

            solg_t v(gate2::state::memristor j, const gate2::state_t& s);

            solg_t dxdt(gate2::state::memristor j, const gate2::state_t& s);

            solg_t mem_1(const gate2::state_t& s);

            solg_t mem_2(const gate2::state_t& s);

            solg_t mem_3(const gate2::state_t& s);

            solg_t mem_1_i_off(const gate2::state_t& s);

            solg_t mem_2_i_off(const gate2::state_t& s);

            solg_t mem_3_i_off(const gate2::state_t& s);

            solg_t dv1dt(const gate2::state_t& s);

            solg_t dv2dt(const gate2::state_t& s);

            solg_t dv3dt(const gate2::state_t& s);

            void voltage_update(const gate2::state_t& s1, gate2::state_t& s2);

            solg_t current_i1();
            solg_t current_i2();

            void current_update() {
                I[0] = current_i1();
                I[1] = current_i2();
            }

            void rk4_1(const solg_t i1, const solg_t i2);
//            void rk4_2();
//            void rk4_3();
//            void rk4_4();

            void rk4(const solg_t i1, const solg_t i2) {
                rk4_1(i1, i2);
//                rk4_2();
//                rk4_3();
//                rk4_4();
            }
        };
    }


    namespace gate3 {
        enum class component : uint8_t {
            terminal, memristor, SIZE
        };

        enum class logic : uint8_t { _and, _or, SIZE };

        namespace state {
//            enum class terminal : uint8_t {
//                v1, v2, v3, v4, V1_R, V2_R, V3_R, V4_R, V1_M, V2_M, V3_M, V4_M, i1, i2, i3, i4, SIZE
//            };
            enum class terminal : uint8_t {
                v1, v2, v3, v4, V1_R, V2_R, V3_R, V1_M, V2_M, V3_M, i1, i2, i3, SIZE
            };
            enum class memristor : uint8_t {
                x1, x2, x3, x4, x5, x6, SIZE
            };
        }

        typedef struct {
            solg_t term_state[(int)gate3::state::terminal::SIZE];
            solg_t mem_state[(int)gate3::state::memristor::SIZE];
            solg_t m1, s2, s3;
        } state_t;


        class Gate {

            const solg_t step = 1e-7;
            const solg_t C = 1e-5;
            const solg_t R = 1;
            const solg_t alpha = 60;
            const solg_t R_on = .01;
            const solg_t R_off = 1;

            const gate3::logic gate;
        public:
            gate3::state_t state, state_prev;//, k1, k2, k3, k4, k1_step, k2_step, k3_step, k4_step;
            solg_t I[3];

            Gate(const gate3::logic& gate) : gate(gate), I{0,0,0} {
                memset(&state, 0, sizeof(state));
                memset(&state_prev, 0, sizeof(state_prev));
//                memset(&k1, 0, sizeof(k1));
//                memset(&k2, 0, sizeof(k2));
//                memset(&k3, 0, sizeof(k3));
//                memset(&k4, 0, sizeof(k4));
//
//                memset(&k1_step, 0, sizeof(k1_step));
//                memset(&k2_step, 0, sizeof(k2_step));
//                memset(&k3_step, 0, sizeof(k3_step));
//                memset(&k4_step, 0, sizeof(k4_step));

                std::uniform_real_distribution<double> unif(0,1);
                std::default_random_engine re;
                for(int i = 0;i < (int)gate3::state::memristor::SIZE; i++) {
                    state.mem_state[i] = unif(re);
                }
                state.term_state[(int)gate3::state::terminal::v4] = 1.0;
            }

            solg_t theta(solg_t y);

            solg_t h(solg_t x, solg_t v);

            solg_t g(gate3::state::memristor j, const gate3::state_t& s);

            solg_t v(gate3::state::memristor j, const gate3::state_t& s);

            solg_t dxdt(gate3::state::memristor j, const gate3::state_t& s);

            solg_t mem_1(const gate3::state_t& s);

            solg_t mem_2(const gate3::state_t& s);

            solg_t mem_3(const gate3::state_t& s);

            solg_t mem_1_i_off(const gate3::state_t& s);

            solg_t mem_2_i_off(const gate3::state_t& s);

            solg_t mem_3_i_off(const gate3::state_t& s);

            solg_t dv1dt(const gate3::state_t& s);

            solg_t dv2dt(const gate3::state_t& s);

            solg_t dv3dt(const gate3::state_t& s);

            void dv123dt(gate3::state_t& s);

            void voltage_update(const gate3::state_t& s1, gate3::state_t& s2);

            solg_t current_i1();
            solg_t current_i2();
            solg_t current_i3();

            void current_update() {
                I[0] = current_i1();
                I[1] = current_i2();
                I[2] = current_i3();
            }

            void rk4_1(const solg_t i1, const solg_t i2, const solg_t i3);
//            void rk4_2();
//            void rk4_3();
//            void rk4_4();

            void rk4(const solg_t i1, const solg_t i2, const solg_t i3) {
                rk4_1(i1, i2, i3);
//                rk4_2();
//                rk4_3();
//                rk4_4();
            }
        };
    }

}


#endif //PCB_SAT_SOLG_H

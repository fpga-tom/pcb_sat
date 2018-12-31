//
// Created by tomas on 12/26/18.
//

#include <cmath>
#include <cstdlib>
#include "Solg.h"

namespace solg::gate2 {

    using namespace solg::gate2::state;

solg_t solg::gate2::Gate::theta(const solg_t y) {
    if(y > (solg_t )1) return 1;
    if(y < (solg_t )0) return 0;
    solg_t ret = 3*std::pow(y,(solg_t )2) -  2*std::pow(y,(solg_t )3);
    return ret;

}
solg_t solg::gate2::Gate::h(const solg_t x, const solg_t v) {
    return (1-std::exp((solg_t )-2*x)) * theta(v/((solg_t )(2*.1))) + ((1-std::exp(-2*(1-x))) * theta(-v/(solg_t )(2*.1)));
}
solg_t solg::gate2::Gate::g(const memristor j,  const gate2::state_t& s) {
    solg_t x = s.mem_state[(int)j];
    return ((solg_t )1.)/((R_off - R_on)*x + R_on);
}
solg_t solg::gate2::Gate::v(const memristor j,  const gate2::state_t& s) {
    switch (j) {
        case memristor ::x1:
            return (s.term_state[(int)terminal ::V1_M] - s.term_state[(int)terminal ::v1]);
        case memristor ::x2:
            return ( s.term_state[(int)terminal  ::V2_M] - s.term_state[(int)terminal ::v2] );
        case memristor ::x3:
            return ( s.term_state[(int)terminal  ::v3] - s.term_state[(int)terminal ::V3_M] );
        case memristor ::x4:
            return (s.term_state[(int)terminal  ::v1] - s.term_state[(int)terminal ::v3]);
        case memristor ::x5:
            return (s.term_state[(int)terminal  ::v2] - s.term_state[(int)terminal ::v3]);
        default:
            std::abort();
    }
}
solg_t solg::gate2::Gate::dxdt(const memristor j, const gate2::state_t& s) {
    solg_t f = (gate == gate2::logic ::_and ? (solg_t )1. : (solg_t )-1.);
    solg_t vm = v(j, s) * f;
    return -alpha*h(s.mem_state[(int) j], vm) * g(j, s)*vm;

}

solg_t solg::gate2::Gate::mem_1_i_off(const gate2::state_t& s) {
    return  v(memristor::x4, s) * g(memristor::x4, s) +
           -v(memristor::x1, s) * g(memristor::x1, s) +
           (s.term_state[(int)terminal ::v1] - s.term_state[(int)terminal ::V1_R])/R;

}

solg_t solg::gate2::Gate::mem_1(const gate2::state_t& s) {
    return -s.term_state[(int)terminal::i1] + mem_1_i_off(s);
}

solg_t solg::gate2::Gate::mem_2_i_off(const gate2::state_t& s) {
    return v(memristor ::x5, s) * g(memristor::x5, s) +
           -v(memristor::x2, s) * g(memristor::x2, s) +
           (s.term_state[(int)terminal ::v2] - s.term_state[(int)terminal ::V2_R])/R;
}

solg_t solg::gate2::Gate::mem_2(const gate2::state_t& s) {
    return -s.term_state[(int)terminal::i2] + mem_2_i_off(s);
}

solg_t solg::gate2::Gate::mem_3_i_off(const gate2::state_t& s) {
    return v(memristor::x4, s) * g(memristor::x4, s) +
           (s.term_state[(int)terminal ::V3_R] - s.term_state[(int)terminal ::v3])/R +
           v(memristor::x5, s) * g(memristor::x5, s) +
           -v(memristor::x3, s) * g(memristor::x3, s);
}

solg_t solg::gate2::Gate::mem_3(const gate2::state_t& s) {
    return -s.term_state[(int)terminal::i3] + mem_3_i_off(s);
}


solg_t solg::gate2::Gate::dv1dt(const gate2::state_t& s) {
#if defined V1 || defined V2
    return ((solg_t )1./C)*
           (
                   ((solg_t)-2) * (mem_1(s)) +

                   mem_3(s)
           );
#else

//    return ((solg_t )1./C/3)*
//           (
//                   ((solg_t)-2) * mem_1(s) +
//
//                   mem_2(s)
//           );
    return -mem_1(s)/3;
#endif
}

solg_t solg::gate2::Gate::dv2dt(const gate2::state_t& s) {
#if defined V1 || defined V2
    return ((solg_t )1./C)*
           (
                   ((solg_t)-2 * mem_2(s)) +

                   mem_3(s)
           );
#else
//    return ((solg_t )1./C/3)*
//           (
//                   ((solg_t)-2) * mem_2(s) +
//
//                   mem_1(s)
//           );
    return -mem_2(s)/3;
#endif
}
solg_t solg::gate2::Gate::dv3dt(const gate2::state_t& s) {
#ifndef V3
    return ((solg_t ).5/C)*
           (
                   ((solg_t)-3.) * (
#if defined V1
                           mem_2(s)) +
#elif defined V2
                           mem_1(s)) +
#endif

                           ((solg_t)2.)     *  (mem_3(s))
                   );
#else
    return 0;
#endif

}

//void solg::gate2::Gate::voltage_update(const gate2::state_t& s1, gate2::state_t& s2) {
//    s2.term_state[(int)terminal ::V1_M] =  0 * s1.term_state[(int)terminal ::v1] - 1 * s1.term_state[(int)terminal ::v2] + 1 * s1.term_state[(int)terminal ::v3] - 1;// * f;
//    s2.term_state[(int)terminal ::V1_R] =  3 * s1.term_state[(int)terminal ::v1] + 1 * s1.term_state[(int)terminal ::v2] - 2 * s1.term_state[(int)terminal ::v3] + 1;// * f;
//    s2.term_state[(int)terminal ::V2_M] = -1 * s1.term_state[(int)terminal ::v1] + 0 * s1.term_state[(int)terminal ::v2] + 1 * s1.term_state[(int)terminal ::v3] - 1;// * f;
//    s2.term_state[(int)terminal ::V2_R] =  1 * s1.term_state[(int)terminal ::v1] + 3 * s1.term_state[(int)terminal ::v2] - 2 * s1.term_state[(int)terminal ::v3] + 1;// * f;
//    s2.term_state[(int)terminal ::V3_M] =  2 * s1.term_state[(int)terminal ::v1] + 2 * s1.term_state[(int)terminal ::v2] - 1 * s1.term_state[(int)terminal ::v3] + 2;// * f;
//    s2.term_state[(int)terminal ::V3_R] = -3 * s1.term_state[(int)terminal ::v1] - 3 * s1.term_state[(int)terminal ::v2] + 5 * s1.term_state[(int)terminal ::v3] - 2;// * f;
//}

    void solg::gate2::Gate::voltage_update(const gate2::state_t& s1, gate2::state_t& s2) {
        s2.term_state[(int)terminal ::V1_M] =  -1 * s1.term_state[(int)terminal ::v1];// - 1;
        s2.term_state[(int)terminal ::V1_R] =  -1 * s1.term_state[(int)terminal ::v1];// - 1;
        s2.term_state[(int)terminal ::V2_M] =  -1 * s1.term_state[(int)terminal ::v2];// - 1;
        s2.term_state[(int)terminal ::V2_R] =  -1 * s1.term_state[(int)terminal ::v2];// - 1;
        s2.term_state[(int)terminal ::V3_M] =  -1 * s1.term_state[(int)terminal ::v3];// - 1;
        s2.term_state[(int)terminal ::V3_R] =  -1 * s1.term_state[(int)terminal ::v3];// - 1;
    }


void solg::gate2::Gate::rk4_1(const solg_t i1, const solg_t i2) {

    state_t k1, k2, k3, k4, k1_step, k2_step, k3_step, k4_step;

    memset(&k1, 0, sizeof(k1));
    memset(&k2, 0, sizeof(k2));
    memset(&k3, 0, sizeof(k3));
    memset(&k4, 0, sizeof(k4));

    memset(&k1_step, 0, sizeof(k1_step));
    memset(&k2_step, 0, sizeof(k2_step));
    memset(&k3_step, 0, sizeof(k3_step));
    memset(&k4_step, 0, sizeof(k4_step));

    state_prev = state;
    state.term_state[(int)terminal::i1] = i1;
    state.term_state[(int)terminal::i2] = i2;
    voltage_update(state, state);

    // rk4 step k1
    k1 = state;
    for (int i = 0; i < (int) memristor::SIZE; i++) {
        k1.mem_state[i] = dxdt((memristor) i, state) * step;
    }
    for (int i = 0; i < (int) memristor::SIZE; i++) {
        k1_step.mem_state[i] = k1.mem_state[i] / 2 + state.mem_state[i];

    }

#if defined(V2) || defined(V3)
    k1.term_state[(int) terminal::v1] = dv1dt(state) * step;
#endif
#if defined(V1) || defined(V3)
    k1.term_state[(int) terminal::v2] = dv2dt(state) * step;
#endif
#ifndef V3
    k1.term_state[(int)terminal ::v3] = dv3dt(state) * step;
#endif
    k1_step = k1;

#if defined(V2) || defined(V3)
    k1_step.term_state[(int) terminal::v1] =
            k1.term_state[(int) terminal::v1] / 2 + state.term_state[(int) terminal::v1];
#endif
#if defined(V1) || defined(V3)
    k1_step.term_state[(int) terminal::v2] =
            k1.term_state[(int) terminal::v2] / 2 + state.term_state[(int) terminal::v2];
#endif
#ifndef V3
    k1_step.term_state[(int)terminal ::v3] = k1.term_state[(int)terminal ::v3] / 2 + state.term_state[(int)terminal ::v3];
#endif
//}
//
//void solg::gate2::Gate::rk4_2() {

    voltage_update(k1_step, k1_step);
    k2 = k1_step;
    // rk4 step k2
    for (int i = 0; i < (int) memristor::SIZE; i++) {
        k2.mem_state[i] = dxdt((memristor) i, k1_step) * step;
    }
    for (int i = 0; i < (int) memristor::SIZE; i++) {
        k2_step.mem_state[i] = k2.mem_state[i] / 2 + state.mem_state[i];

    }

#if defined(V2) || defined(V3)
    k2.term_state[(int) terminal::v1] = dv1dt(k1_step) * step;
#endif
#if defined(V1) || defined(V3)
    k2.term_state[(int) terminal::v2] = dv2dt(k1_step) * step;
#endif
#ifndef V3
    k2.term_state[(int)terminal ::v3] = dv3dt(k1_step) * step;
#endif
    k2_step = k2;

#if defined(V2) || defined(V3)
    k2_step.term_state[(int) terminal::v1] =
            k2.term_state[(int) terminal::v1] / 2 + state.term_state[(int) terminal::v1];
#endif
#if defined(V1) || defined(V3)
    k2_step.term_state[(int) terminal::v2] =
            k2.term_state[(int) terminal::v2] / 2 + state.term_state[(int) terminal::v2];
#endif
#ifndef V3
    k2_step.term_state[(int)terminal ::v3] = k2.term_state[(int)terminal ::v3] / 2 + state.term_state[(int)terminal ::v3];
#endif
//}
//
//void solg::gate2::Gate::rk4_3() {
    // rk4 step k3

    voltage_update(k2_step, k2_step);
    k3 = k2_step;
    // rk4 step k2

    for (int i = 0; i < (int) memristor::SIZE; i++) {
        k3.mem_state[i] = dxdt((memristor) i, k2_step) * step;
    }
    for (int i = 0; i < (int) memristor::SIZE; i++) {
        k3_step.mem_state[i] = k3.mem_state[i] / 2 + state.mem_state[i];

    }

#if defined(V2) || defined(V3)
    k3.term_state[(int) terminal::v1] = dv1dt(k2_step) * step;
#endif
#if defined(V1) || defined(V3)
    k3.term_state[(int) terminal::v2] = dv2dt(k2_step) * step;
#endif
#ifndef V3
    k3.term_state[(int)terminal ::v3] = dv3dt(k2_step) * step;
#endif
    k3_step = k3;

#if defined(V2) || defined(V3)
    k3_step.term_state[(int) terminal::v1] =
            k3.term_state[(int) terminal::v1] + state.term_state[(int) terminal::v1];
#endif
#if defined(V1) || defined(V3)
    k3_step.term_state[(int) terminal::v2] =
            k3.term_state[(int) terminal::v2] + state.term_state[(int) terminal::v2];
#endif
#ifndef V3
    k3_step.term_state[(int)terminal ::v3] = k3.term_state[(int)terminal ::v3] + state.term_state[(int)terminal ::v3];
#endif

//}
//
//void solg::gate2::Gate::rk4_4() {
        // rk4 step k4

    voltage_update(k3_step, k3_step);
    k4 = k3_step;
    for (int i = 0; i < (int) memristor::SIZE; i++) {
        k4.mem_state[i] = dxdt((memristor) i, k3_step) * step;
    }
    for (int i = 0; i < (int) memristor::SIZE; i++) {
        k4_step.mem_state[i] = k4.mem_state[i]  + state.mem_state[i];

    }


#if defined(V2) || defined(V3)
    k4.term_state[(int) terminal::v1] = dv1dt(k3_step) * step;
#endif
#if defined(V1) || defined(V3)
    k4.term_state[(int) terminal::v2] = dv2dt(k3_step) * step;
#endif
#ifndef V3
    k4.term_state[(int)terminal ::v3] = dv3dt(k3_step) * step;
#endif
    k4_step = k4;

#if defined(V2) || defined(V3)
    k4_step.term_state[(int) terminal::v1] =
            k4.term_state[(int) terminal::v1] + state.term_state[(int) terminal::v1];
#endif
#if defined(V1) || defined(V3)
    k4_step.term_state[(int) terminal::v2] =
            k4.term_state[(int) terminal::v2] + state.term_state[(int) terminal::v2];
#endif
#ifndef V3
    k4_step.term_state[(int)terminal ::v3] = k4.term_state[(int)terminal ::v3] + state.term_state[(int)terminal ::v3];
#endif

#ifndef V3
    voltage_k4[(int)voltage_t::v3] = dv3dt(voltage_step, mem_state_step) * step;
#endif

    for(int i = 0 ; i < (int)memristor::SIZE; i++) {
        state.mem_state[i] += (k1.mem_state[i] + 2 * k2.mem_state[i] + 2* k3.mem_state[i] + k4.mem_state[i]) / 6.;
    }

#ifdef V2
    for(int i = 0 ; i <= (int)voltage_t::v3; i++) {
                if(i == (int)voltage_t::v2) continue;
                voltage[i] += (voltage_k1[i] + 2 * voltage_k2[i] + 2* voltage_k3[i] + voltage_k4[i]) / 6.;
            }
#elif defined(V1)
    for(int i = 1 ; i <= (int)voltage_t::v3; i++) {
            voltage[i] += (voltage_k1[i] + 2 * voltage_k2[i] + 2* voltage_k3[i] + voltage_k4[i]) / 6.;
        }
#endif
#ifdef V3
    for(int i = 0 ; i <= (int)terminal ::v2; i++) {
        state.term_state[i] += (k1.term_state[i] + 2 * k2.term_state[i] + 2* k3.term_state[i] + k4.term_state[i]) / 6.;
    }
#endif


}

solg_t solg::gate2::Gate::current_i1() {
//    return -(C*(-2*(state.term_state[(int)terminal::v1] - state_prev.term_state[(int)terminal::v1])
//    -1*(state.term_state[(int)terminal::v2] - state_prev.term_state[(int)terminal::v2])) - mem_1_i_off(state));
    return (C*((1./3)*(state.term_state[(int)terminal::v1] - state_prev.term_state[(int)terminal::v1]))
                 + mem_1_i_off(state));
}

solg_t solg::gate2::Gate::current_i2() {
//    return -(C*(-1*(state.term_state[(int)terminal::v1] - state_prev.term_state[(int)terminal::v1])
//              -2*(state.term_state[(int)terminal::v2] - state_prev.term_state[(int)terminal::v2])) - mem_2_i_off(state));
    return (C*((1./3)*(state.term_state[(int)terminal::v2] - state_prev.term_state[(int)terminal::v2])) + mem_2_i_off(state));
}

}



namespace solg::gate3 {

    using namespace solg::gate3::state;

    solg_t solg::gate3::Gate::theta(const solg_t y) {
        if(y > (solg_t )1) return 1;
        if(y < (solg_t )0) return 0;
        solg_t ret = 3*std::pow(y,(solg_t )2) -  2*std::pow(y,(solg_t )3);
        return ret;

    }
    solg_t solg::gate3::Gate::h(const solg_t x, const solg_t v) {
        return (1-std::exp((solg_t )-2*x)) * theta(v/((solg_t )(2*.1))) + ((1-std::exp(-2*(1-x))) * theta(-v/(solg_t )(2*.1)));
    }
    solg_t solg::gate3::Gate::g(const memristor j,  const gate3::state_t& s) {
        solg_t x = s.mem_state[(int)j];
        return ((solg_t )1.)/((R_off - R_on)*x + R_on);
    }
    solg_t solg::gate3::Gate::v(const memristor j,  const gate3::state_t& s) {
        switch (j) {
            case memristor ::x1:
                return (s.term_state[(int)terminal ::V1_M] - s.term_state[(int)terminal ::v1]);
            case memristor ::x2:
                return ( s.term_state[(int)terminal  ::V2_M] - s.term_state[(int)terminal ::v2] );
            case memristor ::x3:
                return ( s.term_state[(int)terminal  ::V3_M] - s.term_state[(int)terminal ::v3] );
            case memristor ::x4:
                return (s.term_state[(int)terminal  ::v3] - s.term_state[(int)terminal ::v4]);
            case memristor ::x5:
                return (s.term_state[(int)terminal  ::v1] - s.term_state[(int)terminal ::v4]);
            case memristor ::x6:
                return (s.term_state[(int)terminal  ::v2] - s.term_state[(int)terminal ::v4]);
            default:
                std::abort();
        }
    }
    solg_t solg::gate3::Gate::dxdt(const memristor j, const gate3::state_t& s) {
        solg_t f = (gate == gate3::logic ::_and ? (solg_t )1. : (solg_t )-1.);
        solg_t vm = v(j, s) * f;
        return -alpha*h(s.mem_state[(int) j], vm) * g(j, s)*vm;

    }

    solg_t solg::gate3::Gate::mem_1_i_off(const gate3::state_t& s) {
        return  v(memristor::x5, s) * g(memristor::x5, s) +
                -v(memristor::x1, s) * g(memristor::x1, s) +
                (s.term_state[(int)terminal ::v1] - s.term_state[(int)terminal ::V1_R])/R;

    }

    solg_t solg::gate3::Gate::mem_1(const gate3::state_t& s) {
        return -s.term_state[(int)terminal::i1] + mem_1_i_off(s);
    }

    solg_t solg::gate3::Gate::mem_2_i_off(const gate3::state_t& s) {
        return v(memristor ::x6, s) * g(memristor::x6, s) +
               -v(memristor::x2, s) * g(memristor::x2, s) +
               (s.term_state[(int)terminal ::v2] - s.term_state[(int)terminal ::V2_R])/R;
    }

    solg_t solg::gate3::Gate::mem_2(const gate3::state_t& s) {
        return -s.term_state[(int)terminal::i2] + mem_2_i_off(s);
    }

    solg_t solg::gate3::Gate::mem_3_i_off(const gate3::state_t& s) {
        return v(memristor ::x4, s) * g(memristor::x4, s) +
               -v(memristor::x3, s) * g(memristor::x3, s) +
               (s.term_state[(int)terminal ::v3] - s.term_state[(int)terminal ::V3_R])/R;
    }

    solg_t solg::gate3::Gate::mem_3(const gate3::state_t& s) {
        return -s.term_state[(int)terminal::i3] + mem_3_i_off(s);
    }

    void solg::gate3::Gate::dv123dt(gate3::state_t& s) {
//        solg_t m1 = mem_1(s);
//        solg_t s3 = mem_3(s) + m1;
//        solg_t s2 = mem_2(s) - .5*m1;
//        m1 += s3;
//        s2 += 1.5 * s3;
//        m1 += .4 * s2;
//        s3 -= .8 * s2;
//
//        s.m1 = m1;
//        s.s2 = s2;
//        s.s3 = s3;
    }

    solg_t solg::gate3::Gate::dv1dt(const gate3::state_t& s) {
//        return -.5*s.m1;
        return -mem_1(s)/3;
    }

    solg_t solg::gate3::Gate::dv2dt(const gate3::state_t& s) {
//        return s.s3;
        return -mem_2(s)/3;

    }
    solg_t solg::gate3::Gate::dv3dt(const gate3::state_t& s) {
//        return -.4*s.s2;
        return -mem_3(s)/3;

    }

//    void solg::gate3::Gate::voltage_update(const gate3::state_t& s1, gate3::state_t& s2) {
//        s2.term_state[(int)terminal ::V1_M] =  0 * s1.term_state[(int)terminal ::v1] - 1 * s1.term_state[(int)terminal ::v2] + 1 * s1.term_state[(int)terminal ::v3] - 1 * s1.term_state[(int)terminal::v4];
//        s2.term_state[(int)terminal ::V1_R] =  3 * s1.term_state[(int)terminal ::v1] + 1 * s1.term_state[(int)terminal ::v2] - 2 * s1.term_state[(int)terminal ::v3] + 1 * s1.term_state[(int)terminal::v4];
//        s2.term_state[(int)terminal ::V2_M] = -1 * s1.term_state[(int)terminal ::v1] + 0 * s1.term_state[(int)terminal ::v2] + 1 * s1.term_state[(int)terminal ::v3] - 1 * s1.term_state[(int)terminal::v4];
//        s2.term_state[(int)terminal ::V2_R] =  1 * s1.term_state[(int)terminal ::v1] + 3 * s1.term_state[(int)terminal ::v2] - 2 * s1.term_state[(int)terminal ::v3] + 1 * s1.term_state[(int)terminal::v4];
//        s2.term_state[(int)terminal ::V3_M] =  2 * s1.term_state[(int)terminal ::v1] + 2 * s1.term_state[(int)terminal ::v2] - 1 * s1.term_state[(int)terminal ::v3] + 2 * s1.term_state[(int)terminal::v4];
//        s2.term_state[(int)terminal ::V3_R] = -3 * s1.term_state[(int)terminal ::v1] - 3 * s1.term_state[(int)terminal ::v2] + 5 * s1.term_state[(int)terminal ::v3] - 2 * s1.term_state[(int)terminal::v4];
//    }

    void solg::gate3::Gate::voltage_update(const gate3::state_t& s1, gate3::state_t& s2) {
        s2.term_state[(int)terminal ::V1_M] = -1 * s1.term_state[(int)terminal ::v1];// - 1;
        s2.term_state[(int)terminal ::V1_R] = -1 * s1.term_state[(int)terminal ::v1];// - 1;
        s2.term_state[(int)terminal ::V2_M] = -1 * s1.term_state[(int)terminal ::v2];// - 1;
        s2.term_state[(int)terminal ::V2_R] = -1 * s1.term_state[(int)terminal ::v2];// - 1;
        s2.term_state[(int)terminal ::V3_M] = -1 * s1.term_state[(int)terminal ::v3];// - 1;
        s2.term_state[(int)terminal ::V3_R] = -1 * s1.term_state[(int)terminal ::v3];// - 1;
    }


    void solg::gate3::Gate::rk4_1(const solg_t i1, const solg_t i2, const solg_t i3) {
        state_t k1, k2, k3, k4, k1_step, k2_step, k3_step, k4_step;

        memset(&k1, 0, sizeof(k1));
        memset(&k2, 0, sizeof(k2));
        memset(&k3, 0, sizeof(k3));
        memset(&k4, 0, sizeof(k4));

        memset(&k1_step, 0, sizeof(k1_step));
        memset(&k2_step, 0, sizeof(k2_step));
        memset(&k3_step, 0, sizeof(k3_step));
        memset(&k4_step, 0, sizeof(k4_step));

        state_prev = state;
        state.term_state[(int)terminal::i1] = i1;
        state.term_state[(int)terminal::i2] = i2;
        state.term_state[(int)terminal::i3] = i3;
        voltage_update(state, state);
        k1 = state;
        // rk4 step k1
        for (int i = 0; i < (int) memristor::SIZE; i++) {
            k1.mem_state[i] = dxdt((memristor) i, state) * step;
        }
        for (int i = 0; i < (int) memristor::SIZE; i++) {
            k1_step.mem_state[i] = k1.mem_state[i] / 2 + state.mem_state[i];

        }

        dv123dt(state);
        k1.term_state[(int) terminal::v1] = dv1dt(state) * step;
        k1.term_state[(int) terminal::v2] = dv2dt(state) * step;
        k1.term_state[(int) terminal ::v3] = dv3dt(state) * step;
        k1_step = k1;

        k1_step.term_state[(int) terminal::v1] =
                k1.term_state[(int) terminal::v1] / 2 + state.term_state[(int) terminal::v1];
        k1_step.term_state[(int) terminal::v2] =
                k1.term_state[(int) terminal::v2] / 2 + state.term_state[(int) terminal::v2];
        k1_step.term_state[(int)terminal ::v3] = k1.term_state[(int)terminal ::v3] / 2 + state.term_state[(int)terminal ::v3];
//    }
//
//    void solg::gate3::Gate::rk4_2() {

        voltage_update(k1_step, k1_step);
        k2 = k1_step;
        // rk4 step k2
        for (int i = 0; i < (int) memristor::SIZE; i++) {
            k2.mem_state[i] = dxdt((memristor) i, k1_step) * step;
        }
        for (int i = 0; i < (int) memristor::SIZE; i++) {
            k2_step.mem_state[i] = k2.mem_state[i] / 2 + state.mem_state[i];

        }

        dv123dt(k1_step);
        k2.term_state[(int) terminal::v1] = dv1dt(k1_step) * step;
        k2.term_state[(int) terminal::v2] = dv2dt(k1_step) * step;
        k2.term_state[(int)terminal ::v3] = dv3dt(k1_step) * step;
        k2_step = k2;

        k2_step.term_state[(int) terminal::v1] =
                k2.term_state[(int) terminal::v1] / 2 + state.term_state[(int) terminal::v1];
        k2_step.term_state[(int) terminal::v2] =
                k2.term_state[(int) terminal::v2] / 2 + state.term_state[(int) terminal::v2];
        k2_step.term_state[(int)terminal ::v3] = k2.term_state[(int)terminal ::v3] / 2 + state.term_state[(int)terminal ::v3];
//    }
//
//    void solg::gate3::Gate::rk4_3() {
        // rk4 step k3
        voltage_update(k2_step, k2_step);
        k3 = k2_step;
        // rk4 step k2
        for (int i = 0; i < (int) memristor::SIZE; i++) {
            k3.mem_state[i] = dxdt((memristor) i, k2_step) * step;
        }
        for (int i = 0; i < (int) memristor::SIZE; i++) {
            k3_step.mem_state[i] = k3.mem_state[i] / 2 + state.mem_state[i];

        }

        dv123dt(k2_step);
        k3.term_state[(int) terminal::v1] = dv1dt(k2_step) * step;
        k3.term_state[(int) terminal::v2] = dv2dt(k2_step) * step;
        k3.term_state[(int)terminal ::v3] = dv3dt(k2_step) * step;
        k3_step = k3;

        k3_step.term_state[(int) terminal::v1] =
                k3.term_state[(int) terminal::v1] + state.term_state[(int) terminal::v1];
        k3_step.term_state[(int) terminal::v2] =
                k3.term_state[(int) terminal::v2] + state.term_state[(int) terminal::v2];
        k3_step.term_state[(int)terminal ::v3] = k3.term_state[(int)terminal ::v3] + state.term_state[(int)terminal ::v3];

//    }
//
//    void solg::gate3::Gate::rk4_4() {
        // rk4 step k4

        voltage_update(k3_step, k3_step);
        k4 = k3_step;
        for (int i = 0; i < (int) memristor::SIZE; i++) {
            k4.mem_state[i] = dxdt((memristor) i, k3_step) * step;
        }
        for (int i = 0; i < (int) memristor::SIZE; i++) {
            k4_step.mem_state[i] = k4.mem_state[i]  + state.mem_state[i];

        }


        dv123dt(k3_step);
        k4.term_state[(int) terminal::v1] = dv1dt(k3_step) * step;
        k4.term_state[(int) terminal::v2] = dv2dt(k3_step) * step;
        k4.term_state[(int)terminal ::v3] = dv3dt(k3_step) * step;
        k4_step = k4;

        k4_step.term_state[(int) terminal::v1] =
                k4.term_state[(int) terminal::v1] + state.term_state[(int) terminal::v1];
        k4_step.term_state[(int) terminal::v2] =
                k4.term_state[(int) terminal::v2] + state.term_state[(int) terminal::v2];
        k4_step.term_state[(int)terminal ::v3] = k4.term_state[(int)terminal ::v3] + state.term_state[(int)terminal ::v3];


        for(int i = 0 ; i < (int)memristor::SIZE; i++) {
            state.mem_state[i] += (k1.mem_state[i] + 2 * k2.mem_state[i] + 2* k3.mem_state[i] + k4.mem_state[i]) / 6.;
        }

        for(int i = 0 ; i <= (int)terminal ::v3; i++) {
            state.term_state[i] += (k1.term_state[i] + 2 * k2.term_state[i] + 2* k3.term_state[i] + k4.term_state[i]) / 6.;
        }


    }

    solg_t solg::gate3::Gate::current_i1() {
//        return (C*(-2*(state.term_state[(int)terminal::v1] - state_prev.term_state[(int)terminal::v1])
//                  -1*(state.term_state[(int)terminal::v2] - state_prev.term_state[(int)terminal::v2])
//                  +1*(state.term_state[(int)terminal::v3] - state_prev.term_state[(int)terminal::v3])
//                  ) - mem_1_i_off(state));

        return (C*((1./3)*(state.term_state[(int)terminal::v1] - state_prev.term_state[(int)terminal::v1]))
         + mem_1_i_off(state));
    }

    solg_t solg::gate3::Gate::current_i2() {
//        return (C*(-1*(state.term_state[(int)terminal::v1] - state_prev.term_state[(int)terminal::v1])
//                  -2*(state.term_state[(int)terminal::v2] - state_prev.term_state[(int)terminal::v2])
//                  +1*(state.term_state[(int)terminal::v3] - state_prev.term_state[(int)terminal::v3])
//                  ) - mem_2_i_off(state));
        return (C*((1./3)*(state.term_state[(int)terminal::v2] - state_prev.term_state[(int)terminal::v2]))
                + mem_2_i_off(state));
    }

    solg_t solg::gate3::Gate::current_i3() {
//        return (C*(+2*(state.term_state[(int)terminal::v1] - state_prev.term_state[(int)terminal::v1])
//                  +2*(state.term_state[(int)terminal::v2] - state_prev.term_state[(int)terminal::v2])
//                  -3*(state.term_state[(int)terminal::v3] - state_prev.term_state[(int)terminal::v3])
//                  ) - mem_3_i_off(state));
        return (C*((1./3)*(state.term_state[(int)terminal::v3] - state_prev.term_state[(int)terminal::v3]))
                + mem_3_i_off(state));
    }

}
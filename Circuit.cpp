//
// Created by tomas on 12/26/18.
//

#include "Circuit.h"

namespace solg {

    namespace circuit {

        bool operator==(const edge_t &a, const edge_t &b) {
            return a == b;
        }

        bool operator==(const node_t &a, const node_t &b) {
            return a == b;
        }

        bool operator<(const node_t &a, const node_t &b) {
            return a < b;
        }

        edge_t operator+(const node_t &a, const node_t &b) {
            return edge_t(new edge_s(a, b));
        }

    }
}
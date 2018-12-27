//
// Created by tomas on 12/26/18.
//

#ifndef PCB_SAT_CIRCUIT_H
#define PCB_SAT_CIRCUIT_H

#include <vector>
#include <memory>
#include <map>
#include <variant>
#include <cryptominisat5/cryptominisat.h>
#include "Solg.h"

namespace solg {

    namespace circuit {

        using namespace CMSat;

        enum class component {
            wire, resistor, memristor, inverter, SIZE
        };

        enum class power {
            none, vcc, gnd, VM, SIZE
        };

        struct edge_s;

        struct node_s {
            std::shared_ptr<edge_s> edge;
        };

        typedef std::shared_ptr<node_s> node_t;

        struct edge_s {
            component comp;
            node_t oper1, oper2;
            edge_s(const node_t& oper1, const node_t& oper2) : oper1(oper1), oper2(oper2), comp(component::wire) {}


        };

        typedef std::shared_ptr<edge_s> edge_t;

        bool operator==(const edge_t &a, const edge_t &b);

        edge_t operator+(const node_t &a, const node_t &b);

        typedef std::variant<solg::gate2::Gate, solg::gate3::Gate> gate_t;

        struct Circuit {
            std::vector<node_t> nodes;
            std::vector<edge_t> edges;
            std::vector<std::shared_ptr<gate_t>> gates;
            std::map<node_t, std::vector<std::pair<std::shared_ptr<gate_t>, uint32_t >>> terminals;
            std::map<std::pair<std::shared_ptr<gate_t>, uint32_t >, node_t> terminals_reverse;



            Circuit(const std::vector<uint64_t >& linear, const std::vector<std::vector<CMSat::Lit>>& clauses) {
                std::map<CMSat::Lit , node_t > literal_node;

                // create nodes
                for(const std::vector<CMSat::Lit>& clause : clauses) {
                    for(const CMSat::Lit& lit : clause) {
                        if(literal_node.count(lit) == 0) {
                            nodes.emplace_back(new node_s());
                            literal_node.insert(std::make_pair(lit, nodes.back()));
                        }
                    }
                }
                for(auto it = literal_node.begin(); it != literal_node.end(); ++it) {
                    CMSat::Lit a = it->first;

                    if(literal_node.count(~a) > 0 && !a.sign()) {
                        edge_t edge = it->second + literal_node[~a];
                        edge->comp = component ::inverter;
                        it->second->edge = edge;
                        literal_node[~a]->edge = edge;
                        edges.emplace_back(edge);
                    }
                }

                for(const std::vector<CMSat::Lit>& clause : clauses) {
                    if(clause.size() == 2) {
                        std::shared_ptr<gate_t> g(new gate_t(solg::gate2::Gate(solg::gate2::logic::_or)));
                        gates.emplace_back(g);
                        for (int i = 0; i < clause.size(); i++) {
                            const CMSat::Lit &lit = clause[i];
                            node_t v1 = literal_node[lit];
                            if(terminals.count(v1) == 0) {
                                terminals.insert(std::make_pair(v1, std::vector<std::pair<std::shared_ptr<gate_t>, uint32_t >>()));
                            }
                            terminals[v1].emplace_back(std::make_pair(g, i));
                            if(terminals_reverse.count(std::make_pair(g, i)) == 0) {
                                terminals_reverse.insert(std::make_pair(std::make_pair(g, i), v1));
                            }
                        }
                    } else if(clause.size() == 3) {
                        std::shared_ptr<gate_t> g(new gate_t(solg::gate3::Gate(solg::gate3::logic::_or)));
                        gates.emplace_back(g);
                        for (int i = 0; i < clause.size(); i++) {
                            const CMSat::Lit &lit = clause[i];
                            node_t v1 = literal_node[lit];
                            if(terminals.count(v1) == 0) {
                                terminals.insert(std::make_pair(v1, std::vector<std::pair<std::shared_ptr<gate_t>, uint32_t >>()));
                            }
                            terminals[v1].emplace_back(std::make_pair(g, i));
                            if(terminals_reverse.count(std::make_pair(g, i)) == 0) {
                                terminals_reverse.insert(std::make_pair(std::make_pair(g, i), v1));
                            }
                        }
                    }
                }


//                for(const std::vector<CMSat::Lit>& clause : clauses) {
//                    node_t out(new node_s(power::none));
//                    out->pwr = power::vcc;
//                    nodes.emplace_back(out);
//                    for(const CMSat::Lit& lit : clause) {
//                        node_t v1 = literal_node[lit];
//                        connect(v1, out);
//                    }
//                }


            }

//            void connect(const node_t& v1, const node_t& out) {
//                node_t r1_gnd(new node_s(power::VM));
//                node_t x1_gnd(new node_s(power::VM));
//
//                edge_t r1 = r1_gnd + v1;
//                edge_t x1 = x1_gnd + v1;
//                edge_t x4 = out + v1;
//
//                x1->comp = component::memristor;
//                r1->comp = component::resistor;
//                x4->comp = component::memristor;
//
//                edges.emplace_back(r1);
//                edges.emplace_back(x1);
//                edges.emplace_back(x4);
//            }

            node_t other(const node_t& n) {
                edge_t e = n->edge;
                if(e->oper1 == n)
                    return e->oper2;
                if(e->oper2 == n)
                    return e->oper1;
                std::abort();
            }

            struct gate_visitor {
                uint32_t t;
                gate_visitor(const int t) : t(t) {};
                solg_t operator()(const solg::gate2::Gate& g) {
                    return g.I[t];
                }

                solg_t operator()(const solg::gate3::Gate& g) {
                    return g.I[t];
                }
            };


            solg_t current(const std::shared_ptr<gate_t>& gate, int i) {
                auto p = std::make_pair(gate, i);
                if(terminals_reverse.count(p) == 0) return 0;
                node_t n = terminals_reverse[p];
                solg_t ret = 0;
                for(const auto& r : terminals[n]) {
                    if(r.first == gate) continue;
                    ret += std::visit(gate_visitor(i), *r.first);
//                    ret += r.first->I[i];
                }

                if(n->edge.get() != nullptr) {
                    const auto &o = other(n);
                    for (const auto &r : terminals[o]) {
                        if (r.first == gate) continue;
                        ret -= std::visit(gate_visitor(i), *r.first);
                    }

                    if(n->edge->oper1 == o) {
                        ret *= -1;
                    }
                }

                return ret;
            }


            struct rk4_gate_visitor {
                solg_t i1,i2, i3;

                rk4_gate_visitor(const solg_t& i1, const solg_t& i2) : i1(i1), i2(i2), i3(0) {}
                rk4_gate_visitor(const solg_t& i1, const solg_t& i2, const solg_t& i3) : i1(i1), i2(i2), i3(i3) {}
                void operator()(solg::gate2::Gate& g) {
                    g.rk4(i1, i2);
                }

                void operator()(solg::gate3::Gate& g) {
                    g.rk4(i1, i2, i3);
                }
            };

            struct current_gate_visitor {

                void operator()(solg::gate2::Gate& g) {
                    g.current_update();
                }

                void operator()(solg::gate3::Gate& g) {
                    g.current_update();
                }
            };

            void solve() {
                for(int i = 0; i < 10; i++) {
#pragma omp parallel
#pragma omp for
                    for(int j = 0; j < gates.size(); j++) {
                        const std::shared_ptr<gate_t> &g = gates[j];
                        solg_t i1 = current(g, 0);
                        solg_t i2 = current(g, 1);
                        solg_t i3 = current(g, 2);

                        std::visit(rk4_gate_visitor(i1, i2, i3), *g);
//                        g->rk4(i1, i2);
                    }
                    std::cout << "updating current" << std::endl;
#pragma omp parallel
#pragma omp for
                    for(int j = 0; j < gates.size(); j++) {
                        const std::shared_ptr<gate_t> &g = gates[j];
                        std::visit(current_gate_visitor(), *g);
//                        g->current_update();
                    }
                }
            }

        };


    }



}


#endif //PCB_SAT_CIRCUIT_H

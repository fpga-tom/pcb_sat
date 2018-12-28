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
#include <Eigen/Sparse>
#include <Eigen/Dense>
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

        typedef Eigen::SparseMatrix<solg_t, Eigen::ColMajor > SpMat;
        typedef Eigen::VectorXd Vec;
        typedef Eigen::Triplet<solg_t > T;

        struct Circuit {

            std::vector<std::unique_ptr<gate_t>> gates;
            std::map<std::pair<uint64_t , uint64_t >, uint64_t > terminal_id;

            SpMat I_mat;
            Vec I_vec;



            Circuit(const std::vector<uint64_t >& linear, const std::vector<std::vector<CMSat::Lit>>& clauses) {
//                std::vector<node_t> nodes;
//                std::vector<edge_t> edges;

//                std::map<CMSat::Lit, node_t> literal_node;
//                std::map<node_t, std::vector<std::pair<std::shared_ptr<gate_t>, uint64_t >>> terminals;
//                std::map<std::pair<std::shared_ptr<gate_t>, uint64_t>, node_t> terminals_reverse;

                std::map<uint64_t , std::vector<std::pair<uint64_t, bool>> > terminal_map;


                uint64_t terminal_count = 0;
                // create nodes
                for (const std::vector<CMSat::Lit> &clause : clauses) {
                    for (const CMSat::Lit &lit : clause) {
                        terminal_count++;
//                        if (literal_node.count(lit) == 0) {
//                            nodes.emplace_back(new node_s());
//                            literal_node.insert(std::make_pair(lit, nodes.back()));
//                        }
                    }
                }
                I_mat = SpMat(terminal_count, terminal_count);
//                I_mat.reserve(Eigen::VectorXi::Constant(terminal_count, 10));
                I_vec = Vec(terminal_count);
//                for (auto it = literal_node.begin(); it != literal_node.end(); ++it) {
//                    CMSat::Lit a = it->first;
//
//                    if (literal_node.count(~a) > 0 && !a.sign()) {
//                        edge_t edge = it->second + literal_node[~a];
//                        edge->comp = component::inverter;
//                        it->second->edge = edge;
//                        literal_node[~a]->edge = edge;
//                        edges.emplace_back(edge);
//                    }
//                }

                terminal_count = 0;
                for (const std::vector<CMSat::Lit> &clause : clauses) {
                    if (clause.size() == 2) {
                        std::unique_ptr<gate_t> g(new gate_t(solg::gate2::Gate(solg::gate2::logic::_or)));
                        for (uint64_t i = 0; i < clause.size(); i++) {
                            const CMSat::Lit &lit = clause[i];
//                            node_t v1 = literal_node[lit];
//                            if (terminals.count(v1) == 0) {
//                                terminals.insert(std::make_pair(v1,
//                                                                std::vector<std::pair<std::shared_ptr<gate_t>, uint64_t >>()));
//                            }
//                            terminals[v1].emplace_back(std::make_pair(g, i));
//                            if (terminals_reverse.count(std::make_pair(g, i)) == 0) {
//                                terminals_reverse.insert(std::make_pair(std::make_pair(g, i), v1));
//                            }
                            std::pair<uint64_t ,uint64_t > p((uint64_t )g.get(), i);
                            terminal_id.insert(std::make_pair(p, terminal_count++));
                            if(terminal_map.count(lit.var()) == 0) {
                                terminal_map.insert(std::make_pair(lit.var(), std::vector<std::pair<uint64_t, bool>>()));
                            }
                            terminal_map[lit.var()].emplace_back(std::make_pair(terminal_id[p], lit.sign()));

                        }
                        gates.emplace_back(std::move(g));
                    } else if (clause.size() == 3) {
                        std::unique_ptr<gate_t> g(new gate_t(solg::gate3::Gate(solg::gate3::logic::_or)));
                        for (uint64_t i = 0; i < clause.size(); i++) {
                            const CMSat::Lit &lit = clause[i];
//                            node_t v1 = literal_node[lit];
//                            if (terminals.count(v1) == 0) {
//                                terminals.insert(std::make_pair(v1,
//                                                                std::vector<std::pair<std::shared_ptr<gate_t>, uint64_t >>()));
//                            }
//                            terminals[v1].emplace_back(std::make_pair(g, i));
//                            if (terminals_reverse.count(std::make_pair(g, i)) == 0) {
//                                terminals_reverse.insert(std::make_pair(std::make_pair(g, i), v1));
//                              }
                            std::pair<uint64_t ,uint64_t > p((uint64_t )g.get(), i);
                            terminal_id.insert(std::make_pair(p, terminal_count++));
                            if(terminal_map.count(lit.var()) == 0) {
                                terminal_map.insert(std::make_pair(lit.var(), std::vector<std::pair<uint64_t, bool>>()));
                            }
                            terminal_map[lit.var()].emplace_back(std::make_pair(terminal_id[p], lit.sign()));
                        }

                        gates.emplace_back(std::move(g));
                    }
                }

                std::vector<T> triplets;
                for(auto it = terminal_map.begin(); it != terminal_map.end(); ++it) {
                    for(const std::pair<uint64_t, bool>& t : it->second) {
                        triplets.emplace_back(T(it->first, t.first, t.second ? -1 : 1));
                    }
                }

                /*
                for (int j = 0; j < gates.size(); j++) {
                    const std::shared_ptr<gate_t> &g = gates[j];
                    const std::vector<CMSat::Lit> &clause = clauses[j];
                    for (uint64_t i = 0; i < clause.size(); i++) {
                        const CMSat::Lit &lit = clause[i];
                        uint64_t tid = terminal_id[std::make_pair((uint64_t )g.get(), i)];
//                        std::cout << terminal_map[lit.var()].size() << std::endl;
                        for(const uint64_t t : terminal_map[lit.var()]) {
                            if(tid == t) continue;
//                            I_mat.insert(t, tid) = lit.sign() ? -1 : 1;
                            triplets.emplace_back(T(t, tid, lit.sign() ? -1 : 1));
                        }
                    }
                }
                 */


                I_mat.setFromTriplets(triplets.begin(), triplets.end());
                I_mat.makeCompressed();
                std::cout << "Creating triplets done " << I_mat.nonZeros() << "/" << triplets.size() << "/" << terminal_count << "/" << (terminal_count*terminal_count) << std::endl;


//                for (int j = 0; j < gates.size(); j++) {
//                    const std::shared_ptr<gate_t> &g = gates[j];
//                    current(triplets, g, 0, terminals, terminals_reverse, terminal_id);
//                    current(triplets, g, 1, terminals, terminals_reverse, terminal_id);
//                    current(triplets, g, 2, terminals, terminals_reverse, terminal_id);
//                }

            }

            node_t other(const node_t& n) {
                edge_t e = n->edge;
                if(e->oper1 == n)
                    return e->oper2;
                if(e->oper2 == n)
                    return e->oper1;
                std::abort();
            }

            struct gate_visitor {
                uint64_t t;
                gate_visitor(const uint64_t t) : t(t) {};
                solg_t operator()(const solg::gate2::Gate& g) {
                    return g.I[t];
                }

                solg_t operator()(const solg::gate3::Gate& g) {
                    return g.I[t];
                }
            };


            void current(std::vector<T>& triplets, const std::shared_ptr<gate_t>& gate, int i,
                           std::map<node_t, std::vector<std::pair<std::shared_ptr<gate_t>, uint64_t >>>& terminals,
            std::map<std::pair<std::shared_ptr<gate_t>, uint64_t >, node_t>& terminals_reverse,
            std::map<std::pair<std::shared_ptr<gate_t>, uint64_t >, uint64_t >& terminal_id) {
                auto p = std::make_pair(gate, i);
                auto tid = terminal_id[p];
                if(terminals_reverse.count(p) == 0) return;
                node_t n = terminals_reverse[p];
                for(const auto& r : terminals[n]) {
                    if(r.first == gate) continue;
                    auto p1 = std::make_pair(r.first, r.second);
                    triplets.emplace_back(T(tid, terminal_id[p1],(n->edge->oper1 == n ? 1 : -1)));
                }

                if(n->edge.get() != nullptr) {
                    const auto &o = other(n);
                    for (const auto &r : terminals[o]) {
                        if (r.first == gate) continue;
                        auto p1 = std::make_pair(r.first, r.second);
                        triplets.emplace_back(T(tid, terminal_id[p1],(n->edge->oper1 == n ? 1 : -1)));
                    }

                }

            }


            struct rk4_gate_visitor {
                solg_t i1,i2, i3;

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
                for(int i = 0; i < 100; i++) {
                    const Eigen::VectorXd I_vec_new = I_mat * I_vec;
#pragma omp parallel
#pragma omp for
                    for(int j = 0; j < gates.size(); j++) {
                        const std::unique_ptr<gate_t> &g = gates[j];

                        solg_t i1 = 0;
                        solg_t i2 = 0;
                        solg_t i3 = 0;

                        auto g1 = std::make_pair((uint64_t )g.get(), 0);
                        if(terminal_id.count(g1) > 0) {
                            i1 = I_vec_new[terminal_id[g1]] - I_vec[terminal_id[g1]];
                        }
                        auto g2 = std::make_pair((uint64_t )g.get(), 1);
                        if(terminal_id.count(g2) > 0) {
                            i2 = I_vec_new[terminal_id[g2]] - I_vec[terminal_id[g2]];
                        }
                        auto g3 = std::make_pair((uint64_t )g.get(), 2);
                        if(terminal_id.count(g3) > 0) {
                            i3 = I_vec_new[terminal_id[g3]] - I_vec[terminal_id[g3]];
                        }

                        std::visit(rk4_gate_visitor(i1, i2, i3), *g);
                    }
                    if(i % 10 == 0)
                        std::cout << "updating current" << std::endl;
#pragma omp parallel
#pragma omp for
                    for(int j = 0; j < gates.size(); j++) {
                        const std::unique_ptr<gate_t> &g = gates[j];
                        std::visit(current_gate_visitor(), *g);
                    }

#pragma omp parallel
#pragma omp for
                    for(int j = 0; j < gates.size(); j++) {
                        const std::unique_ptr<gate_t> &g = gates[j];

                        auto g1 = std::make_pair((uint64_t )g.get(), 0);
                        if(terminal_id.count(g1) > 0) {
                            I_vec[terminal_id[g1]] = std::visit(gate_visitor(0), *g);
                        }
                        auto g2 = std::make_pair((uint64_t )g.get(), 1);
                        if(terminal_id.count(g2) > 0) {
                            I_vec[terminal_id[g2]] = std::visit(gate_visitor(1), *g);
                        }
                        auto g3 = std::make_pair((uint64_t )g.get(), 2);
                        if(terminal_id.count(g3) > 0) {
                            I_vec[terminal_id[g3]] = std::visit(gate_visitor(2), *g);
                        }
                    }
                }

            }

        };


    }



}


#endif //PCB_SAT_CIRCUIT_H

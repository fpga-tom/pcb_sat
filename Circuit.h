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
#include <xxhash.h>
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

        typedef Eigen::SparseMatrix<solg_t, Eigen::RowMajor > SpMat;
        typedef Eigen::VectorXf Vec;
        typedef Eigen::Triplet<solg_t > T;

        struct Circuit {

            std::vector<std::unique_ptr<gate_t>> gates;
            std::map<std::pair<uint64_t , uint32_t >, uint32_t > terminal_id;

            const std::vector<std::vector<CMSat::Lit>>& clauses;

            uint64_t linear_size;

            SpMat I_mat;
            Vec I_vec;
//            std::vector<solg_t> I_vec;



            Circuit(const std::vector<uint64_t >& linear, const std::vector<std::vector<CMSat::Lit>>& clauses) : clauses(clauses), linear_size(linear.size()) {

                std::map<uint32_t , uint32_t > rt_len;
                std::map<uint64_t , std::pair<uint32_t, bool> > terminal_map;
                std::vector<uint32_t > rt_len_keys;

                std::cout << "sizeof(Gate3) " << sizeof(solg::gate3::Gate) << std::endl;
                std::cout << "sizeof(Gate2) " << sizeof(solg::gate2::Gate) << std::endl;
                uint64_t terminal_count = 0;
                // create nodes
                for (const std::vector<CMSat::Lit> &clause : clauses) {
                    for (const CMSat::Lit &lit : clause) {
                        terminal_count++;

                    }
                }
                I_mat = SpMat(terminal_count, terminal_count);
                I_vec = Vec(terminal_count);
//                I_vec = std::vector<solg_t >(terminal_count);

                terminal_count = 0;
                float run = 0;
                uint32_t count = 0;
                for (const std::vector<CMSat::Lit> &clause : clauses) {
                    count++;
                    if(run++ > 1000) {
                        std::cout << "\r" << 100*(count / (double)clauses.size());
                        run = 0;
                    }
                    if (clause.size() == 2) {
                        std::unique_ptr<gate_t> g(new gate_t(solg::gate2::Gate(solg::gate2::logic::_or)));
                        for (uint64_t i = 0; i < clause.size(); i++) {
                            const CMSat::Lit &lit = clause[i];

                            if(rt_len.count(lit.var()) == 0) {
                                rt_len[lit.var()] = 0;
                                rt_len_keys.emplace_back(lit.var());
                            }
                            rt_len[lit.var()]++;

                            uint64_t hash = lit.var();
                            for(int j = 0; j < rt_len[lit.var()]; j++) {
                                hash = XXH64(&hash, 8, 0);
                            }

                            std::pair<uint64_t ,uint64_t > p((uint64_t )g.get(), i);
                            terminal_id.insert(std::make_pair(p, terminal_count++));
                            if(terminal_map.count(hash) > 0) {
                                std::cout << "collision" << std::endl;
                            }
                            terminal_map.insert(std::make_pair(hash, std::make_pair(terminal_id[p], lit.sign())));

                        }
                        gates.emplace_back(std::move(g));
                    } else if (clause.size() == 3) {
                        std::unique_ptr<gate_t> g(new gate_t(solg::gate3::Gate(solg::gate3::logic::_or)));
                        for (uint64_t i = 0; i < clause.size(); i++) {
                            const CMSat::Lit &lit = clause[i];
                            if(rt_len.count(lit.var()) == 0) {
                                rt_len[lit.var()] = 0;
                                rt_len_keys.emplace_back(lit.var());
                            }
                            rt_len[lit.var()]++;

                            uint64_t hash = lit.var();
                            for(int j = 0; j < rt_len[lit.var()]; j++) {
                                hash = XXH64(&hash, 8, 0);
                            }

                            std::pair<uint64_t ,uint64_t > p((uint64_t )g.get(), i);
                            terminal_id.insert(std::make_pair(p, terminal_count++));
                            if(terminal_map.count(hash) > 0) {
                                std::cout << "collision" << std::endl;
                            }
                            terminal_map.insert(std::make_pair(hash, std::make_pair(terminal_id[p], lit.sign())));
                        }

                        gates.emplace_back(std::move(g));
                    }
                }

                std::vector<T> triplets;
                for(int k = 0; k < rt_len_keys.size(); k++) {
                    auto first = rt_len_keys[k];
                    uint64_t hash = first;
                    for(int j = 0; j < rt_len[first]; j++) {
                        hash = XXH64(&hash, 8, 0);
                        auto& p = terminal_map[hash];
                        triplets.emplace_back(T(first, p.first, p.second ? -1 : 1));
                    }
                }
                /*
                for(auto it = terminal_map.begin(); it != terminal_map.end(); ++it) {
                    for(const std::pair<uint64_t, bool>& t : it->second) {
                        triplets.emplace_back(T(it->first, t.first, t.second ? -1 : 1));
                    }
                }
                 */


                I_mat.setFromTriplets(triplets.begin(), triplets.end());
                I_mat.makeCompressed();
                std::cout << "Creating triplets done " << terminal_id.size() << "/" << terminal_map.size() << "/" << I_mat.nonZeros() << "/" << triplets.size() << "/" << terminal_count << "/" << (terminal_count*terminal_count) << std::endl;

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
//                Eigen::VectorXd I_vec_new(linear_size);
//                std::vector<solg_t > I_vec_new(linear_size);
                for(int i = 0; i < 100; i++) {
                    Vec I_vec_new = I_mat * I_vec;

//#pragma omp parallel
//#pragma omp for
//                    for(int k = 0; k < rt_len_keys.size(); k++) {
//                        auto first = rt_len_keys[k];
//                        uint64_t hash = first;
//                        for(int j = 0; j < rt_len[first]; j++) {
//                            hash = XXH64(&hash, 8, 0);
//                            auto& p = terminal_map[hash];
//                            if(j == 0)
//                                I_vec_new[first] = I_vec[p.first] * (p.second ? -1 : 1);
//                            else
//                                I_vec_new[first] += I_vec[p.first] * (p.second ? -1 : 1);
//                        }
//                    }

#pragma omp parallel
#pragma omp for
                    for(int j = 0; j < gates.size(); j++) {
                        const std::unique_ptr<gate_t> &g = gates[j];
                        const std::vector<CMSat::Lit>& clause = clauses[j];

                        solg_t i1 = 0;
                        solg_t i2 = 0;
                        solg_t i3 = 0;

                        auto g1 = std::make_pair((uint64_t )g.get(), 0);
                        if(terminal_id.count(g1) > 0) {
                            i1 = I_vec_new[clause[0].var()] - I_vec[terminal_id[g1]];
                        }
                        auto g2 = std::make_pair((uint64_t )g.get(), 1);
                        if(terminal_id.count(g2) > 0) {
                            i2 = I_vec_new[clause[1].var()] - I_vec[terminal_id[g2]];
                        }
                        auto g3 = std::make_pair((uint64_t )g.get(), 2);
                        if(terminal_id.count(g3) > 0) {
                            i3 = I_vec_new[clause[2].var()] - I_vec[terminal_id[g3]];
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

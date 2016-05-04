/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*        This file is part of the program PolySCIP                          */
/*                                                                           */
/*    Copyright (C) 2012-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  PolySCIP is distributed under the terms of the ZIB Academic License.     */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with PolySCIP; see the file LICENCE.                               */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "weight_space_polyhedron.h"

#include <algorithm> // std::transform
#include <functional> //std::bind1st
#include <iostream>
#include <list>
#include <memory> // std::shared_ptr
#include <numeric> // std::inner_product
#include <ostream>
#include <utility> // std::move, std::pair
#include <vector>

#include "lemon/list_graph.h"

#include "polyscip.h"
#include "weight_space_facet.h"
#include "weight_space_vertex.h"

using std::inner_product;
using std::list;
using std::make_shared;
using std::ostream;
using std::transform;
using std::vector;

namespace polyscip {


    using OutcomeType = Polyscip::OutcomeType;
    using ValueType = Polyscip::ValueType;
    using WeightType = Polyscip::WeightType;
    using ResContainer = Polyscip::ResultContainer;
    using FacetContainer = WeightSpaceVertex::FacetContainer;

    WeightSpacePolyhedron::WeightSpacePolyhedron(unsigned num_objs,
                                                 OutcomeType point,
                                                 ValueType weighted_obj_val,
                                                 std::pair<bool, WeightType::size_type> unit_weight_info,
                                                 const ResContainer &initial_rays)
            : curr_investigated_vertex_(nullptr),
              skeleton_(),
              nodes_to_vertices_(skeleton_) {
        // create initial boundary facets w_i > 0 for i \in [num_objs]
        auto facets = FacetContainer{};
        for (auto i = 0; i < num_objs; ++i)
            facets.emplace_back(make_shared<const WeightSpaceFacet>(num_objs, i));
        createInitialVertices(num_objs, std::move(point), weighted_obj_val, std::move(facets));
        createInitialSkeleton();
        for (const auto &pair : initial_rays) { //incorporate found unbounded rays
            auto wspoly_changed = updateInitialWeightSpacePolyhedron(pair.first);
            if (!wspoly_changed)
                print(pair.second,
                      {"INITIAL WEIGHT SPACE POLYHEDRON UNCHANGED: no vertex found with weight = "});
        }
        if (unit_weight_info.first)
            setMarkedVertex(unit_weight_info.second);
    }

    WeightSpacePolyhedron::~WeightSpacePolyhedron() {
        for (auto v_ptr : marked_vertices_)
            delete v_ptr;
        for (auto v_ptr : unmarked_vertices_)
            delete v_ptr;
        for (auto v_ptr : obsolete_vertices_)
            delete v_ptr;
        delete curr_investigated_vertex_;
    }

    void WeightSpacePolyhedron::test() {
        auto iter = unmarked_vertices_.front();
        std::cout << "first vertex:\n";
        iter->print(std::cout, true);
        auto iter2 = unmarked_vertices_.back();
        std::cout << "second vertex:\n";
        iter2->print(std::cout, true);
        auto w = WeightSpaceVertex(iter, iter2, {-2,2,2}, true);
        std::cout << "new vertex\n";
        //todo check definiton of facets!
        w.print(std::cout, true);
    }

    bool WeightSpacePolyhedron::hasUntestedWeight() const {
        return !unmarked_vertices_.empty();
    }

    WeightType WeightSpacePolyhedron::getUntestedWeight() {
        assert(!unmarked_vertices_.empty());
        curr_investigated_vertex_ = unmarked_vertices_.front();
        unmarked_vertices_.pop_front();
        return curr_investigated_vertex_->getWeight();
    }

    bool WeightSpacePolyhedron::updateInitialWeightSpacePolyhedron(const OutcomeType &ray) {
        auto obsolete_vertices = vector<WeightSpaceVertex*>{};
        auto non_obsolete_vertices = vector<WeightSpaceVertex*>{};
        // iterate over all vertices
        assert (marked_vertices_.empty());
        for (auto vertex : unmarked_vertices_) {
            if (vertex->isMadeObsolete(ray,0.))
                obsolete_vertices.push_back(vertex);
            else
                non_obsolete_vertices.push_back(vertex);
        }
        assert (!obsolete_vertices.empty() && obsolete_vertices.size() < unmarked_vertices_.size());
        auto new_vertices = vector<WeightSpaceVertex*>{};
        for (auto obs : obsolete_vertices) {
            for (auto non_obs : non_obsolete_vertices) {
                // create new vertices between obsolete and non-obsolete vertices
                //todo implementation
            }
        }




//        for (auto vertex : unmarked_vertices_) {
//            if (vertex->hasSameWeight(ray_weight)) { // find vertex with same weight as ray
//                auto node = vertices_to_nodes_.at(vertex);
//                // iterate over all adjacent nodes
//                for (Graph::IncEdgeIt edge(skeleton_, node); edge != lemon::INVALID; ++edge) {
//                    auto adj_node = skeleton_.oppositeNode(node, edge);
//                    auto adj_vertex = nodes_to_vertices_[adj_node];
//
//
//                    adj_vertex->print(std::cout,true);
//                    auto new_weight = calculateWeight(std::move(adj_vertex->getWeight()),
//                                                      std::move(vertex->getWeight()),
//                                                      ray);
//                    print(new_weight, {"new computed weight ="}, std::cout);
//                    //TODO
//                }
//
//            }
//        }
        return false;

    }



    void WeightSpacePolyhedron::createInitialVertices(unsigned num_objs,
                                                      OutcomeType point,
                                                      ValueType weighted_obj_val,
                                                      FacetContainer boundary_facets) {
        // make facet: point \cdot w >= weighted_obj_val to facets_
        auto point_facet = make_shared<const WeightSpaceFacet>(point,1.0);
        // create initial weight space vertices
        for (auto i = 0; i < num_objs; ++i) {
            auto facets = FacetContainer(begin(boundary_facets), end(boundary_facets));
            facets.at(i) = point_facet; // replace i-th boundary facet
            auto weight = WeightType(num_objs, 0.);
            weight[i] = 1.; // set unit weight
            unmarked_vertices_.push_back(new WeightSpaceVertex(std::move(facets),
                                                               std::move(weight),
                                                               point.at(i)));
        }
    }

    void WeightSpacePolyhedron::createInitialSkeleton() {
        for (auto &vertex : unmarked_vertices_) {
            Node new_node = skeleton_.addNode();
            nodes_to_vertices_[new_node] = vertex; // map from node to vertex
            for (const auto &pair : vertices_to_nodes_) // add edge between all previous nodes and new_node
                skeleton_.addEdge(pair.second, new_node);
            // above for loop needs to come before insertion in order to keep ordering:
            // second node has edge to first node; third node has edge to first node and second node; etc
            vertices_to_nodes_.insert({vertex, new_node}); // map from vertex to node
        }
    }

    void WeightSpacePolyhedron::setMarkedVertex(WeightType::size_type unit_weight_index) {
        for (auto vertex_it = begin(unmarked_vertices_); vertex_it != end(unmarked_vertices_); ++vertex_it) {
            if ((*vertex_it)->hasUnitWeight(unit_weight_index)) {
                marked_vertices_.push_back(*vertex_it);
                unmarked_vertices_.erase(vertex_it);
                return;
            }
        }
    }


    void WeightSpacePolyhedron::incorporateOutcome(const OutcomeType &outcome,
                                                   const WeightType &weight,
                                                   bool outcome_is_ray) {
        assert(curr_investigated_vertex_ != nullptr);
        assert(curr_investigated_vertex_->hasSameWeight(weight));
        if (isNewOutcome(outcome, outcome_is_ray)) {
            incorporateNewOutcome(outcome, outcome_is_ray);
        }
        else {
            assert(!outcome_is_ray);
            incorporateOldOutcome(outcome);
        }
    }


    bool WeightSpacePolyhedron::isNewOutcome(const Polyscip::OutcomeType &outcome,
                                             bool outcome_is_ray) {
        if (outcome_is_ray)
            return curr_investigated_vertex_->isMadeObsolete(outcome, 0.);
        else
            return curr_investigated_vertex_->isMadeObsolete(outcome);
    }

    void WeightSpacePolyhedron::incorporateOldOutcome(const Polyscip::OutcomeType &outcome) {
        marked_vertices_.push_back(curr_investigated_vertex_);
        curr_investigated_vertex_ = nullptr;
    }

    void WeightSpacePolyhedron::incorporateNewOutcome(const OutcomeType &outcome,
                                                      bool outcome_is_ray) {
        // adjacent nodes of obsolete node might also be made obsolete by outcome
        auto unscanned_nodes = list < Node > {vertices_to_nodes_.at(curr_investigated_vertex_)};
        curr_investigated_vertex_ = nullptr;
        auto cut_edges = list<Graph::Edge>(); // container for edges between obsolete and non-obsolete nodes
        while (!unscanned_nodes.empty()) { // find all vertices which are made obsolete
            auto node = unscanned_nodes.front();
            for (Graph::IncEdgeIt edge(skeleton_, node); edge != lemon::INVALID; ++edge) {
                auto adj_node = skeleton_.oppositeNode(node, edge);
                auto adj_vertex = nodes_to_vertices_[adj_node];
                // check whether adjacent vertex is made obsolete by outcome
                auto is_obsolete = (outcome_is_ray) ? adj_vertex->isMadeObsolete(outcome, 0.) :
                                   adj_vertex->isMadeObsolete(outcome);
                if (is_obsolete) {
                    //todo
                    ;
                }
                else {
                    //todo cut_edges
                }
            }
            unscanned_nodes.pop_front(); // discard node
        }
    }

    void WeightSpacePolyhedron::printUnmarkedVertices(ostream &os, bool printFacets) const {
        printVertices(unmarked_vertices_, {"UNMARKED VERTICES:"}, os, printFacets);
    }

    void WeightSpacePolyhedron::printMarkedVertices(ostream &os, bool printFacets) const {
        printVertices(marked_vertices_, {"MARKED VERTICES:"}, os, printFacets);
    }

    void WeightSpacePolyhedron::printObsoleteVertices(ostream &os, bool printFacets) const {
        printVertices(obsolete_vertices_, {"OBSOLETE VERTICES:"}, os, printFacets);
    }

}


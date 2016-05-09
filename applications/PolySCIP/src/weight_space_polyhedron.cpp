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
#include <iostream>
#include <iterator> // std::next
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
using std::pair;
using std::transform;
using std::vector;

namespace polyscip {


    using OutcomeType = Polyscip::OutcomeType;
    using ValueType = Polyscip::ValueType;
    using WeightType = Polyscip::WeightType;
    using RayContainer = Polyscip::RayContainer;
    using FacetContainer = WeightSpaceVertex::FacetContainer;

    WeightSpacePolyhedron::WeightSpacePolyhedron(unsigned num_objs,
                                                 OutcomeType point,
                                                 ValueType weighted_obj_val,
                                                 const RayContainer& initial_rays,
                                                 std::pair<bool, WeightType::size_type> unit_weight_info)
            : curr_investigated_vertex_(nullptr),
              skeleton_(),
              nodes_to_vertices_(skeleton_) {
        // create initial boundary facets w_i > 0 for i \in [num_objs]
        auto facets = FacetContainer{};
        for (auto i = 0; i < num_objs; ++i)
            facets.emplace_back(make_shared<const WeightSpaceFacet>(num_objs, i));
        createInitialVertices(num_objs, std::move(point), weighted_obj_val, std::move(facets));
        createInitialSkeleton();
        // incorporate non-dominated rays computed in initial phase
        for (const auto& ray : initial_rays) {
            auto wspoly_changed = updateInitialWeightSpacePolyhedron(ray);
            if (!wspoly_changed)
                print(ray,
                      {"INITIAL WEIGHT SPACE POLYHEDRON UNCHANGED: no vertex made obsolete by ray: "});
        }
        if (unit_weight_info.first) // mark initial vertex correspoding to non-dominated point
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

    bool WeightSpacePolyhedron::hasUntestedWeight() const {
        return !unmarked_vertices_.empty();
    }

    WeightType WeightSpacePolyhedron::getUntestedWeight() {
        assert(!unmarked_vertices_.empty());
        curr_investigated_vertex_ = unmarked_vertices_.front();
        unmarked_vertices_.pop_front();
        return curr_investigated_vertex_->getWeight();
    }

    void WeightSpacePolyhedron::addNewVerticesToSkeleton(const vector<pair<WeightSpaceVertex*, WeightSpaceVertex*>>& vertex_pairs) {
        for (const auto& pair : vertex_pairs) {
            auto new_node = skeleton_.addNode();
            nodes_to_vertices_[new_node] = pair.first;
            vertices_to_nodes_.insert({pair.first, new_node});
            skeleton_.addEdge(pair.second, new_node);
        }
        // add edges of complete graph of new
        for (auto it=begin(vertex_pairs); it!=end(vertex_pairs); ++it) {
            for (auto it2=std::next(it); it2!=end(vertex_pairs); ++it2) {
                skeleton_.addEdge(getNode(it->first), getNode(it2->first));
            }
        }


    }

    void WeightSpacePolyhedron::deleteVertexFromSkeleton(WeightSpaceVertex *v) {
        skeleton_.erase(getNode(v));
        auto ret = vertices_to_nodes_.erase(v);
        assert (ret == 1);
    }

    void WeightSpacePolyhedron::incorporateObsoleteVertex(WeightSpaceVertex* v) {
        auto size_before = unmarked_vertices_.size();
        unmarked_vertices_.remove(v);
        assert( --size_before == unmarked_vertices_.size());
        obsolete_vertices_.push_back(v);
    }

    bool WeightSpacePolyhedron::updateInitialWeightSpacePolyhedron(const OutcomeType &ray) {
        auto obsolete = vector<WeightSpaceVertex*>{};
        // iterate over all non-obsolete vertices
        assert (marked_vertices_.empty());
        for (auto vertex : unmarked_vertices_) {
            if (vertex->isMadeObsolete(ray,0.))
                obsolete.push_back(vertex);
        }
        if (obsolete.empty()) {
            return false;
        }
        else {
            assert (obsolete.size() < unmarked_vertices_.size());
            auto new_vertices = vector<pair<WeightSpaceVertex*,WeightSpaceVertex*>>{}; // pair: new vertex, adjacent vertex
            // create new vertices
            for (auto obs : obsolete) {
                for (Graph::IncEdgeIt edge(skeleton_, getNode(obs)); edge != lemon::INVALID; ++edge) {
                    auto adj_node = skeleton_.oppositeNode(getNode(obs), edge);
                    auto adj_vertex = getVertex(adj_node);
                    if (!adj_vertex->isMadeObsolete(ray, 0.)) {
                        auto new_vertex = new WeightSpaceVertex(obs, adj_vertex, ray,
                                                                true); // new vertex between obsolete and non-obsolete verts
                        new_vertices.emplace_back({new_vertex, adj_vertex});
                    }
                }
            }
            // incorporate new vertices into 1-skeleton
            for (const auto& pair : new_vertices) {

            }
            // delete obsolete vertices and corresponding nodes

            return true;
        }
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
            auto new_node = skeleton_.addNode();
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


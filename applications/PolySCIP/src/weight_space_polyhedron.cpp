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
#include <cstddef>
#include <iostream>
#include <iterator> // std::next
#include <forward_list>
#include <memory> // std::make_shared
#include <numeric> // std::inner_product
#include <ostream>
#include <utility> // std::move, std::pair
#include <vector>

#include "global_functions.h"
#include "lemon/list_graph.h"
#include "objscip/objscip.h"
#include "polyscip_types.h"
#include "weight_space_facet.h"
#include "weight_space_vertex.h"

using std::inner_product;
using std::make_shared;
using std::ostream;
using std::pair;
using std::shared_ptr;
using std::transform;
using std::vector;

namespace polyscip {

    WeightSpacePolyhedron::WeightSpacePolyhedron(std::size_t num_objs,
                                                 const OutcomeType& point)
            : curr_investigated_vertex_(nullptr),
              skeleton_(),
              nodes_to_vertices_(skeleton_) {
        // create initial boundary facets w_i > 0 for i \in [num_objs]
        auto facets = FacetContainer{};
        for (decltype(num_objs) i=0; i < num_objs; ++i)
            facets.push_back(make_shared<const WeightSpaceFacet>(num_objs, i));
        createInitialVertices(num_objs, point, std::move(facets));
        createInitialSkeleton();
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
        assert (!unmarked_vertices_.empty());
        assert (curr_investigated_vertex_== nullptr);
        curr_investigated_vertex_ = unmarked_vertices_.front();
        unmarked_vertices_.pop_front();
        return curr_investigated_vertex_->getWeight();
    }

    ValueType WeightSpacePolyhedron::getUntestedVertexWOV(const WeightType& untested_weight) const {
        assert (curr_investigated_vertex_ != nullptr);
        assert (curr_investigated_vertex_->hasSameWeight(untested_weight));
        return curr_investigated_vertex_->getCurrentWOV();
    }

    void WeightSpacePolyhedron::addToSkeleton(const vector<WeightSpaceVertex*>& new_vertices,
                                              const vector<pair<WeightSpaceVertex*, Node> >& new_edges) {
        // add new nodes
        for (const auto& v : new_vertices) {
            auto new_node = skeleton_.addNode();
            nodes_to_vertices_[new_node] = v;
            vertices_to_nodes_.insert({v, new_node});
        }
        // add new edges
        for (const auto& edge : new_edges)
            skeleton_.addEdge(getNode(edge.first), edge.second);
    }


    void WeightSpacePolyhedron::deleteFromSkeleton(WeightSpaceVertex *v) {
        skeleton_.erase(getNode(v));
        auto ret = vertices_to_nodes_.erase(v);
        assert (ret == 1);
    }

    void WeightSpacePolyhedron::setStatusToObsolete(WeightSpaceVertex* v) {
        if (v != curr_investigated_vertex_) {
            auto size_before = unmarked_vertices_.size();
            unmarked_vertices_.remove(v);
            assert( --size_before == unmarked_vertices_.size());
        }
        obsolete_vertices_.push_back(v);
    }

    bool WeightSpacePolyhedron::updateInitialWSP(SCIP* scip, std::size_t unit_weight_index, const OutcomeType& outcome, bool outcome_is_ray) {
        auto vertex = marked_vertices_.at(unit_weight_index);
        assert (vertex->hasUnitWeight(unit_weight_index));
        if (isVertexObsolete(scip, vertex, outcome, outcome_is_ray)) {
            auto new_vertices = vector<WeightSpaceVertex*> {};
            auto new_edges = vector< pair<WeightSpaceVertex*, Node> > {}; // pair: new vertex, node of adjacent vertex
            auto node = getNode(vertex);
            for (std::size_t i=0; i<unit_weight_index; ++i) {
                auto adj_vertex = marked_vertices_.at(i);
                assert (adj_vertex->hasUnitWeight(i));
                if (!isVertexObsolete(scip, adj_vertex, outcome, outcome_is_ray)) {
                    auto adj_node = getNode(adj_vertex);
                    auto edge = lemon::findEdge(skeleton_, node, adj_node);
                    skeleton_.erase(edge); // delete edge between node and adj_node
                    auto new_vertex = new WeightSpaceVertex(scip, vertex, adj_vertex, outcome,
                                                            outcome_is_ray); // new vertex between obsolete and non-obsolete verts
                    unmarked_vertices_.push_back(new_vertex); // add to unmarked vertices
                    new_vertices.push_back(new_vertex);
                    new_edges.push_back({new_vertex, node});
                    new_edges.push_back({new_vertex, adj_node});
                }
            }
            addToSkeleton(new_vertices, new_edges);
            return true;
        }
        else {
            return false;
        }

    }


    void WeightSpacePolyhedron::updateWeightSpacePolyhedron(SCIP* scip,
                                                            const vector<WeightSpaceVertex*>& obsolete_vertices,
                                                            const OutcomeType& outcome,
                                                            bool outcome_is_ray) {
        auto new_vertices = vector<WeightSpaceVertex*> {};
        auto new_edges = vector< pair<WeightSpaceVertex*, Node> > {}; // pair: new vertex, node of adjacent vertex
        // create new vertices between obsolete and non-obsolete vertices
        for (auto obs_vertex : obsolete_vertices) {
            for (Graph::IncEdgeIt edge(skeleton_, getNode(obs_vertex)); edge != lemon::INVALID; ++edge) {
                auto adj_node = skeleton_.oppositeNode(getNode(obs_vertex), edge);
                auto adj_vertex = getVertex(adj_node);
                auto adj_is_obsolete = isVertexObsolete(scip, adj_vertex, outcome, outcome_is_ray);
                if (!adj_is_obsolete) {
                    auto new_vertex = new WeightSpaceVertex(scip, obs_vertex, adj_vertex, outcome,
                                                            outcome_is_ray); // new vertex between obsolete and non-obsolete verts
                    unmarked_vertices_.push_back(new_vertex); // add to unmarked vertices //todo double entry
                    new_vertices.push_back(new_vertex);
                    new_edges.push_back({new_vertex, adj_node});
                }
            }
        }
        addToSkeleton(new_vertices, new_edges); // incorporate new vertices into 1-skeleton
        addCliqueEdgesToSkeleton(new_vertices);
        for (auto obs_vertex : obsolete_vertices) {
            if (outcome_is_ray || !obs_vertex->isCorner() || obs_vertex == curr_investigated_vertex_) {
                setStatusToObsolete(obs_vertex);
                deleteFromSkeleton(obs_vertex);
            }
        }
    }

    std::size_t WeightSpacePolyhedron::getNumberOfGraphNodes() const {
        std::size_t size = 0;
        for (Graph::NodeIt n(skeleton_); n != lemon::INVALID; ++n)
            ++size;
        return size;
    }

    bool WeightSpacePolyhedron::isVertexObsolete(SCIP* scip,
                                                 const WeightSpaceVertex* vertex,
                                                 const OutcomeType& outcome,
                                                 bool outcome_is_ray) {
        if (outcome_is_ray)
            return SCIPisLT(scip, vertex->getWeightedOutcome(outcome), 0.) == TRUE;
        else
            return SCIPisLT(scip, vertex->getWeightedOutcome(outcome), vertex->getCurrentWOV()) == TRUE;
    }

    void WeightSpacePolyhedron::createInitialVertices(std::size_t num_objs,
                                                      const OutcomeType& point,
                                                      FacetContainer boundary_facets) {
        // make facet: point \cdot w >= weighted_obj_val to facets_
        auto point_facet = make_shared<const WeightSpaceFacet>(point,1.0);
        // create initial weight space vertices
        for (decltype(num_objs) i = 0; i < num_objs; ++i) {
            auto facets = FacetContainer(boundary_facets);
            facets.at(i) = point_facet; // replace i-th boundary facet
            auto weight = WeightType(num_objs, 0.);
            weight[i] = 1.; // set unit weight
            marked_vertices_.push_back(new WeightSpaceVertex(std::move(facets),
                                                             std::move(weight),
                                                             point.at(i)));
        }
    }

    void WeightSpacePolyhedron::createInitialSkeleton() {
        for (auto vertex : unmarked_vertices_) {
            auto new_node = skeleton_.addNode();
            nodes_to_vertices_[new_node] = vertex;
            for (const auto& v : vertices_to_nodes_) // add edge between all previous nodes and new_node
                skeleton_.addEdge(v.second, new_node);
            // above for-loop needs to come before insertion in order to keep ordering:
            // second node has edge to first node; third node has edge to first node and second node; etc
            vertices_to_nodes_.insert({vertex, new_node});
        }
    }

    vector<WeightSpaceVertex*> WeightSpacePolyhedron::computeObsoleteVertices(SCIP* scip,
                                                                              const OutcomeType& outcome,
                                                                              bool outcome_is_ray) {
        auto obsolete = vector<WeightSpaceVertex*> {curr_investigated_vertex_};
        auto unscanned = std::forward_list<WeightSpaceVertex*> {curr_investigated_vertex_};
        while (!unscanned.empty()) {
            auto obs_vertex = unscanned.front();
            unscanned.pop_front();
            std::cout << "currently considered: ";
            obs_vertex->print(std::cout);
            std::cout << "\n";
            for (Graph::IncEdgeIt edge(skeleton_, getNode(obs_vertex)); edge != lemon::INVALID; ++edge) {
                auto adj_node = skeleton_.oppositeNode(getNode(obs_vertex), edge);
                auto adj_vertex = getVertex(adj_node);
                std::cout << "adjacent verts: ";
                adj_vertex->print(std::cout);
                auto adj_is_obsolete = isVertexObsolete(scip, adj_vertex, outcome, outcome_is_ray);
                if (adj_is_obsolete && std::find(begin(obsolete), end(obsolete), adj_vertex) == end(obsolete)) {
                    std::cout << " is obsolete vertex.\n ";
                    obsolete.push_back(adj_vertex);
                    unscanned.push_front(adj_vertex);
                }
            }
        }
        assert (std::search(begin(marked_vertices_),
                            end(marked_vertices_),
                            begin(obsolete),
                            end(obsolete)) == end(marked_vertices_)); // assert no marked vertex is made obsolete by outcome
        return obsolete;
    }

    vector<WeightSpaceVertex*> WeightSpacePolyhedron::computeObsoleteVerticesWithCompleteLoop(SCIP* scip,
                                                                                              const OutcomeType& outcome,
                                                                                              bool outcome_is_ray) {
        for (auto vertex : marked_vertices_) {
            assert (!isVertexObsolete(scip, vertex, outcome, outcome_is_ray)); // assert no marked vertex is made obsolete by outcome
        }
        auto obsolete = vector<WeightSpaceVertex*> {};
        if (curr_investigated_vertex_ != nullptr)
            obsolete.push_back(curr_investigated_vertex_);
        for (auto vertex : unmarked_vertices_) {
            if (isVertexObsolete(scip, vertex, outcome, outcome_is_ray)) {
                obsolete.push_back(vertex);
                std::cout << "obsolete vertex: ";
                vertex->print(std::cout);
            }
        }
        assert (!obsolete.empty());
        return obsolete;
    }

    void WeightSpacePolyhedron::incorporateNewOutcome(SCIP* scip,
                                                      bool completeLoopForObsolete,
                                                      const WeightType& used_weight,
                                                      const OutcomeType &outcome,
                                                      bool outcome_is_ray) {
        assert (curr_investigated_vertex_ != nullptr);
        assert (curr_investigated_vertex_->hasSameWeight(used_weight));
        std::cout << "obs with complete loop: ";
        auto obsolete = completeLoopForObsolete ? computeObsoleteVerticesWithCompleteLoop(scip, outcome, outcome_is_ray) :
                        computeObsoleteVertices(scip, outcome, outcome_is_ray);
        std::cout << "\n obs without complete loop: ";
        auto obsolete2 = computeObsoleteVertices(scip, outcome, outcome_is_ray);
        std::cout << "size1 = " << obsolete.size();
        std::cout << "\nsize2 = " << obsolete2.size();

        //for_each(begin(obsolete), end(obsolete), [&obsolete2](WeightSpaceVertex* v){assert (std::find(begin(obsolete2), end(obsolete2), v) != end(obsolete2));});
        std::cout << "\n";
        assert (obsolete.size() == obsolete2.size());
        updateWeightSpacePolyhedron(scip, obsolete, outcome, outcome_is_ray);
        resetCurrentInvestigatedVertex(); // requirement: call after updateWeightSpacePolyhedron;
    }


    void WeightSpacePolyhedron::incorporateKnownOutcome(const WeightType& used_weight) {
        assert (curr_investigated_vertex_ != nullptr);
        assert (curr_investigated_vertex_->hasSameWeight(used_weight));
        marked_vertices_.push_back(curr_investigated_vertex_);
        resetCurrentInvestigatedVertex();
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


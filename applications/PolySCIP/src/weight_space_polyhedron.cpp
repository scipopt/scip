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
#include <cstddef> // std::size_t
#include <iostream>
#include <forward_list>
#include <functional> // std::negate
#include <memory> // std::make_shared
#include <numeric> // std::inner_product
#include <ostream>
#include <utility> // std::move, std::pair
#include <vector>

#include "double_description_method.h"
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
using std::size_t;
using std::transform;
using std::vector;

namespace polyscip {

    using DD = DoubleDescription;

    /*WeightSpacePolyhedron::WeightSpacePolyhedron(std::size_t num_objs,
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
    }*/

    WeightSpacePolyhedron::WeightSpacePolyhedron(SCIP* scip,
                                                 DD::V_RepContainer v_rep,
                                                 DD::H_RepContainer h_rep,
                                                 DD::V_RepAdjacencyContainer adjacencies)
            : curr_investigated_vertex_(nullptr),
              skeleton_(),
              nodes_to_vertices_(skeleton_) {
        createInitialFacets(h_rep);
        createInitialVerticesAndSkeleton(scip, v_rep, adjacencies);
    }

    void WeightSpacePolyhedron::createInitialFacets(DD::H_RepContainer h_rep) {
        for (auto& h : h_rep) {
            facets_.emplace_back(make_shared<WeightSpaceFacet>(h.first, h.second));
        }
    }

    void WeightSpacePolyhedron::createInitialVerticesAndSkeleton(SCIP* scip,
                                                                 DD::V_RepContainer v_rep,
                                                                 DD::V_RepAdjacencyContainer adjacencies) {
        auto vrep_to_node = vector<Node> {};
        for (auto& v : v_rep) {
            auto inc_facets = computeIncidentFacets(scip, v);
            auto vertex = new WeightSpaceVertex(std::move(inc_facets),
                                                std::move(v.first),
                                                std::move(v.second));
            std::cout << "new vertex = ";
            vertex->print(std::cout, true);
            auto node = skeleton_.addNode();
            vrep_to_node.push_back(node); // maps v to node; used below for adding adjacent edges
            nodes_to_vertices_[node] = vertex;
            vertices_to_nodes_.insert({vertex, node});
            if (vertex->hasUnitWeight())
                marked_vertices_.push_back(vertex);
            else
                unmarked_vertices_.push_back(vertex);

        }
        // add edge between adjacent vertices/nodes
        assert (v_rep.size() == adjacencies.size());
        for (size_t i=0; i<v_rep.size(); ++i) {
            auto node = vrep_to_node[i];
            for (const auto& adj_ind : adjacencies[i]) {
                assert (i != adj_ind);
                skeleton_.addEdge(node, vrep_to_node[adj_ind]);
            }
        }
    }


    /*void WeightSpacePolyhedron::createInitialSkeleton() {
        for (auto vertex : marked_vertices_) {
            auto new_node = skeleton_.addNode();
            nodes_to_vertices_[new_node] = vertex;
            for (const auto& v : vertices_to_nodes_) // add edge between all previous nodes and new_node
                skeleton_.addEdge(v.second, new_node);
            // above for-loop needs to come before insertion in order to keep ordering:
            // second node has edge to first node; third node has edge to first node and second node; etc
            vertices_to_nodes_.insert({vertex, new_node});
        }
    }*/

    WeightSpacePolyhedron::FacetContainer WeightSpacePolyhedron::computeIncidentFacets(SCIP* scip,
                                                                                       const DD::V_RepT& v) const {
        auto facets = FacetContainer {};
        for (const auto& f : facets_) {
            if (SCIPisEQ(scip, f->getWeightedWeight(v.first), f->getWOVCoeff()*v.second))
                facets.push_back(f);
        }
        return facets;
    }



    /*void WeightSpacePolyhedron::createInitialVertices(std::size_t num_objs,
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
    }*/


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

    /*ValueType WeightSpacePolyhedron::getMarkedVertexWOV(std::size_t unit_weight_index) const {
        assert (marked_vertices_.at(unit_weight_index)->hasUnitWeight(unit_weight_index));
        return marked_vertices_[unit_weight_index]->getCurrentWOV();
    }*/

    /*void WeightSpacePolyhedron::updateInitialWSP(SCIP* scip, std::size_t unit_weight_index, const OutcomeType& outcome, bool outcome_is_ray) {
        for (std::size_t i=0; i<marked_vertices_.size(); ++i)
            assert (marked_vertices_[i]->hasUnitWeight(i));

        auto new_vertices = vector<WeightSpaceVertex *> {};
        auto new_edges = vector<pair<WeightSpaceVertex *, Node> > {}; // pair: new vertex, node of adjacent vertex
        auto obsolete = computeObsoleteVertices(scip, marked_vertices_, outcome, outcome_is_ray);
        for (auto obs_vertex : obsolete) {
            auto obs_node = getNode(obs_vertex);
            for (Graph::IncEdgeIt edge(skeleton_, getNode(obs_vertex)); edge != lemon::INVALID; ++edge) {
                auto adj_node = skeleton_.oppositeNode(getNode(obs_vertex), edge);
                auto adj_vertex = getVertex(adj_node);
                auto adj_is_obsolete = isVertexObsolete(scip, adj_vertex, outcome, outcome_is_ray);
            }
        }

        addToSkeleton(new_vertices, new_edges);
    }*/

    bool WeightSpacePolyhedron::areAdjacent(const WeightSpaceVertex* v, const WeightSpaceVertex* w) const {
        auto inc_facets = FacetContainer {};
        assert (v->incident_facets_.size() == w->incident_facets_.size());
        std::set_intersection(v->incident_facets_.cbegin(),
                              v->incident_facets_.cend(),
                              w->incident_facets_.cbegin(),
                              w->incident_facets_.cend(),
                              std::back_inserter(inc_facets),
                              WeightSpaceFacet::compare_facet_ptr);
        return (inc_facets.size() + 1) == w->incident_facets_.size();
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
                if (!isVertexObsolete(scip, adj_vertex, outcome, outcome_is_ray)) {
                    std::cout << "obs vertex: ";
                    obs_vertex->print(std::cout, false);
                    std::cout << "nonobs vertex :";
                    adj_vertex->print(std::cout);
                    auto new_vertex = new WeightSpaceVertex(scip, obs_vertex, adj_vertex, outcome,
                                                            outcome_is_ray); // new vertex between obsolete and non-obsolete verts
                    std::cout << "new vertex: ";
                    new_vertex->print(std::cout);
                    std::cout << "\n";
                    unmarked_vertices_.push_back(new_vertex); // add to unmarked vertices
                    new_vertices.push_back(new_vertex);
                    new_edges.push_back({new_vertex, adj_node});
                }
            }
        }
        addToSkeleton(new_vertices, new_edges); // incorporate new vertices into 1-skeleton
        addCliqueEdgesToSkeleton(new_vertices);
        for (auto obs_vertex : obsolete_vertices) {
            setStatusToObsolete(obs_vertex);
            deleteFromSkeleton(obs_vertex);
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





    vector<WeightSpaceVertex*> WeightSpacePolyhedron::computeObsoleteVertices(SCIP* scip,
                                                                              const OutcomeType& outcome,
                                                                              bool outcome_is_ray) {
        auto obsolete = vector<WeightSpaceVertex*> {curr_investigated_vertex_};
        auto unscanned = std::forward_list<WeightSpaceVertex*> {curr_investigated_vertex_};
        while (!unscanned.empty()) {
            auto obs_vertex = unscanned.front();
            unscanned.pop_front();
            for (Graph::IncEdgeIt edge(skeleton_, getNode(obs_vertex)); edge != lemon::INVALID; ++edge) {
                auto adj_node = skeleton_.oppositeNode(getNode(obs_vertex), edge);
                auto adj_vertex = getVertex(adj_node);
                auto adj_is_obsolete = isVertexObsolete(scip, adj_vertex, outcome, outcome_is_ray);
                if (adj_is_obsolete && std::find(begin(obsolete), end(obsolete), adj_vertex) == end(obsolete)) {
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

    std::vector<WeightSpaceVertex *> WeightSpacePolyhedron::computeObsoleteVerticesCompleteLoop(SCIP *scip,
                                                                                                const OutcomeType &outcome,
                                                                                                bool outcome_is_ray) {
        auto obsolete = std::vector<WeightSpaceVertex *> {curr_investigated_vertex_};
        for (const auto& mv : marked_vertices_)
            assert (!isVertexObsolete(scip, mv, outcome, outcome_is_ray));
        for (const auto& uv : unmarked_vertices_) {
            if (isVertexObsolete(scip, uv, outcome, outcome_is_ray))
                obsolete.push_back(uv);
        }
        return obsolete;
    };

    void WeightSpacePolyhedron::incorporateNewOutcome(SCIP* scip,
                                                      bool completeLoopForObsolete,
                                                      const WeightType& used_weight,
                                                      const OutcomeType &outcome,
                                                      bool outcome_is_ray) {
        assert (curr_investigated_vertex_ != nullptr);
        assert (curr_investigated_vertex_->hasSameWeight(used_weight));
        /*auto obsolete = completeLoopForObsolete ? computeObsoleteVerticesCompleteLoop(scip, outcome, outcome_is_ray) :
                        computeObsoleteVertices(scip, outcome, outcome_is_ray);*/
        auto obsolete = computeObsoleteVerticesCompleteLoop(scip, outcome, outcome_is_ray);
        //auto obsolete = computeObsoleteVertices(scip, outcome, outcome_is_ray);
        //global::print(obsolete, "obsComplete = ");
        //global::print(obsolete2, "obsNoncomplete = ");
        //assert (obsolete.size() == obsolete2.size());
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


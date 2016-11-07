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
#include <iomanip> // std::setprecision
#include <iostream>
#include <forward_list>
#include <functional> // std::negate
#include <iterator>
#include <memory> // std::make_shared
#include <numeric> // std::inner_product
#include <ostream>
#include <unordered_map>
#include <utility> // std::move, std::pair
#include <vector>

#include "double_description_method.h"
#include "global_functions.h"
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

    WeightSpacePolyhedron::WeightSpacePolyhedron(SCIP* scip,
                                                 size_t dimension,
                                                 V_RepC v_rep,
                                                 H_RepC h_rep)
            : wsp_dimension_(dimension),
              curr_investigated_vertex_(nullptr),
              skeleton_(),
              nodes_to_vertices_(skeleton_) {
        assert (!v_rep.empty());
        assert (!h_rep.empty());
        //createInitialFacets(h_rep);
        createInitialVerticesAndSkeleton(scip, std::move(h_rep), std::move(v_rep));
    }

    /*void WeightSpacePolyhedron::createInitialFacets(polytoperepresentation::H_RepC h_rep) {
        for (auto& h : h_rep) {
            facets_.emplace_back(make_shared<WeightSpaceFacet>(h.first, h.second));  //todo check for std::move in computeWeightSpaceResults...
        }
    }*/

    void WeightSpacePolyhedron::createInitialVerticesAndSkeleton(SCIP* scip,
                                                                 H_RepC h_rep,
                                                                 V_RepC v_rep) {
        auto initial_facets = FacetContainer {};
        for (auto& h : h_rep) {
            initial_facets.emplace_back(make_shared<WeightSpaceFacet>(std::move(h.first), std::move(h.second)));  //todo check for std::move in computeWeightSpaceResults...
        }
        for (auto& v : v_rep) {
            //if (v->hasNonZeroWeight()) { // do not consider (0,0,...,0 | -1)
                auto inc_facets = getIncidentFacets(*v, initial_facets);
                auto vertex = new WeightSpaceVertex(inc_facets,
                                                    v->moveWeight(),
                                                    v->getWov());
                auto node = skeleton_.addNode();
                nodes_to_vertices_[node] = vertex;
                vertices_to_nodes_.insert({vertex, node});
                if (vertex->hasUnitWeight()) {
                    vertex->setStatus(WeightSpaceVertex::VertexStatus::marked);
                    marked_vertices_.push_back(vertex);
                }
                else if (vertex->hasZeroWeight()) {
                    vertex->setStatus(WeightSpaceVertex::VertexStatus::special);
                }
                else {
                    unmarked_vertices_.push_back(vertex);
                }
            }
        //}
        for (Graph::NodeIt node(skeleton_); node != lemon::INVALID; ++node) {
            for (Graph::NodeIt succ_node(node); succ_node != lemon::INVALID; ++succ_node) {
                if (succ_node != node && areAdjacent(getVertex(node), getVertex(succ_node)))
                    skeleton_.addEdge(node, succ_node);
            }
        }
    }

    WeightSpacePolyhedron::FacetContainer WeightSpacePolyhedron::getIncidentFacets(const V_RepT& v,
                                                                                   const FacetContainer& initial_facets) const {
        auto facets = FacetContainer{};
        for (size_t i=0; i<initial_facets.size(); i++) {
            if (v.isZeroSlackIndex(i))
                facets.push_back(initial_facets[i]);
        }
        return facets;
    }

    /*vector<pair<OutcomeType, OutcomeType>> WeightSpacePolyhedron::getConstraintsForUnsupported() const {
        auto pairs = vector<pair<OutcomeType,OutcomeType >>{};
        for (auto v : marked_vertices_) {
            if (v->hasUnitWeight())
                continue;
            v->print(std::cout, true);
            pairs.emplace_back(v->getIncFacetsUpperBounds(), v->getIncFacetsLowerBounds());
        }
        return pairs;
    }*/

    /*WeightSpacePolyhedron::FacetContainer WeightSpacePolyhedron::computeIncidentFacets(SCIP* scip,
                                                                                       const FacetContainer& initial_facets,
                                                                                       const polytoperepresentation::V_RepT& v) const {
        auto facets = FacetContainer {};
        for (const auto& f : initial_facets) {
            if (SCIPisEQ(scip, f->getWeightedWeight(v.getWeight()), f->getWOVCoeff()*v.getWov()))
                facets.push_back(f);
        }
        return facets;
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
        assert (curr_investigated_vertex_ != nullptr);
        return curr_investigated_vertex_->getWeight();
    }

    ValueType WeightSpacePolyhedron::getUntestedVertexWOV(const WeightType& untested_weight) const {
        assert (curr_investigated_vertex_ != nullptr);
        assert (curr_investigated_vertex_->hasSameWeight(untested_weight));
        return curr_investigated_vertex_->getCurrentWOV();
    }


    void WeightSpacePolyhedron::addEdgesOfAdjacentVerticesToSkeleton(const vector<WeightSpaceVertex*>& new_vertices,
                                                                     const vector<WeightSpaceVertex*>& zero_vertices) {

        if (!new_vertices.empty()) {
            for (auto it = new_vertices.cbegin(); it != std::prev(new_vertices.cend()); ++it) {
                for (auto succ_it = std::next(it); succ_it != new_vertices.cend(); ++succ_it) {
                    if (areAdjacent(*it, *succ_it))
                        skeleton_.addEdge(getNode(*it), getNode(*succ_it));
                }
                for (auto zero = zero_vertices.cbegin(); zero!=zero_vertices.cend(); ++zero) {
                    if (areAdjacent(*it, *zero))
                        skeleton_.addEdge(getNode(*it), getNode(*zero));
                }
            }
        }
    }

    void WeightSpacePolyhedron::addToSkeleton(const vector<WeightSpaceVertex*>& new_vertices,
                                              const vector<pair<WeightSpaceVertex*, WeightSpaceVertex*> >& new_edges) {
        for (const auto& v : new_vertices) {  // add new nodes
            auto new_node = skeleton_.addNode();
            nodes_to_vertices_[new_node] = v;
            vertices_to_nodes_.insert({v, new_node});
        }
        for (const auto& edge : new_edges) { // add new edges
            skeleton_.addEdge(getNode(edge.first), getNode(edge.second));
        }
    }


    void WeightSpacePolyhedron::deleteFromSkeleton(WeightSpaceVertex* v) {
        skeleton_.erase(getNode(v));
        auto ret = vertices_to_nodes_.erase(v);
        assert (ret == 1);
    }

    void WeightSpacePolyhedron::removeFrom(WeightSpaceVertex* v) {
        if (v != curr_investigated_vertex_) {
            auto v_status = v->getStatus();
            if (v_status == WeightSpaceVertex::VertexStatus::unmarked) {
                auto pos = std::find(begin(unmarked_vertices_),
                                     end(unmarked_vertices_),
                                     v);
                if (pos == end(unmarked_vertices_)) {
                    v->print(std::cout, true);
                    throw std::runtime_error("Weight space vertex expected in Unmarked Vertex Container.\n");
                }
                    unmarked_vertices_.erase(pos);
            }
            else if (v_status == WeightSpaceVertex::VertexStatus::marked) { // v must be weakly-non-dominated point
                auto pos = std::find(begin(marked_vertices_),
                                     end(marked_vertices_),
                                     v);
                if (pos == end(marked_vertices_)) {
                    v->print(std::cout, true);
                    throw std::runtime_error("Weight space vertex expected in Marked Vertex Container.\n");
                }
                marked_vertices_.erase(pos);
            }
            else {
                throw std::runtime_error("Unexpected weight space vertex status.\n");
            }
            v->setStatus(WeightSpaceVertex::VertexStatus::obsolete);
        }
        assert (std::find(obsolete_vertices_.cbegin(),
                          obsolete_vertices_.cend(),
                          v) == obsolete_vertices_.cend());

        obsolete_vertices_.push_back(v);
    }

    bool WeightSpacePolyhedron::areAdjacent(const WeightSpaceVertex* v, const WeightSpaceVertex* w) {
        auto inc_facets = FacetContainer {};
        std::set_intersection(v->incident_facets_.cbegin(),
                              v->incident_facets_.cend(),
                              w->incident_facets_.cbegin(),
                              w->incident_facets_.cend(),
                              std::back_inserter(inc_facets),
                              WeightSpaceFacet::Compare());
        return inc_facets.size() >= wsp_dimension_-1;
    }

    bool WeightSpacePolyhedron::hasValidSkeleton(size_t dim) const {
        for (Graph::NodeIt node(skeleton_); node!=lemon::INVALID; ++node) {
            size_t adj_count{0};
            for (Graph::IncEdgeIt edge(skeleton_, node); edge != lemon::INVALID; ++edge) {
                ++adj_count;
            }
            if (adj_count < dim) {
                auto v = nodes_to_vertices_[node];
                v->print(std::cout, true);
                std::cout << "No of adjacent verts: " << adj_count << " Expected: " << dim << "\n";
                return false;
            }
        }
        return true;
    }

    void WeightSpacePolyhedron::updateWeightSpacePolyhedron(SCIP* scip,
                                                            const OutcomeType& outcome,
                                                            bool outcome_is_ray) {

        auto obs_nonobs_pairs = vector<pair<WeightSpaceVertex*, WeightSpaceVertex*>> {}; // pair: obsolete vertex , adjacent non-obsolete vertex
        assert (SCIPisNegative(scip, curr_investigated_vertex_->computeSlack(outcome, outcome_is_ray)));


        auto tmp_plus = vector<WeightSpaceVertex*> {};
        auto tmp_minus = vector<WeightSpaceVertex*> {};
        auto tmp_zero = vector<WeightSpaceVertex*> {};

        for (const auto& elem : vertices_to_nodes_) {
            if (elem.first->getStatus() == WeightSpaceVertex::VertexStatus::obsolete)
                continue;
            auto facet_partition = getFacetPartition(scip, elem.first, outcome, outcome_is_ray);
            if (facet_partition == FacetPartition::plus) {
                tmp_plus.push_back(elem.first);
            }
            else if (facet_partition == FacetPartition::minus) {
                tmp_minus.push_back(elem.first);
            }
            else if (facet_partition == FacetPartition::zero) {
                tmp_zero.push_back(elem.first);
            }
            else {
                throw std::runtime_error("Unexpected Facet Partition Status\n");
            }
        }
        for (const auto& non_obs : tmp_plus) {
            for (const auto& obs : tmp_minus) {
                if (areAdjacent(non_obs, obs)) {
                    obs_nonobs_pairs.push_back({obs, non_obs});
                }
            }
        }

        curr_investigated_vertex_->setStatus(WeightSpaceVertex::VertexStatus::obsolete);
        /*auto unscanned = std::list<WeightSpaceVertex*> {curr_investigated_vertex_};
        auto obsolete = std::unordered_map<WeightSpaceVertex*, bool>({{curr_investigated_vertex_, true}});
        auto zeros = std::unordered_map<WeightSpaceVertex*, bool>{};*/
        auto new_facet = outcome_is_ray ? make_shared<const WeightSpaceFacet>(outcome, 0) :
                         make_shared<const WeightSpaceFacet>(outcome, 1);
        /*while (!unscanned.empty()) {
            auto obs_vertex = unscanned.front();
            unscanned.pop_front();
            for (Graph::IncEdgeIt edge(skeleton_, getNode(obs_vertex)); edge != lemon::INVALID; ++edge) {
                auto adj_node = skeleton_.oppositeNode(getNode(obs_vertex), edge);
                auto adj_vertex = getVertex(adj_node);
                auto facet_partition = getFacetPartition(scip, adj_vertex, outcome, outcome_is_ray);
                if (facet_partition == FacetPartition::plus && obsolete.count(obs_vertex) == 1) {
                    obs_nonobs_pairs.push_back({obs_vertex, adj_vertex});
                }
                else if (facet_partition == FacetPartition::minus) {
                    if (obsolete.count(adj_vertex) == 0) {
                        obsolete[adj_vertex] = true;
                        unscanned.push_back(adj_vertex);
                    }
                }
                else if (facet_partition == FacetPartition::zero) {
                    if (zeros.count(adj_vertex) == 0) {
                        zeros[adj_vertex] = true;
                        unscanned.push_back(adj_vertex);
                    }
                }
                else {
                    std::runtime_error("Unexpected Facet Partiton Status\n");
                }
            }
        }*/

       /* std::cout << "size of obs_nonobs_pairs: " << obs_nonobs_pairs.size() << "\n";
        std::cout << "size of tmp_obs_non_pairs: " << tmp_nonobs_obs.size() << "\n";
        std::cout << "size of zeros: " << zeros.size() << "\n";
        for (const auto& zero : zeros)
            zero.first->print(std::cout, true);
        std::cout << "size of tmp_zeros: " << tmp_zero.size() << "\n";
        for (const auto& zero : tmp_zero)
            zero->print(std::cout, true);
        std::cout << "okay\n";

        assert (zeros.size() == tmp_zero.size());
        assert (obs_nonobs_pairs.size() == tmp_nonobs_obs.size());

        assert (obs_nonobs_pairs.size() + zeros.size() > 0);*/

        auto new_vertices = vector<WeightSpaceVertex*> {};
        auto new_edges = vector< pair<WeightSpaceVertex*, WeightSpaceVertex*> > {}; // pair: new vertex, adjacent vertex
        for (const auto& v_pair : obs_nonobs_pairs) { // create new vertices between obsolete and non-obsolete vertices
            auto obs_vertex = v_pair.first;
            auto non_obs_vertex = v_pair.second;
            double obs_coeff = non_obs_vertex->computeSlack(outcome, outcome_is_ray);
            double non_obs_coeff = obs_vertex->computeSlack(outcome, outcome_is_ray);
            assert (SCIPisPositive(scip,obs_coeff));
            assert (SCIPisNegative(scip, non_obs_coeff));
            auto new_vertex = new WeightSpaceVertex(obs_coeff, non_obs_coeff,
                                                    obs_vertex, non_obs_vertex,
                                                    new_facet, wsp_dimension_);
            unmarked_vertices_.push_back(new_vertex);
            new_vertices.push_back(new_vertex);
            new_edges.push_back({new_vertex, non_obs_vertex});
            /*}
            else {
                throw std::runtime_error("Unexpected convex combination of obsolete and non-obsolete weight space vertices\n");
            }*/
        }

        addToSkeleton(new_vertices, new_edges);
        addEdgesOfAdjacentVerticesToSkeleton(new_vertices, tmp_zero);
        for (auto obs : tmp_minus) {
            removeFrom(obs);
            deleteFromSkeleton(obs);
        }
        std::cout << "NO UNMARKED VERTICES: " << unmarked_vertices_.size() << "\n";
        std::cout << "NO MARKED VERTICES: " << marked_vertices_.size() << "\n";
    }

    double WeightSpacePolyhedron::calculateConvexCombValue(const WeightSpaceVertex* obs,
                                                           const WeightSpaceVertex* non_obs,
                                                           const OutcomeType& outcome,
                                                           bool outcome_is_ray) {
        auto wov_obs = outcome_is_ray ? 0. : obs->getCurrentWOV();
        auto wov_non_obs = outcome_is_ray ? 0. : non_obs->getCurrentWOV();
        double numerator = wov_obs - std::inner_product(begin(obs->weight_),
                                                             end(obs->weight_),
                                                             begin(outcome),
                                                             0.);
        double denominator = numerator - wov_non_obs + std::inner_product(begin(non_obs->weight_),
                                                                               end(non_obs->weight_),
                                                                               begin(outcome),
                                                                               0.);
        assert (denominator != 0.);
        return numerator / denominator;
    }

    WeightSpacePolyhedron::FacetPartition WeightSpacePolyhedron::getFacetPartition(SCIP* scip,
                                                                                   const WeightSpaceVertex *vertex,
                                                                                   const OutcomeType &outcome,
                                                                                   bool outcome_is_ray) const {
        auto result = vertex->computeSlack(outcome, outcome_is_ray);

        if (SCIPisNegative(scip, result)) {
            return FacetPartition::minus;
        }
        else if (SCIPisPositive(scip, result)) {
            return FacetPartition::plus;
        }
        else {
            assert(SCIPisZero(scip, result));
            return FacetPartition::zero;
        }
    }

    void WeightSpacePolyhedron::incorporateNewOutcome(SCIP* scip,
                                                      const WeightType& used_weight,
                                                      const OutcomeType &outcome,
                                                      bool outcome_is_ray) {
        assert (curr_investigated_vertex_ != nullptr);
        assert (curr_investigated_vertex_->hasSameWeight(used_weight));
        updateWeightSpacePolyhedron(scip, outcome, outcome_is_ray);
        resetCurrentInvestigatedVertex(); // requirement: call after updateWeightSpacePolyhedron;
    }


    void WeightSpacePolyhedron::incorporateKnownOutcome(const WeightType& used_weight) {
        assert (curr_investigated_vertex_ != nullptr);
        assert (curr_investigated_vertex_->hasSameWeight(used_weight));
        curr_investigated_vertex_->setStatus(WeightSpaceVertex::VertexStatus::marked);
        marked_vertices_.push_back(curr_investigated_vertex_);
        resetCurrentInvestigatedVertex();
    }

    void WeightSpacePolyhedron::printUnmarkedVertices(ostream &os, bool printFacets) const {
        printVertices(unmarked_vertices_, "UNMARKED VERTICES:", os, printFacets);
    }

    void WeightSpacePolyhedron::printMarkedVertices(ostream &os, bool printFacets) const {
        printVertices(marked_vertices_, "MARKED VERTICES:", os, printFacets);
    }

    void WeightSpacePolyhedron::printObsoleteVertices(ostream &os, bool printFacets) const {
        printVertices(obsolete_vertices_, "OBSOLETE VERTICES:", os, printFacets);
    }


}


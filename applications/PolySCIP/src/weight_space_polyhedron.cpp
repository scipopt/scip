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
              curr_investigated_vertex_(nullptr) {
        assert (!v_rep.empty());
        assert (!h_rep.empty());
        createInitialVerticesAndSkeleton(scip, std::move(h_rep), std::move(v_rep));
    }


    void WeightSpacePolyhedron::createInitialVerticesAndSkeleton(SCIP* scip,
                                                                 H_RepC h_rep,
                                                                 V_RepC v_rep) {
        auto initial_facets = FacetContainer {};
        for (auto& h : h_rep) {
            initial_facets.emplace_back(make_shared<WeightSpaceFacet>(std::move(h.first), std::move(h.second)));  //todo check for std::move in computeWeightSpaceResults...
        }
        for (auto& v : v_rep) {
            auto inc_facets = getIncidentFacets(*v, initial_facets);
            auto vertex = new WeightSpaceVertex(inc_facets,
                                                v->moveWeight(),
                                                v->getWov());

            if (vertex->hasUnitWeight()) {
                vertex->setStatus(WeightSpaceVertex::VertexStatus::marked);
                marked_and_special_vertices_.push_back(vertex);
            }
            else if (vertex->hasZeroWeight()) {
                vertex->setStatus(WeightSpaceVertex::VertexStatus::special);
                marked_and_special_vertices_.push_back(vertex);
            }
            else {
                unmarked_vertices_.push_back(vertex);  // VertexStatus::unmarked
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


    WeightSpacePolyhedron::~WeightSpacePolyhedron() {
        for (auto v_ptr : marked_and_special_vertices_)
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


    void WeightSpacePolyhedron::makeObsolete(WeightSpaceVertex *v) {
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
                auto pos = std::find(begin(marked_and_special_vertices_),
                                     end(marked_and_special_vertices_),
                                     v);
                if (pos == end(marked_and_special_vertices_)) {
                    v->print(std::cout, true);
                    throw std::runtime_error("Weight space vertex expected in Marked Vertex Container.\n");
                }
                marked_and_special_vertices_.erase(pos);
            }
            else {
                throw std::runtime_error("Unexpected weight space vertex status.\n");
            }
        }
        assert (std::find(obsolete_vertices_.cbegin(),
                          obsolete_vertices_.cend(),
                          v) == obsolete_vertices_.cend());
        v->setStatus(WeightSpaceVertex::VertexStatus::obsolete);
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

    void WeightSpacePolyhedron::addVertexToCorrespondingPartition(SCIP* scip,
                                                                  WeightSpaceVertex* vertex,
                                                                  const OutcomeType& outcome,
                                                                  bool outcome_is_ray,
                                                                  vector<WeightSpaceVertex*>& minus_vertices,
                                                                  vector<WeightSpaceVertex*>& plus_vertices,
                                                                  vector<WeightSpaceVertex*>& zero_vertices) const {
        auto result = vertex->computeSlack(outcome, outcome_is_ray);

        if (SCIPisNegative(scip, result)) {
            minus_vertices.push_back(vertex);
        }
        else if (SCIPisPositive(scip, result)) {
            plus_vertices.push_back(vertex);
        }
        else {
            assert(SCIPisZero(scip, result));
            zero_vertices.push_back(vertex);
        }
    }


    void WeightSpacePolyhedron::updateWeightSpacePolyhedron(SCIP* scip,
                                                            const OutcomeType& outcome,
                                                            bool outcome_is_ray) {

        auto obs_nonobs_pairs = vector<pair<WeightSpaceVertex*, WeightSpaceVertex*>> {}; // pair: obsolete vertex , adjacent non-obsolete vertex
        assert (SCIPisNegative(scip, curr_investigated_vertex_->computeSlack(outcome, outcome_is_ray)));


        auto plus = vector<WeightSpaceVertex*> {};
        auto minus = vector<WeightSpaceVertex*> {};
        auto zero = vector<WeightSpaceVertex*> {};

        for (auto vertex : marked_and_special_vertices_) {
            addVertexToCorrespondingPartition(scip,
                                              vertex,
                                              outcome,
                                              outcome_is_ray,
                                              minus,
                                              plus,
                                              zero);
        }
        for (auto vertex : unmarked_vertices_) {
            addVertexToCorrespondingPartition(scip,
                                              vertex,
                                              outcome,
                                              outcome_is_ray,
                                              minus,
                                              plus,
                                              zero);
        }
        minus.push_back(curr_investigated_vertex_);

        for (const auto& non_obs : plus) {
            for (const auto& obs : minus) {
                if (areAdjacent(non_obs, obs)) {
                    obs_nonobs_pairs.push_back({obs, non_obs});
                }
            }
        }


        auto new_facet = outcome_is_ray ? make_shared<const WeightSpaceFacet>(outcome, 0) :
                         make_shared<const WeightSpaceFacet>(outcome, 1);

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
        }

        for (auto vertex : minus) {
            makeObsolete(vertex);
        }

        std::cout << "NO UNMARKED VERTICES: " << unmarked_vertices_.size() << "\n";
        std::cout << "NO MARKED VERTICES: " << marked_and_special_vertices_.size() << "\n";
    }


    void WeightSpacePolyhedron::incorporateNewOutcome(SCIP* scip,
                                                      const WeightType& used_weight,
                                                      const OutcomeType& outcome,
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
        marked_and_special_vertices_.push_back(curr_investigated_vertex_);
        resetCurrentInvestigatedVertex();
    }

    void WeightSpacePolyhedron::printUnmarkedVertices(ostream &os, bool printFacets) const {
        printVertices(unmarked_vertices_, "UNMARKED VERTICES:", os, printFacets);
    }

    void WeightSpacePolyhedron::printMarkedVertices(ostream &os, bool printFacets) const {
        printVertices(marked_and_special_vertices_, "MARKED VERTICES:", os, printFacets);
    }

    void WeightSpacePolyhedron::printObsoleteVertices(ostream &os, bool printFacets) const {
        printVertices(obsolete_vertices_, "OBSOLETE VERTICES:", os, printFacets);
    }


}

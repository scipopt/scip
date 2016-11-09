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

/** @brief  The (partial) weight space polyhedron
 *
 * This class represents the (partial) weight space polyhedron P =
 * {(w,a) \in \Lambda \times R : w \cdot y >= a \forall y \in Y'}
 * where Y' is the set of non-dominated points computed so far and
 * \Lambda is the set of normalized weights
 */

#ifndef POLYSCIP_SRC_WEIGHT_SPACE_POLYHEDRON_H_INCLUDED
#define POLYSCIP_SRC_WEIGHT_SPACE_POLYHEDRON_H_INCLUDED

#include <cstddef>
#include <deque>
#include <list>
#include <iomanip> //std::set_precision
#include <iostream>
#include <iterator>
#include <map>
#include <memory> // std::shared_ptr
#include <ostream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility> // std::pair
#include <vector>

#include "objscip/objscip.h"
#include "polyscip_types.h"
#include "weight_space_facet.h"
#include "double_description_method.h"

namespace polyscip {

    using V_RepT = doubledescription::V_RepT;
    using V_RepC = doubledescription::V_RepC;
    using H_RepC = doubledescription::H_RepC;
    class WeightSpaceVertex;

    /** 1-skeleton of the (partial) weight space polyhedron. */
    class WeightSpacePolyhedron {
    public:
        using FacetContainer = std::vector<std::shared_ptr<const WeightSpaceFacet>>;

        /** Creates the skeleton of the initial (partial) weight
        * space polyhedron P = {(w,a) \in \Lambda \times R : w \cdot y^1
        * >= a}
        * @param num_objs number of objectives of given problem
        * @param point first computed (weakly non-dominated) point
        * @param point_weighted_obj_val weighted objective value of first point
        * @param initial_rays rays which were computed while computing the first non-dominated point
        * @param unit_weight_info pair indicating whether first point was computed by using
        * a unit weight; if unit_weight_info.first is true, then unit_weight_info.second contains the
        * index with value 1; note: first index is 0
        */
        explicit WeightSpacePolyhedron(SCIP* scip,
                                       std::size_t wsp_dimension,
                                       V_RepC v_rep,
                                       H_RepC h_rep);

        /** Destructor */
        ~WeightSpacePolyhedron();

        /** Checks whether there is an unmarked weight space vertex with an untested weight
         *  @return true if there is an unmarked weight space vertex with untested weight; false otherwise
         */
        bool hasUntestedWeight() const;

        /** Returns an untested weight; Should only be called after hasUnmarkedVertex()
         * returned true
         *  @return an untested weight
         */
        WeightType getUntestedWeight();

        ValueType getUntestedVertexWOV(const WeightType& untested_weight) const;

        /** Incorporates an newly found unbounded non-dominated ray
         * into the (partial) weight space polyhedron
         //todo
         */
        void incorporateNewOutcome(SCIP* scip,
                                   const WeightType& used_weight,
                                   const OutcomeType& outcome,
                                   bool outcome_is_ray = false);

        void incorporateKnownOutcome(const WeightType& weight);

        /** Prints unmarked vertices to output stream
         *  @param printFacets if true, facet information of unmarked vertices is also printed
         */
        void printUnmarkedVertices(std::ostream& os = std::cout, bool printFacets = false) const;

        /** Prints marked vertices to standard output
         *  @param printFacets if true, facet information of marked vertices is also printed
         */
        void printMarkedVertices(std::ostream& os = std::cout, bool printFacets = false) const;

        /** Prints obsolete vertices to output stream
         * @param printFacets if true, facet information of obsolete vertices is also printed
         */
        void printObsoleteVertices(std::ostream& os = std::cout, bool printFacets = false) const;

        double calculateConvexCombValue(const WeightSpaceVertex* obsolete,
                                        const WeightSpaceVertex* non_obsolete,
                                        const OutcomeType& outcome,
                                        bool outcome_is_ray);

        bool areAdjacent(const WeightSpaceVertex* v, const WeightSpaceVertex* w);

    private:

        using MarkedVertexContainer = std::vector<WeightSpaceVertex*>;
        /** Container used to store the unmarked weight space vertices
         *  Needs to support: empty(), size(), valid iterators of objects after erasing another object
         */
        using UnmarkedVertexContainer = std::list<WeightSpaceVertex*>;


        using ObsoleteVertexContainer = std::vector<WeightSpaceVertex*>;


        void createInitialFacets(H_RepC h_rep) = delete;

        void createInitialVerticesAndSkeleton(SCIP* scip,
                                              H_RepC h_rep,
                                              V_RepC v_rep);

        FacetContainer computeIncidentFacets(SCIP* scip,
                                             const FacetContainer& initial_facets,
                                             const V_RepT& v) const = delete;

        FacetContainer getIncidentFacets(const V_RepT& v,
                                         const FacetContainer& initial_facets) const;

        void updateWeightSpacePolyhedron(SCIP* scip,
                                         const OutcomeType& outcome,
                                         bool outcome_is_ray);


        void resetCurrentInvestigatedVertex() {
            curr_investigated_vertex_ = nullptr;
        }

        void addVertexToCorrespondingPartition(SCIP* scip,
                                               WeightSpaceVertex* vertex,
                                               const OutcomeType& outcome,
                                               bool outcome_is_ray,
                                               std::vector<WeightSpaceVertex*>& minus_vertices,
                                               std::vector<WeightSpaceVertex*>& plus_vertices,
                                               std::vector<WeightSpaceVertex*>& zero_vertices) const;

        void makeObsolete(WeightSpaceVertex *v);

        /** Template function to print vertices; is used by public print{Marked,Obsolete,Unmarked}Vertices
         * functions
         * @param container Container of vertices which are to be printed
         * @param description Initial string describing what is going to be printed
         * @param os output stream
         * @param printFacets if true, also facet information of vertices is printed
         */
        template <typename Container>
        void printVertices(const Container& container,
                           std::string description,
                           std::ostream& os,
                           bool printFacets) const;

        std::size_t wsp_dimension_;

        /**< all marked weight space vertices */
        MarkedVertexContainer marked_and_special_vertices_;
        /**< all unmarked weight space vertices */
        UnmarkedVertexContainer unmarked_vertices_;
        /**< all obsolete weight space vertices */
        ObsoleteVertexContainer obsolete_vertices_;
        /**< unmarked weight space vertex currently investigated */
        WeightSpaceVertex* curr_investigated_vertex_;
        /**< 1-skeleton of the weight space polyhedron */
    };



    template <typename Container>
    void WeightSpacePolyhedron::printVertices(const Container& container,
                                              std::string description,
                                              std::ostream &os,
                                              bool printFacets) const {
        os << description << "\n";
        for (const auto& elem : container) {
            elem->print(os, printFacets);
            os << "\n";
        }
        os << "\n";
    }

}

#endif //POLYSCIP_SRC_WEIGHT_SPACE_POLYHEDRON_H_INCLUDED

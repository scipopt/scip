/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*        This file is part of the program PolySCIP                          */
/*                                                                           */
/*    Copyright (C) 2012-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  PolySCIP is distributed under the terms of the ZIB Academic License.     */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with PolySCIP; see the file LICENCE.                               */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**
 * @file weight_space_polyhedron.h
 * @brief Class representing the 1-skeleton of the weight space polyhedron
 * @author Sebastian Schenker
 *
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

    using V_RepT = doubledescription::V_RepT; ///< abbreviation
    using V_RepC = doubledescription::V_RepC; ///< abbreviation
    using H_RepC = doubledescription::H_RepC; ///< abbreviation
    class WeightSpaceVertex;

    /**
     * @class WeightSpacePolyhedron
     * @brief Class representing the 1-skeleton of the weight space polyhedron
     * @details The (partial) weight space polyhedron
     * P = {(w,a) \\in W \\times R : w \\cdot y >= a \\forall y \\in Y'}
     * where Y' is the set of non-dominated extreme points computed so far and
     * W = {w \\in R^k: w_i >= 0, w_1 + ... + w_k = 1} is the set of non-negative normalized weights
     */
    class WeightSpacePolyhedron {
    public:
        using FacetContainer = std::vector<std::shared_ptr<const WeightSpaceFacet>>; ///< Container for facets

        /**
         * Default constructor
         * @param wsp_dimension Dimension of the weight space polyhedron
         * @param v_rep v-representation of the initial polyhedron
         * @param h_rep h-representation of the initial polyhedron
         */
        explicit WeightSpacePolyhedron(std::size_t wsp_dimension,
                                       V_RepC v_rep,
                                       H_RepC h_rep);

        /**
         * Destructor
         */
        ~WeightSpacePolyhedron();

        /**
         * Indicates whether there is an unmarked vertex in the current polyhedron
         * @return true if there is an unmarked vertex; false otherwise
         */
        bool hasUntestedWeight() const;

        /**
         * Get corresponding weight of unmarked vertex of the polyhedron
         * @attention Should only be called after hasUnmarkedVertex() returned true
         * @return Weight vector
         */
        WeightType getUntestedWeight();

        /**
         * Get weighted objective value of unmarked vertex
         * @param untested_weight Corresponding weight of currently investigated vertex
         * @return Weighted objective value of currently investigated vertex
         */
        ValueType getUntestedVertexWOV(const WeightType& untested_weight) const;

        /**
         * Incorporate new outcome (yielding new vertex) into weight space polyhedron
         * @param scip SCIP pointer
         * @param used_weight Weight used to compute outcome
         * @param outcome Newly computed outcome
         * @param outcome_is_ray Indicates whether newly computed outcome is a ray
         */
        void incorporateNewOutcome(SCIP* scip,
                                   const WeightType& used_weight,
                                   const OutcomeType& outcome,
                                   bool outcome_is_ray = false);

        /**
         * Incorporate weight not yielding new vertex into weight space polyhedron
         * @param used_weight Used weight vector
         */
        void incorporateKnownOutcome(const WeightType& used_weight);

        /**
         * Print function for unmarked vertices
         * @param os Output stream to write to
         * @param printFacets Indicates whether corresponding facets should be printed
         */
        void printUnmarkedVertices(std::ostream& os = std::cout,
                                   bool printFacets = false) const;

        /**
         * Print function for marked vertices
         * @param os Output stream to write to
         * @param printFacets Indicates whether corresponding facets should be printed
         */
        void printMarkedVertices(std::ostream& os = std::cout,
                                 bool printFacets = false) const;

        /**
         * Print function for obsolete vertices
         * @param os Output stream to write to
         * @param printFacets Indicates whether corresponding facets should be printed
         */
        void printObsoleteVertices(std::ostream& os = std::cout,
                                   bool printFacets = false) const;

        /**
         * Indicates whether two weight space vertices are adjacent in weight space polyhedron
         * @param v First vertex
         * @param w Second vertex
         * @return true if first vertex is adjacent to second vertex; false otherwise
         * @todo Const qualification
         */
        bool areAdjacent(const WeightSpaceVertex* v,
                         const WeightSpaceVertex* w);

    private:

        using MarkedVertexContainer = std::vector<WeightSpaceVertex*>; ///< Container for marked vertices @attention needs to support empty(), size(), valid iterators after erasing objects
        using UnmarkedVertexContainer = std::list<WeightSpaceVertex*>; ///< Container for unmarked vertices
        using ObsoleteVertexContainer = std::vector<WeightSpaceVertex*>; ///< Container for obsolete vertices

        /**
         * Create initial weight space vertices and corresponding 1-skeleton of weight space polyhedron
         * @param h_rep h-representation of initial weight space polyhedron
         * @param v_rep v-representation of initial weight space polyhedron
         */
        void createInitialVerticesAndSkeleton(H_RepC h_rep,
                                              V_RepC v_rep);

        /**
         * Get incident facets
         * @param v Corresponding vertex
         * @param initial_facets Facets
         * @return Incident facets of given vertex with respect to given facets
         */
        FacetContainer getIncidentFacets(const V_RepT& v,
                                         const FacetContainer& initial_facets) const;

        /**
         * Update 1-skeleton of weight space polyhedron
         * @param scip SCIP pointer
         * @param outcome New outcome yielding a new vertex in weight space polyhedron
         * @param outcome_is_ray Indicates whether outcome corresponds to ray
         */
        void updateWeightSpacePolyhedron(SCIP* scip,
                                         const OutcomeType& outcome,
                                         bool outcome_is_ray);

        /**
         * Nullify currently investigated vertex
         */
        void resetCurrentInvestigatedVertex() {
            curr_investigated_vertex_ = nullptr;
        }

        /**
         * Add vertex to corresponding partition (minus, zero, plus) with respect to the given outcome
         * @param scip SCIP pointer
         * @param vertex Weight space vertex
         * @param outcome Outcome
         * @param outcome_is_ray Indicates whether given outcome is a ray
         * @param minus_vertices Container storing vertices in minus partition
         * @param plus_vertices Container storing vertices in plus partition
         * @param zero_vertices Container storing vertices in zero partition
         */
        void addVertexToCorrespondingPartition(SCIP* scip,
                                               WeightSpaceVertex* vertex,
                                               const OutcomeType& outcome,
                                               bool outcome_is_ray,
                                               std::vector<WeightSpaceVertex*>& minus_vertices,
                                               std::vector<WeightSpaceVertex*>& plus_vertices,
                                               std::vector<WeightSpaceVertex*>& zero_vertices) const;

        /**
         * Set status of given vertex to obsolete status
         * @param v Weight space vertex
         */
        void makeObsolete(WeightSpaceVertex *v);


        /**
         * Print function for vertex container
         * @tparam Container Container type
         * @param container Container with elements to be printed
         * @param description Description to be printed before elements
         * @param os Output stream to write to
         * @param printFacets Indicates whether corresponding facets should be printed
         */
        template <typename Container>
        void printVertices(const Container& container,
                           std::string description,
                           std::ostream& os,
                           bool printFacets) const;

        std::size_t wsp_dimension_; ///< Corresponding dimension of weight space polyhedron

        MarkedVertexContainer marked_and_special_vertices_; ///< Marked weight space vertices
        UnmarkedVertexContainer unmarked_vertices_; ///< Unmarked weight spaces vertices
        ObsoleteVertexContainer obsolete_vertices_; ///< Obsolete weight space vertices
        WeightSpaceVertex* curr_investigated_vertex_; ///< Weight space vertex currently investigated
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

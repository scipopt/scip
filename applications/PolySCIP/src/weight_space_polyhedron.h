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


#include <iostream>
#include <list>
#include <memory> // std::shared_ptr
#include <ostream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility> // std::pair
#include <vector>

#undef GCC_VERSION /* lemon/core.h redefines GCC_VERSION additionally to scip/def.h */
#include "lemon/list_graph.h"

#include "polyscip.h"
#include "weight_space_facet.h"
#include "weight_space_vertex.h"

namespace polyscip {

    /** 1-skeleton of the (partial) weight space polyhedron. */
    class WeightSpacePolyhedron {
    public:

        /** Creates the skeleton of the initial (partial) weight
        * space polyhedron P = {(w,a) \in \Lambda \times R : w \cdot y^1
        * >= a}
        * @param num_objs number of objectives of given problem
        * @param point first computed (weakly non-dominated) point
        * @param point_weighted_obj_val weighted objective value of first point
        * @param unit_weight_info pair indicating whether first point was computed by using
        * a unit weight; if unit_weight_info.first is true, then unit_weight_info.second contains the
        * index with value 1; note: first index is 0
        * @param initial_rays rays which were computed while computing the first non-dominated point
        */
        explicit WeightSpacePolyhedron(unsigned num_objs,
                                       Polyscip::OutcomeType point,
                                       Polyscip::ValueType point_weighted_obj_val,
                                       std::pair<bool, Polyscip::WeightType::size_type> unit_weight_info,
                                       const Polyscip::ResultContainer& initial_rays);

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
        Polyscip::WeightType getUntestedWeight();

        /** Incorporates an newly found unbounded non-dominated ray
         * into the (partial) weight space polyhedron
         * @param new_ray the newly found non-dominated ray that was
         * computed by considering the weight and weighted objective
         * value given by old_vertex
         *  @param old_vertex the vertex (yielding weight and weight
         * objective value) that was considered in last computation
         */
        void incorporateOutcome(const Polyscip::OutcomeType& outcome,
                                const Polyscip::WeightType& weight,
                                bool outcome_is_ray);

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

    private:
        using Graph = lemon::ListGraph;
        using Node = Graph::Node;
        /** Container used to store the marked weight space vertices
         *  Needs to support: push_back, size()
         */
        using MarkedVertexContainer = std::vector<WeightSpaceVertex*>;
        /** Container used to store the unmarked weight space vertices
         *  Needs to support: empty(), size(), valid iterators of objects after erasing another object
         */
        using UnmarkedVertexContainer = std::list<WeightSpaceVertex*>;

        using ObsoleteVertexContainer = std::vector<WeightSpaceVertex*>;

        using NodeMap = Graph::NodeMap<WeightSpaceVertex*>;
        using VertexMap = std::unordered_map<WeightSpaceVertex*, Node>;

        /** Checks whether outcome is newly found (non-dominated) outcome, i.e,
         * whether weighted objective value of outcome with respect to the weight of the currently
         * considered weight space vertex is less than previously known weighted objective value
         * @param outcome computed outcome
         * @param outcome_is_ray true if computed outcome corresponds to unbounded ray; false otherwise
         */
        bool isNewOutcome(const Polyscip::OutcomeType& outcome, bool outcome_is_ray);

        /** Incorporates a newly found (non-dominated) outcome [point or ray] into the partial
         * weight space polyhedron; is only called after fundction isNewOutcome returned true
         * @param outcome newly found point or ray
         * @param outcome_is_ray true if computed outcome is an unbounded ray
         */
        void incorporateNewOutcome(const Polyscip::OutcomeType& outcome, bool outcome_is_ray);

        /** Incorporates an already known outcome [point], i.e, the current investigated vertex is marked;
         * is only called after function isNewOutcome returned false
         * @param outcome computed outcome
         */
        void incorporateOldOutcome(const Polyscip::OutcomeType& outcome);

        /** Creates initial weight space vertices
         *  @param num_objs number of objectives of given problem
         *  @param point first computed (weakly non-dominated) point
         *	@param weighted_obj_val weighted objective value of given point
         *  @param boundary_facets initial boundary facets of the weight space polyhedron
        */
        void createInitialVertices(unsigned num_objs,
                                   Polyscip::OutcomeType point,
                                   Polyscip::ValueType weighted_obj_val,
                                   WeightSpaceVertex::FacetContainer boundary_facets);

        /** Creates initial 1-skeleton of complete graph with number of
         * objectives many vertices
         */
        void createInitialSkeleton();

        /** Makes unmarked vertex with unit weight (1 in unit weight is at
         *  unit_weight_index) an marked vertex
         *  @param unit_weight_index index of 1 in unit weight
         */
        void setMarkedVertex(Polyscip::WeightType::size_type unit_weight_index);

        /** The initial weight space vertex v* having weight which coincides with
         * given ray_weight is changed, i.e., for each adjacent vertex of v* a new vertex/node with weight
         * slightly leaning towards the adjacent vertex is added and connected via an edge in the
         * 1-Skeleton; then edges among all new nodes are added in the 1-Skeleton such that the subgraph
         * of the new nodes is a complete graph; finally the initial weight space vertex v* is made
         * obsolete and its corresponding node is deleted from the 1-Skeleton graph
         * @param ray computed ray
         * @param ray_weight associated weight which yields ray
         * @return true if initial weight space vertex with weight coinciding with ray_weight was found
         * (and the weight space polyhedron changed); false otherwiseo
         */
        bool updateInitialWeightSpacePolyhedron(const Polyscip::OutcomeType& ray,
                                                const Polyscip::WeightType& ray_weight);

        /** Returns weight w that fulfills the following equations:
         * 1) w = h1 * weight1 + h2 * weight2 [with h1,h2 >= 0 and h1 + h2 = 1]
         * 2) w \cdot outcome = 0
         * weight w is calculated by insertion of 1) into 2)
         * @param weight1 weight of vertex
         * @param weight2 weight of another vertex
         * @param outcome computed outcome
         * @param h1_shift_value //TODO explaination
         * @return convex combination w of weight1 and weight2 fulfilling w \cdot outcome = 0
         */
        Polyscip::WeightType calculateWeight(const Polyscip::WeightType& weight1,
                                             const Polyscip::WeightType& weight2,
                                             const Polyscip::OutcomeType& outcome,
                                             double h1_shift_value = 0.);

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

        /**< all marked weight space vertices */
        MarkedVertexContainer marked_vertices_;
        /**< all unmarked weight space vertices  */
        UnmarkedVertexContainer unmarked_vertices_;
        ObsoleteVertexContainer obsolete_vertices_;
        /**< unmarked weight space vertex currently investigated */
        WeightSpaceVertex* curr_investigated_vertex_;
        /**< 1-skeleton of the weight space polyhedron */
        Graph skeleton_;
        /**< maps nodes to vertices */
        NodeMap nodes_to_vertices_;
        /**< maps vertices to nodes */
        VertexMap vertices_to_nodes_;

    };

    template <typename Container>
    void WeightSpacePolyhedron::printVertices(const Container& container,
                                              std::string description,
                                              std::ostream &os,
                                              bool printFacets) const {
        os << description << "\n";
        for (const auto& elem : container)
            elem->print(os, printFacets);
        os << "\n";
    }

}

#endif // POLYSCIP_SRC_WEIGHT_SPACE_POLYHEDRON_H_INCLUDED 

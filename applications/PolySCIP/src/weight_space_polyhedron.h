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

#undef GCC_VERSION /* lemon/core.h redefines GCC_VERSION additionally to scip/def.h */
#include "lemon/list_graph.h"
#include "objscip/objscip.h"
#include "polyscip_types.h"
#include "weight_space_facet.h"

namespace polyscip {

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
        explicit WeightSpacePolyhedron(std::size_t num_objs, const OutcomeType& point);

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

        ValueType getMarkedVertexWOV(std::size_t unit_weight_index) const;

        /** Incorporates an newly found unbounded non-dominated ray
         * into the (partial) weight space polyhedron
         //todo
         */
        void incorporateNewOutcome(SCIP* scip,
                                   bool completeLoopForObsolete,
                                   const WeightType& used_weight,
                                   const OutcomeType& outcome,
                                   bool outcome_is_ray = false);

        void incorporateKnownOutcome(const WeightType& weight);

        /** Makes unmarked vertex with unit weight (1 in unit weight is at
         *  unit_weight_index) an marked vertex
         *  @param unit_weight_index index of 1 in unit weight
         */
        void markVertex(std::size_t unit_weight_index) = delete;

        //TODO adjust documentation
        /** The initial weight space vertex v* having weight which coincides with
         * given ray_weight is changed, i.e., for each adjacent vertex of v* a new vertex/node with weight
         * slightly leaning towards the adjacent vertex is added and connected via an edge in the
         * 1-Skeleton; then edges among all new nodes are added in the 1-Skeleton such that the subgraph
         * of the new nodes is a complete graph; finally the initial weight space vertex v* is made
         * obsolete and its corresponding node is deleted from the 1-Skeleton graph
         * @param ray computed ray
         * @return true if weight space polyhedron changed
         */
        void updateInitialWSP(SCIP* scip, std::size_t unit_weight_index, const OutcomeType& outcome, bool outcome_is_ray = false);

        void addCliqueEdgesToSkeleton() {addCliqueEdgesToSkeleton(unmarked_vertices_);};

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



        std::size_t getNumberOfGraphNodes() const;

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

        using WeightMap = std::map<OutcomeType, WeightType>;

        /** Checks whether outcome is newly found (non-dominated) outcome, i.e,
         * whether weighted objective value of outcome with respect to the weight of the currently
         * considered weight space vertex is less than previously known weighted objective value
         * @param outcome computed outcome
         * @param outcome_is_ray true if computed outcome corresponds to unbounded ray; false otherwise
         */
        //bool isNewOutcome(const OutcomeType& outcome, bool outcome_is_ray);

        bool isVertexObsolete(SCIP* scip,
                              const WeightSpaceVertex* vertex,
                              const OutcomeType& outcome,
                              bool outcome_is_ray = false);

        /** Creates initial weight space vertices
         *  @param num_objs number of objectives of given problem
         *  @param point first computed (weakly non-dominated) point
         *	@param weighted_obj_val weighted objective value of given point
         *  @param boundary_facets initial boundary facets of the weight space polyhedron
        */
        void createInitialVertices(std::size_t num_objs,
                                   const OutcomeType& point,
                                   FacetContainer boundary_facets);

        /** Creates initial 1-skeleton of complete graph with number of
         * objectives many vertices
         */
        void createInitialSkeleton();

        void updateWeightSpacePolyhedron(SCIP* scip, const std::vector<WeightSpaceVertex*>& obsolete_vertices,
                                         const OutcomeType& outcome,
                                         bool outcome_is_ray);

        template <typename Container>
        std::vector<WeightSpaceVertex*> computeObsoleteVertices(SCIP* scip,
                                                                const Container& container,
                                                                const OutcomeType& outcome,
                                                                bool outcome_is_ray);

        std::vector<WeightSpaceVertex*> computeObsoleteVertices(SCIP* scip,
                                                                const OutcomeType& outcome,
                                                                bool outcome_is_ray);

        void resetCurrentInvestigatedVertex() {
            curr_investigated_vertex_ = nullptr;
        }

        void addToSkeleton(const std::vector<WeightSpaceVertex*>& new_vertices,
                           const std::vector< std::pair<WeightSpaceVertex*, Node> >& new_edges);

        template <typename Container>
        void addCliqueEdgesToSkeleton(const Container& vertices);

        void deleteFromSkeleton(WeightSpaceVertex* v);

        void setStatusToObsolete(WeightSpaceVertex* v);

        Node getNode(WeightSpaceVertex* vertex) {return vertices_to_nodes_.at(vertex);};

        WeightSpaceVertex* getVertex(Node n) {return nodes_to_vertices_[n];};

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
        /**< all unmarked weight space vertices */
        UnmarkedVertexContainer unmarked_vertices_;
        /**< all obsolete weight space vertices */
        ObsoleteVertexContainer obsolete_vertices_;
        /**< unmarked weight space vertex currently investigated */
        WeightSpaceVertex* curr_investigated_vertex_;
        /**< 1-skeleton of the weight space polyhedron */
        Graph skeleton_;
        /**< maps nodes to vertices */
        NodeMap nodes_to_vertices_;
        /**< maps vertices to nodes */
        VertexMap vertices_to_nodes_;
        /**< maps outcomes to weights */
        WeightMap outcomes_to_weights_;

    };

    template <typename Container>
    std::vector<WeightSpaceVertex*> WeightSpacePolyhedron::computeObsoleteVertices(SCIP* scip,
                                                                                   const Container& container,
                                                                                   const OutcomeType& outcome,
                                                                                   bool outcome_is_ray) {
        auto obsolete = std::vector<WeightSpaceVertex*> {};
        for (auto vertex : container) {
            if (isVertexObsolete(scip, vertex, outcome, outcome_is_ray)) {
                obsolete.push_back(vertex);
            }
        }
        return obsolete;
    };

    template <typename Container>
    void WeightSpacePolyhedron::addCliqueEdgesToSkeleton(const Container& vertices) {
        for (auto it = std::begin(vertices); it!=std::end(vertices); ++it ) {
            for (auto succ_it = std::begin(vertices); succ_it!=std::end(vertices); ++succ_it) {
                skeleton_.addEdge(getNode(*it), getNode(*succ_it));
            }
        }
    };

    template <typename Container>
    void WeightSpacePolyhedron::printVertices(const Container& container,
                                              std::string description,
                                              std::ostream &os,
                                              bool printFacets) const {
        os << std::setprecision(9) << description << "\n";
        for (const auto& elem : container)
            elem->print(os, printFacets);
        os << "\n";
    };

}

#endif //POLYSCIP_SRC_WEIGHT_SPACE_POLYHEDRON_H_INCLUDED

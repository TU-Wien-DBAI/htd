/* 
 * File:   BucketEliminationGraphDecompositionAlgorithm.hpp
 * 
 * Author: ABSEHER Michael (abseher@dbai.tuwien.ac.at)
 * 
 * Copyright 2015-2016, Michael Abseher
 *    E-Mail: <abseher@dbai.tuwien.ac.at>
 * 
 * This file is part of htd.
 * 
 * htd is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free 
 * Software Foundation, either version 3 of the License, or (at your 
 * option) any later version.
 * 
 * htd is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public 
 * License for more details.

 * You should have received a copy of the GNU General Public License
 * along with htd.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HTD_HTD_BUCKETELIMINATIONGRAPHDECOMPOSITIONALGORITHM_HPP
#define	HTD_HTD_BUCKETELIMINATIONGRAPHDECOMPOSITIONALGORITHM_HPP

#include <htd/Globals.hpp>
#include <htd/IGraphDecompositionAlgorithm.hpp>

#include <htd/ILabelingFunction.hpp>
#include <htd/IMutableGraphDecomposition.hpp>
#include <htd/IGraphDecompositionManipulationOperation.hpp>

#include <vector>
#include <unordered_set>

namespace htd
{
    class BucketEliminationGraphDecompositionAlgorithm : public virtual htd::IGraphDecompositionAlgorithm
    {
        public:
            BucketEliminationGraphDecompositionAlgorithm(void);

            BucketEliminationGraphDecompositionAlgorithm(const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations);

            ~BucketEliminationGraphDecompositionAlgorithm();
            
            htd::IGraphDecomposition * computeDecomposition(const htd::IHypergraph & graph) const HTD_OVERRIDE;

            htd::IGraphDecomposition * computeDecomposition(const htd::IHypergraph & graph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations) const HTD_OVERRIDE;

            htd::IGraphDecomposition * computeDecomposition(const htd::IHypergraph & graph, int manipulationOperationCount, ...) const;

            void setManipulationOperations(const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations) HTD_OVERRIDE;

            BucketEliminationGraphDecompositionAlgorithm * clone(void) const HTD_OVERRIDE;

        protected:
            BucketEliminationGraphDecompositionAlgorithm & operator=(const BucketEliminationGraphDecompositionAlgorithm &) { return *this; }

        private:
            std::vector<htd::ILabelingFunction *> labelingFunctions_;

            std::vector<htd::IGraphDecompositionManipulationOperation *> postProcessingOperations_;

            htd::IMutableGraphDecomposition * computeMutableDecomposition(const htd::IHypergraph & graph) const;

            htd::vertex_t getMinimumVertex(const std::vector<htd::vertex_t> & vertices, const std::vector<htd::index_t> & vertexIndices) const;

            htd::vertex_t getMinimumVertex(const std::vector<htd::vertex_t> & vertices, const std::vector<htd::index_t> & vertexIndices, htd::vertex_t excludedVertex) const;
    };
}

#endif /* HTD_HTD_BUCKETELIMINATIONGRAPHDECOMPOSITIONALGORITHM_HPP */
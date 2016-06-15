/*
 * File:   TreeDecompositionOptimizationOperation.cpp
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

#ifndef HTD_HTD_TREEDECOMPOSITIONOPTIMIZATIONOPERATION_CPP
#define	HTD_HTD_TREEDECOMPOSITIONOPTIMIZATIONOPERATION_CPP

#include <htd/Globals.hpp>
#include <htd/Helpers.hpp>
#include <htd/TreeDecompositionOptimizationOperation.hpp>
#include <htd/ExhaustiveVertexSelectionStrategy.hpp>
#include <htd/CompressionOperation.hpp>

#include <algorithm>

htd::TreeDecompositionOptimizationOperation::TreeDecompositionOptimizationOperation(const htd::ITreeDecompositionFitnessFunction & fitnessFunction) : strategy_(new htd::ExhaustiveVertexSelectionStrategy()), fitnessFunction_(fitnessFunction.clone()), manipulationOperations_()
{

}

htd::TreeDecompositionOptimizationOperation::~TreeDecompositionOptimizationOperation()
{
    delete strategy_;

    delete fitnessFunction_;

    for (htd::ITreeDecompositionManipulationOperation * operation : manipulationOperations_)
    {
        delete operation;
    }
}

void htd::TreeDecompositionOptimizationOperation::apply(htd::IMutableTreeDecomposition & decomposition) const
{
    apply(decomposition, std::vector<htd::ILabelingFunction *>());
}

void htd::TreeDecompositionOptimizationOperation::apply(htd::IMutableTreeDecomposition & decomposition, const std::vector<htd::vertex_t> & relevantVertices) const
{
    HTD_UNUSED(relevantVertices)

    apply(decomposition, std::vector<htd::ILabelingFunction *>());
}

//TODO Remove
#include <htd/PreOrderTreeTraversal.hpp>

void debug(const htd::ITreeDecomposition & decomposition)
{
    htd::PreOrderTreeTraversal traversal;

    traversal.traverse(decomposition, [&](htd::vertex_t vertex, htd::vertex_t parent, std::size_t distanceToRoot)
    {
        HTD_UNUSED(parent)

        for (htd::index_t index = 0; index < distanceToRoot; ++index)
        {
            std::cout << "   ";
        }

        std::cout << "NODE " << vertex << ": " << decomposition.bagContent(vertex) << std::endl;
    });

    std::cout << std::endl;
}

void htd::TreeDecompositionOptimizationOperation::apply(htd::IMutableTreeDecomposition & decomposition, const std::vector<htd::ILabelingFunction *> & labelingFunctions) const
{
    if (decomposition.vertexCount() > 0)
    {
        const htd::ITreeDecompositionFitnessFunction & fitnessFunction = *fitnessFunction_;

        if (!manipulationOperations_.empty())
        {
            htd::CompressionOperation compressionOperation;

            compressionOperation.apply(decomposition);

            for (const htd::ITreeDecompositionManipulationOperation * operation : manipulationOperations_)
            {
                operation->apply(decomposition, labelingFunctions);
            }
        }

        htd::vertex_t initialRoot = decomposition.root();

        htd::vertex_t optimalRoot = initialRoot;

        double optimalFitness = fitnessFunction.fitness(decomposition);

        //TODO
        debug(decomposition);

        std::cout << "INITIAL FITNESS:     " << optimalFitness << "   (ROOT: " << optimalRoot << ")" << std::endl << std::endl << std::endl << std::endl;

        std::vector<htd::vertex_t> candidates;

        strategy_->selectVertices(decomposition, candidates);

        for (htd::vertex_t vertex : candidates)
        {
            if (vertex != initialRoot)
            {
                if (!manipulationOperations_.empty())
                {
                    htd::vertex_t currentVertex = vertex;

                    std::vector<htd::vertex_t> affectedVertices;

                    while (!decomposition.isRoot(currentVertex))
                    {
                        affectedVertices.push_back(currentVertex);

                        currentVertex = decomposition.parent(currentVertex);
                    }

                    affectedVertices.push_back(currentVertex);

                    std::cout << "AFFECTED VERTICES: " << affectedVertices << std::endl << std::endl;

                    decomposition.makeRoot(vertex);

                    htd::CompressionOperation compressionOperation;

                    compressionOperation.apply(decomposition, affectedVertices);

                    for (const htd::ITreeDecompositionManipulationOperation * operation : manipulationOperations_)
                    {
                        operation->apply(decomposition, affectedVertices, labelingFunctions);
                    }
                }
                else
                {
                    decomposition.makeRoot(vertex);
                }

                //TODO
                debug(decomposition);

                double currentFitness = fitnessFunction.fitness(decomposition);

                std::cout << "CURRENT FITNESS:     " << currentFitness << "   (ROOT: " << vertex << ")" << std::endl;

                if (currentFitness > optimalFitness)
                {
                    optimalFitness = currentFitness;

                    optimalRoot = vertex;

                    std::cout << "NEW OPTIMAL FITNESS: " << optimalFitness << "   (ROOT: " << optimalRoot << ")" << std::endl;
                }

                std::cout << std::endl << std::endl << std::endl;
            }
        }

        decomposition.makeRoot(optimalRoot);
    }
}

void htd::TreeDecompositionOptimizationOperation::apply(htd::IMutableTreeDecomposition & decomposition, const std::vector<htd::vertex_t> & relevantVertices, const std::vector<htd::ILabelingFunction *> & labelingFunctions) const
{
    HTD_UNUSED(relevantVertices)

    apply(decomposition, labelingFunctions);
}

void htd::TreeDecompositionOptimizationOperation::setManipulationOperations(const std::vector<htd::ITreeDecompositionManipulationOperation *> & manipulationOperations)
{
    manipulationOperations_.clear();

    std::copy(manipulationOperations.begin(), manipulationOperations.end(), std::back_inserter(manipulationOperations_));
}

void htd::TreeDecompositionOptimizationOperation::addManipulationOperation(htd::ITreeDecompositionManipulationOperation * manipulationOperation)
{
    manipulationOperations_.push_back(manipulationOperation);
}

void htd::TreeDecompositionOptimizationOperation::addManipulationOperations(const std::vector<htd::ITreeDecompositionManipulationOperation *> & manipulationOperations)
{
    std::copy(manipulationOperations.begin(), manipulationOperations.end(), std::back_inserter(manipulationOperations_));
}

void htd::TreeDecompositionOptimizationOperation::setVertexSelectionStrategy(htd::IVertexSelectionStrategy * strategy)
{
    HTD_ASSERT(strategy != nullptr)

    delete strategy_;

    strategy_ = strategy;
}

htd::TreeDecompositionOptimizationOperation * htd::TreeDecompositionOptimizationOperation::clone(void) const
{
    return new htd::TreeDecompositionOptimizationOperation(*fitnessFunction_);
}

#endif /* HTD_HTD_TREEDECOMPOSITIONOPTIMIZATIONOPERATION_CPP */

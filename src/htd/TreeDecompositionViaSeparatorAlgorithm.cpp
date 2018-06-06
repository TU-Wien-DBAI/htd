/*
* File:   TreeDecompositionViaSeparatorAlgorithm.cpp
*
* Author: MILAKOVIC Andrea
*
* Copyright 2015-2017, Michael Abseher
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

#ifndef HTD_HTD_TREEDECOMPOSITIONVIASEPARATORALGORITHM_CPP
#define HTD_HTD_TREEDECOMPOSITIONVIASEPARATORALGORITHM_CPP

#include <htd/Globals.hpp>
#include <htd/Helpers.hpp>

#include <htd/TreeDecompositionViaSeparatorAlgorithm.hpp>
#include <htd/ConnectedComponentAlgorithmFactory.hpp>
#include <htd/MultiHypergraph.hpp>
#include <htd/TreeDecompositionFactory.hpp>
#include <htd/GraphLabeling.hpp>
#include <htd/ILabelingFunction.hpp>
#include <htd/IMutableTreeDecomposition.hpp>
#include <htd/GraphPreprocessorFactory.hpp>
#include <htd/IGraphPreprocessor.hpp>
#include <htd/GraphSeparatorAlgorithmFactory.hpp>
#include <htd/ITreeDecompositionManipulationOperation.hpp>
#include <htd/TrivialTreeDecompositionAlgorithm.hpp>
#include <htd/WidthReductionOperation.hpp>

#include <algorithm>
#include <cstdarg>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

/**
*  Private implementation details of class htd::TreeDecompositionViaSeparatorAlgorithm.
*/
struct htd::TreeDecompositionViaSeparatorAlgorithm::Implementation
{
	/**
	*  Constructor for the implementation details structure.
	*
	*  @param[in] manager  The management instance to which the current object instance belongs.
	*/
	Implementation(const htd::LibraryInstance * const manager) : managementInstance_(manager), labelingFunctions_(), postProcessingOperations_()
	{

	}

	/**
	*  Copy constructor for the implementation details structure.
	*
	*  @param[in] original The original implementation details structure.
	*/
	Implementation(const Implementation & original) : managementInstance_(original.managementInstance_), labelingFunctions_(), postProcessingOperations_()
	{
		for (htd::ILabelingFunction * labelingFunction : original.labelingFunctions_)
		{
#ifndef HTD_USE_VISUAL_STUDIO_COMPATIBILITY_MODE
			labelingFunctions_.push_back(labelingFunction->clone());
#else
			labelingFunctions_.push_back(labelingFunction->cloneLabelingFunction());
#endif
		}

		for (htd::ITreeDecompositionManipulationOperation * postProcessingOperation : original.postProcessingOperations_)
		{
#ifndef HTD_USE_VISUAL_STUDIO_COMPATIBILITY_MODE
			postProcessingOperations_.push_back(postProcessingOperation->clone());
#else
			postProcessingOperations_.push_back(postProcessingOperation->cloneTreeDecompositionManipulationOperation());
#endif
		}
	}

	virtual ~Implementation()
	{		
		for (auto & labelingFunction : labelingFunctions_)
		{
			delete labelingFunction;
		}

		for (auto & postProcessingOperation : postProcessingOperations_)
		{
			delete postProcessingOperation;
		}
	}

	/**
	*  The management instance to which the current object instance belongs.
	*/
	const htd::LibraryInstance * managementInstance_;

	/**
	*  The labeling functions which are applied after the decomposition was computed.
	*/
	std::vector<htd::ILabelingFunction *> labelingFunctions_;

	/**
	*  The manipuation operations which are applied after the decomposition was computed.
	*/
	std::vector<htd::ITreeDecompositionManipulationOperation *> postProcessingOperations_;
		
};

htd::TreeDecompositionViaSeparatorAlgorithm::TreeDecompositionViaSeparatorAlgorithm(const htd::LibraryInstance * const manager) : implementation_(new Implementation(manager))
{

}

htd::TreeDecompositionViaSeparatorAlgorithm::TreeDecompositionViaSeparatorAlgorithm(const htd::LibraryInstance * const manager, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations) : implementation_(new Implementation(manager))
{
	setManipulationOperations(manipulationOperations);
}

htd::TreeDecompositionViaSeparatorAlgorithm::TreeDecompositionViaSeparatorAlgorithm(const htd::TreeDecompositionViaSeparatorAlgorithm & original) : implementation_(new Implementation(*(original.implementation_)))
{

}

htd::TreeDecompositionViaSeparatorAlgorithm::~TreeDecompositionViaSeparatorAlgorithm()
{

}

htd::ITreeDecomposition * htd::TreeDecompositionViaSeparatorAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph) const
{
	return computeDecomposition(graph, std::vector<htd::IDecompositionManipulationOperation *>());
}

htd::ITreeDecomposition * htd::TreeDecompositionViaSeparatorAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations) const
{
	htd::IGraphPreprocessor * preprocessor = implementation_->managementInstance_->graphPreprocessorFactory().createInstance();

	htd::IPreprocessedGraph * preprocessedGraph = preprocessor->prepare(graph);

	htd::ITreeDecomposition * ret = computeDecomposition(graph, *preprocessedGraph, manipulationOperations);

	delete preprocessedGraph;
	delete preprocessor;

	return ret;
}

htd::ITreeDecomposition * htd::TreeDecompositionViaSeparatorAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph) const
{
	return computeDecomposition(graph, preprocessedGraph, std::vector<htd::IDecompositionManipulationOperation *>());
}

htd::ITreeDecomposition * htd::TreeDecompositionViaSeparatorAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations) const
{ 
	//DECLARATION
	htd::IMutableTreeDecomposition * ret;
	htd::IGraphSeparatorAlgorithm * separatorAlgorithm;
	htd::IConnectedComponentAlgorithm * conAlg;
	MultiHypergraph * graphWithoutSeparator ;	
	htd::MultiHypergraph referrenceGraph(graph);	
	
	std::vector<std::vector<htd::vertex_t>> components;
	std::vector<std::vector<htd::vertex_t>> neighbors;

	//******* SEPARATOR and COMPONENTS *******//
	
	std::vector<vertex_t> * separator = separatorAlgorithm->computeSeparator(graph);	
	
	graphWithoutSeparator->assign(* graph.cloneMultiHypergraph());
	for (int i = 0; i < separator->size(); i++)
	{		
		graphWithoutSeparator->removeVertex(separator->at(i));
	}	

	conAlg->determineComponents(graph, components);

	//neighbors
	int counter = 0;
	for (htd::vertex_t v : *separator)
	{		
		for (std::vector<htd::vertex_t> comp : components)
		{
			for (htd::vertex_t n : graph.neighbors(v))
			{
				if (std::find(comp.begin(), comp.end(), n) != comp.end())
					comp.push_back(n);
			}			
		}
	}
	
	//******* TREE DECOMPOSITION *********//
	htd::vertex_t root = ret->insertRoot(* separator, graph.hyperedgesAtPositions(* separator));
	ret->addChild(root /* , vector<vertices>  */); 
	
	

	
	return ret;
}

htd::ITreeDecomposition * htd::TreeDecompositionViaSeparatorAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, int manipulationOperationCount, ...) const
{
	va_list arguments;

	va_start(arguments, manipulationOperationCount);

	std::vector<htd::IDecompositionManipulationOperation *> manipulationOperations;

	for (int manipulationOperationIndex = 0; manipulationOperationIndex < manipulationOperationCount; manipulationOperationIndex++)
	{
		manipulationOperations.push_back(va_arg(arguments, htd::IDecompositionManipulationOperation *));
	}

	va_end(arguments);

	return computeDecomposition(graph, manipulationOperations);
}

htd::ITreeDecomposition * htd::TreeDecompositionViaSeparatorAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, int manipulationOperationCount, ...) const
{
	va_list arguments;

	va_start(arguments, manipulationOperationCount);

	std::vector<htd::IDecompositionManipulationOperation *> manipulationOperations;

	for (int manipulationOperationIndex = 0; manipulationOperationIndex < manipulationOperationCount; manipulationOperationIndex++)
	{
		manipulationOperations.push_back(va_arg(arguments, htd::IDecompositionManipulationOperation *));
	}

	va_end(arguments);

	return computeDecomposition(graph, preprocessedGraph, manipulationOperations);
}

void htd::TreeDecompositionViaSeparatorAlgorithm::setManipulationOperations(const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations)
{
	for (auto & labelingFunction : implementation_->labelingFunctions_)
	{
		delete labelingFunction;
	}

	for (auto & postProcessingOperation : implementation_->postProcessingOperations_)
	{
		delete postProcessingOperation;
	}

	implementation_->labelingFunctions_.clear();

	implementation_->postProcessingOperations_.clear();

	addManipulationOperations(manipulationOperations);
}

void htd::TreeDecompositionViaSeparatorAlgorithm::addManipulationOperation(htd::IDecompositionManipulationOperation * manipulationOperation)
{
	bool assigned = false;

	htd::ILabelingFunction * labelingFunction = dynamic_cast<htd::ILabelingFunction *>(manipulationOperation);

	if (labelingFunction != nullptr)
	{
		implementation_->labelingFunctions_.emplace_back(labelingFunction);

		assigned = true;
	}

	htd::ITreeDecompositionManipulationOperation * newManipulationOperation = dynamic_cast<htd::ITreeDecompositionManipulationOperation *>(manipulationOperation);

	if (newManipulationOperation != nullptr)
	{
		implementation_->postProcessingOperations_.emplace_back(newManipulationOperation);

		assigned = true;
	}

	if (!assigned)
	{
		delete manipulationOperation;
	}
}

void htd::TreeDecompositionViaSeparatorAlgorithm::addManipulationOperations(const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations)
{
	for (htd::IDecompositionManipulationOperation * operation : manipulationOperations)
	{
		addManipulationOperation(operation);
	}
}

bool htd::TreeDecompositionViaSeparatorAlgorithm::isSafelyInterruptible(void) const
{
	return false;
}

const htd::LibraryInstance * htd::TreeDecompositionViaSeparatorAlgorithm::managementInstance(void) const HTD_NOEXCEPT
{
	return implementation_->managementInstance_;
}

void htd::TreeDecompositionViaSeparatorAlgorithm::setManagementInstance(const htd::LibraryInstance * const manager)
{
	HTD_ASSERT(manager != nullptr)

		implementation_->managementInstance_ = manager;
}



#endif /* HTD_HTD_SEPARATORBASEDTREEDECOMPOSITIONALGORITHM_CPP */

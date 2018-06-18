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
#include <htd/MinimumSeparatorAlgorithm.hpp>
#include <htd/WidthReductionOperation.hpp>
#include <htd/TrivialTreeDecompositionAlgorithm.hpp>

#include <algorithm>
#include <cstdarg>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>
#include <math.h> 

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
	Implementation(const htd::LibraryInstance * const manager) : managementInstance_(manager), separatorAlgorithm_(manager->graphSeparatorAlgorithmFactory().createInstance()), labelingFunctions_(), postProcessingOperations_()
	{

	}

	/**
	*  Copy constructor for the implementation details structure.
	*
	*  @param[in] original The original implementation details structure.
	*/
	Implementation(const Implementation & original) : managementInstance_(original.managementInstance_), separatorAlgorithm_(original.separatorAlgorithm_->clone()), labelingFunctions_(), postProcessingOperations_()
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
		delete separatorAlgorithm_;

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
	*  The underlying graph separator algorithm.
	*/
	htd::IGraphSeparatorAlgorithm * separatorAlgorithm_;

	/**
	*  The labeling functions which are applied after the decomposition was computed.
	*/
	std::vector<htd::ILabelingFunction *> labelingFunctions_;

	/**
	*  A boolean flag indicating whether the hyperedges induced by a respective bag shall be computed.
	*/
	bool computeInducedEdges_;

	/**
	*  Compute a new mutable tree decompostion of the given graph.
	*
	*  @param[in] graph                The graph which shall be decomposed.
	*  @param[in] preprocessedGraph    The input graph in preprocessed format.
	*
	*  @return A mutable tree decompostion of the given graph.
	*/
	void computeMutableDecomposition( htd::IMutableTreeDecomposition * ret, const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, vertex_t index, const int counter, const int max) const;


	void computeMutableDecompositionBagSize(htd::IMutableTreeDecomposition * ret, const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, vertex_t index, const int limit) const;

	void computeMutableDecompositionSeparators(htd::IMutableTreeDecomposition * ret, const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, vertex_t index, std::vector<vertex_t> oldSeparator) const;

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
	htd::IMutableTreeDecomposition * ret = dynamic_cast<htd::IMutableTreeDecomposition *>(implementation_ ->managementInstance_ ->treeDecompositionFactory().createInstance());
	
	int numberOfSteps = ceil(log(graph.vertexCount()));
	int countSteps = 0;
	
	

	//implementation_  -> computeMutableDecompositionBagSize(ret, graph, preprocessedGraph, NULL, 7);
	implementation_->computeMutableDecomposition(ret, graph, preprocessedGraph, NULL, countSteps, numberOfSteps);
	
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

void htd::TreeDecompositionViaSeparatorAlgorithm::setGraphSeparatorAlgorithm(htd::IGraphSeparatorAlgorithm * algorithm)
{
	HTD_ASSERT(algorithm != nullptr)

		delete implementation_->separatorAlgorithm_;

	implementation_->separatorAlgorithm_ = algorithm;
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

void htd::TreeDecompositionViaSeparatorAlgorithm::setComputeInducedEdgesEnabled(bool computeInducedEdgesEnabled)
{
	implementation_->computeInducedEdges_ = computeInducedEdgesEnabled;
}

bool htd::TreeDecompositionViaSeparatorAlgorithm::isComputeInducedEdgesEnabled(void) const
{
	return implementation_->computeInducedEdges_;
}

htd::TreeDecompositionViaSeparatorAlgorithm * htd::TreeDecompositionViaSeparatorAlgorithm::clone(void) const
{
	return new TreeDecompositionViaSeparatorAlgorithm(*this); // new htd::TreeDecompositionViaSeparatorAlgorithm(/*this*/);
}

void htd::TreeDecompositionViaSeparatorAlgorithm::Implementation::computeMutableDecomposition(htd::IMutableTreeDecomposition * ret, const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, vertex_t index, const int counter, const int max) const
{	
	MultiHypergraph  graphWithoutSeparator = MultiHypergraph(managementInstance_);
	std::vector<std::vector<htd::vertex_t>> components = std::vector<std::vector<htd::vertex_t>>();
	components.reserve(graph.vertexCount());  //TODO reserve - the correct number of places
	std::vector<std::vector<htd::vertex_t>> neighbors = std::vector<std::vector<htd::vertex_t>>();
	neighbors.reserve(graph.vertexCount());
	std::vector<htd::vertex_t> neighborsOneComponent = std::vector<htd::vertex_t>();
	neighborsOneComponent.reserve(graph.vertexCount());
	std::vector<htd::FilteredHyperedgeCollection> neighborsEdges = std::vector<htd::FilteredHyperedgeCollection>();
	neighborsEdges.reserve(graph.vertexCount());//TODO : initialize neighborsEdges

	std::vector<vertex_t> separator = *separatorAlgorithm_->computeSeparator(graph);  //separatorAlgorithm_->computeSeparator(graph);		

	if (index != NULL)
	{
		htd::vertex_t sep = ret->addChild(index, separator, graph.hyperedgesAtPositions(separator));
		index = sep;
	}
	else
	{
		htd::vertex_t root = ret->insertRoot(separator, graph.hyperedgesAtPositions(separator));
		index = root;
	}


	graphWithoutSeparator.assign(*graph.cloneMultiHypergraph());
	for (int i = 0; i < separator.size(); i++)
	{
		graphWithoutSeparator.removeVertex(separator.at(i));
	}
	managementInstance_->connectedComponentAlgorithmFactory().createInstance()->determineComponents(graphWithoutSeparator, components);

	for (std::vector<htd::vertex_t> comp : components)
	{
		neighborsOneComponent.clear();
		for (htd::vertex_t v : separator)
		{			
			for (htd::vertex_t n : graph.neighbors(v))
			{
				if (std::find(comp.begin(), comp.end(), n) != comp.end())
				{
					if (!(std::find(neighborsOneComponent.begin(), neighborsOneComponent.end(), n) != neighborsOneComponent.end()))
					{
						neighborsOneComponent.push_back(n);
					}					
					if (!(std::find(neighborsOneComponent.begin(), neighborsOneComponent.end(), v) != neighborsOneComponent.end()))
					{
						neighborsOneComponent.push_back(v);
					}				
				}
			}			
		}
		neighbors.push_back(neighborsOneComponent);
	}

	//******* TREE DECOMPOSITION *********//
	
	for (unsigned int i = 0; i < neighbors.size(); i++)
	{
		htd::vertex_t w = ret->addChild(index, neighbors.at(i), graph.hyperedgesAtPositions(neighbors.at(i)));
		if (counter + 1 <= max && (components.at(i).size()> neighbors.at(i).size() || (components.at(i).size() == neighbors.at(i).size() && components.at(i) != neighbors.at(i)) ) )  // TODO treba ustvari provjeriti da li su razliciti ili ne?
		{
			MultiHypergraph g = *graph.cloneMultiHypergraph();
			for (vertex_t v : graph.vertices())
			{
				if (!(std::find(components.at(i).begin(), components.at(i).end(), v) != components.at(i).end()))
					g.removeVertex(v);
			}
		
			computeMutableDecomposition(ret, g, preprocessedGraph, w, counter + 1, max);
			//if (components.at(i) != neighbors.at(i) && components.at(i).size() >= neighbors.at(i).size())  //TODO sortirati vektore
			//	ret->addChild(w, components.at(i), graph.hyperedgesAtPositions(components.at(i)));
		}
		else
		{
			if (components.at(i) != neighbors.at(i) && components.at(i).size() >= neighbors.at(i).size())  //TODO sortirati vektore
				ret->addChild(w, components.at(i), graph.hyperedgesAtPositions(components.at(i)));
		}
		//if (components.at(i) != neighbors.at(i) && components.at(i).size() >= neighbors.at(i).size())  //TODO sortirati vektore
		//		ret->addChild(w, components.at(i), graph.hyperedgesAtPositions(components.at(i)));
		
	}
}

void htd::TreeDecompositionViaSeparatorAlgorithm::Implementation::computeMutableDecompositionBagSize(htd::IMutableTreeDecomposition * ret, const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, vertex_t index, const int limit) const
{
	MultiHypergraph  graphWithoutSeparator = MultiHypergraph(managementInstance_);
	std::vector<std::vector<htd::vertex_t>> components = std::vector<std::vector<htd::vertex_t>>();
	components.reserve(graph.vertexCount());  //TODO reserve - the correct number of places
	std::vector<std::vector<htd::vertex_t>> neighbors = std::vector<std::vector<htd::vertex_t>>();
	neighbors.reserve(graph.vertexCount());
	std::vector<htd::vertex_t> neighborsOneComponent = std::vector<htd::vertex_t>();
	neighborsOneComponent.reserve(graph.vertexCount());
	std::vector<htd::FilteredHyperedgeCollection> neighborsEdges = std::vector<htd::FilteredHyperedgeCollection>();
	neighborsEdges.reserve(graph.vertexCount());//TODO : initialize neighborsEdges

	std::vector<vertex_t>  separator = *separatorAlgorithm_->computeSeparator(graph);  //separatorAlgorithm_->computeSeparator(graph);		

	if (index != NULL)
	{
		htd::vertex_t sep = ret->addChild(index, separator, graph.hyperedgesAtPositions(separator));
		index = sep;
	}
	else
	{
		htd::vertex_t root = ret->insertRoot(separator, graph.hyperedgesAtPositions(separator));
		index = root;
	}

	graphWithoutSeparator.assign(*graph.cloneMultiHypergraph());
	for (int i = 0; i < separator.size(); i++)
	{
		graphWithoutSeparator.removeVertex(separator.at(i));
	}
	managementInstance_->connectedComponentAlgorithmFactory().createInstance()->determineComponents(graphWithoutSeparator, components);  //determineComponents(graphWithoutSeparator, components); 

	for (std::vector<htd::vertex_t> comp : components)
	{
		neighborsOneComponent.clear();
		for (htd::vertex_t v : separator)
		{
			for (htd::vertex_t n : graph.neighbors(v))
			{
				if (std::find(comp.begin(), comp.end(), n) != comp.end())
				{
					if (!(std::find(neighborsOneComponent.begin(), neighborsOneComponent.end(), n) != neighborsOneComponent.end()))
					{
						neighborsOneComponent.push_back(n);
					}
					if (!(std::find(neighborsOneComponent.begin(), neighborsOneComponent.end(), v) != neighborsOneComponent.end()))
					{
						neighborsOneComponent.push_back(v);
					}
				}
			}
		}
		neighbors.push_back(neighborsOneComponent);
	}

	//******* TREE DECOMPOSITION *********//

	for (unsigned int i = 0; i < neighbors.size(); i++)
	{
		htd::vertex_t w = ret->addChild(index, neighbors.at(i), graph.hyperedgesAtPositions(neighbors.at(i)));
		if (components.at(i).size()> limit && components.at(i).size()> neighbors.at(i).size())
		{
			MultiHypergraph g = *graph.cloneMultiHypergraph();
			for (vertex_t v : graph.vertices())
			{
				if (!(std::find(components.at(i).begin(), components.at(i).end(), v) != components.at(i).end()))
					g.removeVertex(v);
			}
			computeMutableDecompositionBagSize(ret, g, preprocessedGraph, w, limit);

		}
		else
		{
			if (neighbors.at(i).size() <= components.at(i).size())
				ret->addChild(w, components.at(i), graph.hyperedgesAtPositions(components.at(i)));
		}

	}
}

void htd::TreeDecompositionViaSeparatorAlgorithm::Implementation::computeMutableDecompositionSeparators(htd::IMutableTreeDecomposition * ret, const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, vertex_t index, std::vector<vertex_t> oldSeparator) const
{
	MultiHypergraph  graphWithoutSeparator = MultiHypergraph(managementInstance_);
	std::vector<std::vector<htd::vertex_t>> components = std::vector<std::vector<htd::vertex_t>>();
	components.reserve(graph.vertexCount());  //TODO reserve - the correct number of places
	
	std::vector<vertex_t> separator = *separatorAlgorithm_->computeSeparator(graph);  //separatorAlgorithm_->computeSeparator(graph);		
	std::vector<vertex_t> separators = separator;

	if (index != NULL)
	{
		htd::vertex_t sep = ret->addChild(index, separator, graph.hyperedgesAtPositions(separator));
		index = sep;
	}
	else
	{
		for (vertex_t v : oldSeparator)
		{
			if (!(std::find(separator.begin(), separator.end(), v) != separator.end()))
				separators.push_back(v);
		}
		htd::vertex_t root = ret->insertRoot(separators, graph.hyperedgesAtPositions(separators));
		index = root;
	}

	if (separator.size() > 1) 
	{
		graphWithoutSeparator.assign(*graph.cloneMultiHypergraph());
		for (int i = 0; i < separator.size(); i++)
		{
			graphWithoutSeparator.removeVertex(separator.at(i));
		}
		managementInstance_->connectedComponentAlgorithmFactory().createInstance()->determineComponents(graphWithoutSeparator, components);

		for (std::vector<htd::vertex_t> comp : components)
		{
			MultiHypergraph g = *graph.cloneMultiHypergraph();
			for (vertex_t v : graph.vertices())
			{
				if (!(std::find(comp.begin(), comp.end(), v) != comp.end()))
					g.removeVertex(v);
			}

			computeMutableDecompositionSeparators(ret, g, preprocessedGraph, index, separator);
		}

	}
	else
	{
		MultiHypergraph g = *graph.cloneMultiHypergraph();
		ret->addChild(index, g.vertexVector() , g.hyperedgesAtPositions(g.vertexVector()));
	}	
}

#endif /* HTD_HTD_SEPARATORBASEDTREEDECOMPOSITIONALGORITHM_CPP */

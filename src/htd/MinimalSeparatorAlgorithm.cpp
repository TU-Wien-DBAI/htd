#ifndef HTD_HTD_MINIMALSEPARATORALGORITHM_CPP
#define HTD_HTD_MINIMALSEPARATORALGORITHM_CPP

#include <htd/MinimalSeparatorAlgorithm.hpp>

#include <htd/FlowNetworkStructure.hpp>
#include <htd/DinitzMaxFlowAlgorithm.hpp>

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include "..\..\include\htd\MinimalSeparatorAlgorithm.hpp"

/**
*  Private implementation details of class htd::MinimalSeparatorAlgorithm.
*/
struct htd::MinimalSeparatorAlgorithm::Implementation
{
	/**
	*  Constructor for the implementation details structure.
	*
	*  @param[in] manager   The management instance to which the current object instance belongs.
	*/
	Implementation(const htd::LibraryInstance * const manager) : managementInstance_(manager) 
	{

	}

	virtual ~Implementation()
	{

	}

	/**
	*  The management instance to which the current object instance belongs.
	*/
	const htd::LibraryInstance * managementInstance_;
};

htd::MinimalSeparatorAlgorithm::MinimalSeparatorAlgorithm(const htd::LibraryInstance * const manager) : implementation_(new Implementation(manager))
{

}

htd::MinimalSeparatorAlgorithm::~MinimalSeparatorAlgorithm(void)
{

}

std::vector<htd::vertex_t> getKeys(std::unordered_map<htd::vertex_t, bool> list)
{
	std::vector<htd::vertex_t> keys;
	
	keys.reserve(list.size());
	
	for (auto kv : list)
	{
		keys.push_back(kv.first);
	}

	return keys;
}

void removeVertex(const htd::IGraphStructure & graph, htd::vertex_t vertex)
{
	// TODO
	htd::ConstCollection<htd::vertex_t> selectedNeighborhood = graph.neighbors(vertex);

	for (htd::vertex_t neighbor : selectedNeighborhood)
	{
		htd::ConstCollection<htd::vertex_t> & currentNeighborhood = graph.neighbors(neighbor);

	}


}

std::vector<htd::vertex_t> * htd::MinimalSeparatorAlgorithm::computeSeparator(const htd::IGraphStructure & graph) const
{
	std::size_t k = 5; // Implement k

	std::vector<htd::vertex_t> * separator = new std::vector<htd::vertex_t>();

	std::size_t counter = 0;

	int graphSize = graph.vertexCount();
	
	std::unordered_map<htd::vertex_t, bool>* list = new std::unordered_map<htd::vertex_t, bool>[graph.vertexCount];
	
	for (htd::vertex_t vertex : graph.vertices())
	{
		std::size_t degree = graph.neighborCount(vertex);
	
		list[degree].emplace(vertex, true);
	}

	htd::index_t current = graph.vertexCount - 1;

	while (counter < k)
	{
		while (getKeys(list[current]).size == 0)
		{
			current--;
		}
		
		std::vector<htd::vertex_t> keys = getKeys(list[current]);

		while (keys.size > 0)
		{
			htd::vertex_t choice = keys.at(0);

			separator->push_back(choice);

			counter++;

			htd::ConstCollection<htd::vertex_t> neighbors = graph.neighbors(choice);
			
			for (auto n : neighbors)
			{
				std:size_t neighborDegree = graph.neighborCount(n);
				list[neighborDegree].erase(n);
				list[neighborDegree - 1].emplace(n, true);
			}

			removeVertex(graph, choice);
			
			list[current].erase(choice);

			if (separator->size == k) 
			{
				break;
			}
			keys = getKeys(list[current]);
		}
	}

	return separator;
}


const htd::LibraryInstance * htd::MinimalSeparatorAlgorithm::managementInstance(void) const HTD_NOEXCEPT
{
	return implementation_->managementInstance_;
}

void htd::MinimalSeparatorAlgorithm::setManagementInstance(const htd::LibraryInstance * const manager)
{
	HTD_ASSERT(manager != nullptr)

		implementation_->managementInstance_ = manager;
}

htd::MinimalSeparatorAlgorithm * htd::MinimalSeparatorAlgorithm::clone(void) const
{
	return new htd::MinimalSeparatorAlgorithm(managementInstance());
}

#endif /* HTD_HTD_MINIMALSEPARATORALGORITHM_CPP */

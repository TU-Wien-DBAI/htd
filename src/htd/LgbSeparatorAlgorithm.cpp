#ifndef HTD_HTD_LGBSEPARATORALGORITHM_CPP
#define HTD_HTD_LGBSEPARATORALGORITHM_CPP

#include <htd/LgbSeparatorAlgorithm.hpp>
#include <htd/MultiGraph.hpp>
#include <htd/LabeledMultiGraph.hpp>
#include <htd/DepthFirstConnectedComponentAlgorithm.hpp>

#include <algorithm>
#include <unordered_map>
#include <iostream>

/**
*  Private implementation details of class htd::LgbSeparatorAlgorithm.
*/
struct htd::LgbSeparatorAlgorithm::Implementation
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

htd::LgbSeparatorAlgorithm::LgbSeparatorAlgorithm(const htd::LibraryInstance * const manager) : implementation_(new Implementation(manager))
{

}

htd::LgbSeparatorAlgorithm::~LgbSeparatorAlgorithm(void)
{

}

std::vector<htd::vertex_t> * htd::LgbSeparatorAlgorithm::computeSeparator(const htd::IGraphStructure & graph) const
{
	htd::MultiGraph * newGraph = new htd::MultiGraph(managementInstance(), graph.vertexCount());

	for (auto v1 : graph.vertices())
	{
		for (auto v2 : graph.neighbors(v1))
		{
			if (!newGraph->isEdge(v1, v2))
			{
				newGraph->addEdge(v1, v2);
			}
		}
	}

	return computeSeparator(newGraph);
}

std::vector<htd::vertex_t> * htd::LgbSeparatorAlgorithm::computeSeparator(htd::MultiGraph * graph) const
{
	std::vector<htd::vertex_t> * separator = new std::vector<htd::vertex_t>();

	htd:


	return separator;
}

const htd::LibraryInstance * htd::LgbSeparatorAlgorithm::managementInstance(void) const HTD_NOEXCEPT
{
	return implementation_->managementInstance_;
}

void htd::LgbSeparatorAlgorithm::setManagementInstance(const htd::LibraryInstance * const manager)
{
	HTD_ASSERT(manager != nullptr)

		implementation_->managementInstance_ = manager;
}

htd::LgbSeparatorAlgorithm * htd::LgbSeparatorAlgorithm::clone(void) const
{
	return new htd::LgbSeparatorAlgorithm(managementInstance());
}

#endif /* HTD_HTD_HIGHESTDEGREESEPARATORALGORITHM_CPP */

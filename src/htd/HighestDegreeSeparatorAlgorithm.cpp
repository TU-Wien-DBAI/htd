#ifndef HTD_HTD_HIGHESTDEGREESEPARATORALGORITHM_CPP
#define HTD_HTD_HIGHESTDEGREESEPARATORALGORITHM_CPP

#include <htd/HighestDegreeSeparatorAlgorithm.hpp>
#include <htd/MultiGraph.hpp>
#include <htd/DepthFirstConnectedComponentAlgorithm.hpp>

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <chrono>
#include <ctime>
#include <random>
/**
*  Private implementation details of class htd::HighestDegreeAlgorithm.
*/
struct htd::HighestDegreeSeparatorAlgorithm::Implementation
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

htd::HighestDegreeSeparatorAlgorithm::HighestDegreeSeparatorAlgorithm(const htd::LibraryInstance * const manager) : implementation_(new Implementation(manager))
{

}

htd::HighestDegreeSeparatorAlgorithm::~HighestDegreeSeparatorAlgorithm(void)
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

std::vector<htd::vertex_t> * htd::HighestDegreeSeparatorAlgorithm::computeSeparator(const htd::IGraphStructure & graph) const
{
	htd::MultiGraph * newGraph = new htd::MultiGraph(managementInstance(), graph.vertexCount());
	
	const std::clock_t begin_time = std::clock();
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
	std::cout << "Algorithm execution time: " << float(clock() - begin_time) / CLOCKS_PER_SEC << " second(s)." << std::endl;

	return computeSeparator(newGraph);
}

std::vector<htd::vertex_t> * htd::HighestDegreeSeparatorAlgorithm::computeSeparator(htd::MultiGraph * graph) const
{
	const std::clock_t begin_time = std::clock();
	
	htd::DepthFirstConnectedComponentAlgorithm conectedComponentsAlgorithm(managementInstance());

	std::vector<std::vector<htd::vertex_t>> target = std::vector<std::vector<htd::vertex_t>>();

	std::vector<htd::vertex_t> * separator = new std::vector<htd::vertex_t>();

	std::unordered_map<htd::vertex_t, bool>* list = new std::unordered_map<htd::vertex_t, bool>[graph->vertexCount()];

	conectedComponentsAlgorithm.determineComponents(*graph, target);

	bool connected = target.size() < 2;

	if (!connected) 
	{
		std::cout << "The input graph is not connected" << std::endl;
		return separator;
	}

	target.clear();

	std::size_t counter = 0;

	for (htd::vertex_t vertex : graph->vertices())
	{
		std::size_t degree = graph->neighborCount(vertex);

		list[degree].emplace(vertex, true);
	}

	htd::index_t current = graph->vertexCount() - 1;

	while (connected)
	{
		while (getKeys(list[current]).size() == 0)
		{
			current--;
		}

		std::vector<htd::vertex_t> keys = getKeys(list[current]);

		while (keys.size() > 0)
		{

			std::random_device random_device;
			std::mt19937 engine{ random_device() };
			std::uniform_int_distribution<int> dist(0, keys.size() - 1);
			htd::vertex_t choice = keys[dist(engine)];

			separator->push_back(choice);

			counter++;

			htd::ConstCollection<htd::vertex_t> neighbors = graph->neighbors(choice);

			for (auto n : neighbors)
			{
				std::size_t neighborDegree = graph->neighborCount(n);

				list[neighborDegree].erase(n);

				list[neighborDegree - 1].emplace(n, true);
			}

			graph->removeVertex(choice);

			list[current].erase(choice);

			conectedComponentsAlgorithm.determineComponents(*graph, target);

			connected = target.size() < 2;

			target.clear();

			if (!connected)
			{
				break;
			}
			keys = getKeys(list[current]);
		}
	}
	std::cout << "Algorithm execution time: " << float(clock() - begin_time) / CLOCKS_PER_SEC << " second(s)." << std::endl;

	conectedComponentsAlgorithm.determineComponents(*graph, target);

	std::cout << "The input graph is divided in: " << target.size() << " components" << std::endl;

	return separator;
}

const htd::LibraryInstance * htd::HighestDegreeSeparatorAlgorithm::managementInstance(void) const HTD_NOEXCEPT
{
	return implementation_->managementInstance_;
}

void htd::HighestDegreeSeparatorAlgorithm::setManagementInstance(const htd::LibraryInstance * const manager)
{
	HTD_ASSERT(manager != nullptr)

		implementation_->managementInstance_ = manager;
}

htd::HighestDegreeSeparatorAlgorithm * htd::HighestDegreeSeparatorAlgorithm::clone(void) const
{
	return new htd::HighestDegreeSeparatorAlgorithm(managementInstance());
}

#endif /* HTD_HTD_HIGHESTDEGREESEPARATORALGORITHM_CPP */

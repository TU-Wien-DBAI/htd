#ifndef HTD_HTD_LGBSEPARATORALGORITHM_CPP
#define HTD_HTD_LGBSEPARATORALGORITHM_CPP

#include <htd/LgbSeparatorAlgorithm.hpp>
#include <htd/MultiGraph.hpp>
#include <htd/LabeledMultiGraph.hpp>
#include <htd/DepthFirstConnectedComponentAlgorithm.hpp>
#include <algorithm>
#include <unordered_map>
#include <iostream>
#include <random>
#include <numeric>
#include <cmath> 

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
	htd::LabeledMultiGraph * newGraph = new htd::LabeledMultiGraph(managementInstance(), graph.vertexCount());

	int label = 1;

	for (auto v1 : graph.vertices())
	{
		for (auto v2 : graph.neighbors(v1))
		{
			if (!newGraph->isEdge(v1, v2) && !newGraph->isEdge(v2, v1))
			{

				htd::id_t edgeId = newGraph->addEdge(v1, v2);

				newGraph->setEdgeLabel("" + label, edgeId, new htd::Label<int>(label));

				label++;
 			}
		}
	}

	return computeSeparator(newGraph);
}

struct pair {
	struct bin *nij;
	const htd::ILabel *lij;
};

struct bin {
	htd::vertex_t node;
	int degree;
	int partition;

	int gain;
	bool locked;
	
	const htd::ILabel *l1n;
	const htd::ILabel *l2n;

	int l1s;
	int l2s;

	int l1o;
	int l2o;

	std::vector<struct pair> *pairs;

	struct bin *next, *prev;
};

struct binGainPair
{
	int gain;
	struct bin *bin;
};

struct move
{
	int seq;
	int gain;
	int loadDiff;
	int cutSize;
	struct bin * bin;
};

void initBins(std::vector<struct bin*> * binspack, std::vector<binGainPair *> * binGainPairsPartition1, std::vector<binGainPair *> * binGainPairsPartition2)
{
	int i = -2;
	while (i <= 2) 
	{
		struct binGainPair *bg1 = new struct binGainPair;
		struct binGainPair *bg2 = new struct binGainPair;
		bg1->bin = nullptr;
		bg1->gain = i;
		bg2->bin = nullptr;
		bg2->gain = i;
		binGainPairsPartition1->push_back(bg1);
		binGainPairsPartition2->push_back(bg2);
		i++;
	}
	struct bin * prev = new struct bin;
	
	prev->gain = 0;
	prev->node = -10;
	binspack->at(0)->prev = prev;
	
	bool isPartition1filled = false;
	bool isPartition2filled = false;

	for (int i = 0; i < binspack->size(); i++)
	{
		if (i != binspack->size() - 1)
		{
			binspack->at(i)->next = binspack->at(i + 1);
		}

		if (i != 0)
		{
			binspack->at(i)->prev = binspack->at(i - 1);
		}

		if (binspack->at(i)->partition == 1 && !isPartition1filled)
		{
			binGainPairsPartition1->at(2)->bin = binspack->at(i);
			isPartition1filled = true;
		}
		
		if (binspack->at(i)->partition == 2 && !isPartition2filled)
		{
			binGainPairsPartition2->at(2)->bin = binspack->at(i);
			isPartition2filled = true;
		}
	}
	
	return;
}

void printBins(std::vector<struct bin*> * binspack)
{
	if (binspack->size() == 0)
	{
		printf("Empty List!! \n");
		return;
	}

	for (int i = 0; i < binspack->size(); i++)
	{
		struct bin * temp = binspack->at(i);
		printf("Node: %d", temp->node);
		printf(" Degree: %d", temp->degree);
		printf(" Partition: %d", temp->partition);

		printf(" ;Gain: %d", temp->gain);
		printf(" Locked: %s", temp->locked ? "yes" : "no");

		std::ostream stream(nullptr);
		stream.rdbuf(std::cout.rdbuf());
		printf(" ;L1n ");
		temp->l1n->print(stream);
		printf(" ;L1s %d", temp->l1s);
		printf(" ;L1o %d", temp->l1o);
		
		printf(" L2n ");
		temp->l2n->print(stream);
		printf(" ;L2s %d", temp->l2s);
		printf(" ;L2o %d", temp->l2o);

		printf("\n");
	}
	printf("\n");
}

void fillLineGraph(htd::LabeledMultiGraph * lineGraph, htd::LabeledMultiGraph * graph) {

	for (auto v : graph->vertices())
	{
		graph->setVertexLabel("edge", v, new htd::Label<int>((int)v));
		if (graph->neighborCount(v) == 1)
		{

			htd::vertex_t artificial = lineGraph->addVertex();

			lineGraph->setVertexLabel("Artificial" + artificial, artificial, new htd::Label<std::string>("Artificial"));

			htd::id_t edgeId = graph->associatedEdgeIds(v, graph->neighborAtPosition(v, 0)).empty() ? *graph->associatedEdgeIds(graph->neighborAtPosition(v, 0), v).begin() : *graph->associatedEdgeIds(v, graph->neighborAtPosition(v, 0)).begin();

			htd::id_t newEdgeId = lineGraph->addEdge(edgeId, artificial);

			lineGraph->setEdgeLabel("edge", newEdgeId, new htd::Label<int>((int)v));
		}
		else
		{
			for (auto neigh1 : graph->neighbors(v))
			{
				htd::id_t edgeId1 = graph->associatedEdgeIds(v, neigh1).empty() ? *graph->associatedEdgeIds(neigh1, v).begin() : *graph->associatedEdgeIds(v, neigh1).begin();
			
				for (auto neigh2 : graph->neighbors(v))
				{
					if (neigh1 != neigh2)
					{
						htd::id_t edgeId2 = graph->associatedEdgeIds(v, neigh2).empty() ? *graph->associatedEdgeIds(neigh2, v).begin() : *graph->associatedEdgeIds(v, neigh2).begin();
						
						htd::id_t newEdgeId = edgeId1 < edgeId2 ? lineGraph->addEdge(edgeId1, edgeId2) : lineGraph->addEdge(edgeId2, edgeId1);
				
						lineGraph->setEdgeLabel("edge", newEdgeId, new htd::Label<int>((int)v));
					}
				}
			}
		}
	}

	return;
}

void initializePartitions(htd::LabeledMultiGraph * lineGraph, std::size_t numberOfnodes, std::vector<htd::vertex_t> *partition1, std::vector<htd::vertex_t> *partition2, std::vector<struct bin*> * binspack)
{
	std::vector<int> seq(numberOfnodes);
	std::iota(seq.begin(), seq.end(), 1);

	std::random_device rand_dev;
	std::mt19937 twister(rand_dev());
	std::shuffle(seq.begin(), seq.end(), twister);

	for (int i = 1; i <= std::floor(numberOfnodes / 2); i++)
	{
		htd::vertex_t vertex = seq.at(i-1);
		partition1->push_back(vertex);
		
		struct bin *bin = new struct bin;
		bin->node = vertex;
		bin->degree = lineGraph->neighborCount(vertex);
		bin->gain = 0;
		bin->locked = false;
		bin->partition = 1;
		binspack->push_back(bin);
	}
	for (int i = static_cast<int>(std::floor(numberOfnodes / 2)) + 1; i <= numberOfnodes; i++)
	{
		htd::vertex_t vertex = seq.at(i-1);
		partition2->push_back(vertex);
		
		struct bin *bin = new struct bin;
		bin->node = vertex;
		bin->degree = lineGraph->neighborCount(vertex);
		bin->gain = 0;
		bin->locked = false;
		bin->partition = 2;
		binspack->push_back(bin);
	}

	return;
}

void computeLs(htd::LabeledMultiGraph * lineGraph, std::size_t numberOfnodes, std::vector<htd::vertex_t> *partition1, std::vector<htd::vertex_t> *partition2, std::vector<struct bin*> * binspack)
{
	for (auto v : lineGraph->vertices())
	{
		if (v <= numberOfnodes)
		{
			int l1s = 0, l1o = 0, l2s = 0, l2o = 0;
			const htd::ILabel *l1n = nullptr, *l2n = nullptr;
			std::vector<struct pair> *pairs = new std::vector<struct pair>();

			for (auto n : lineGraph->neighbors(v))
			{
				htd::id_t edgeId = v < n ? *lineGraph->associatedEdgeIds(v, n).begin() : *lineGraph->associatedEdgeIds(n, v).begin();
				if (l1n != nullptr && !(*l1n == lineGraph->edgeLabel("edge", edgeId)))
				{
					l2n = &lineGraph->edgeLabel("edge", edgeId);
					break;
				}
				l1n = &lineGraph->edgeLabel("edge", edgeId);
			}

			for (auto n : lineGraph->neighbors(v))
			{
				if (n <= numberOfnodes)
				{
					struct pair * pair = new struct pair;
					pair->nij = binspack->at(n - 1);
					htd::id_t edgeId = v < n ? *lineGraph->associatedEdgeIds(v, n).begin() : *lineGraph->associatedEdgeIds(n, v).begin();
					pair->lij = &lineGraph->edgeLabel("edge", edgeId);
					pairs->push_back(*pair);
				}
			}

			for (auto n : lineGraph->neighbors(v))
			{
				htd::id_t edgeId = v < n ? *lineGraph->associatedEdgeIds(v, n).begin() : *lineGraph->associatedEdgeIds(n, v).begin();

				bool inTheSamePartition = binspack->at(v - 1)->partition == 1 ? (std::find(partition1->begin(), partition1->end(), n) != partition1->end()) : (std::find(partition2->begin(), partition2->end(), n) != partition2->end());

				if (lineGraph->edgeLabel("edge", edgeId) == *l1n && (inTheSamePartition || n > numberOfnodes))
				{
					l1s++;
				}
				else if (lineGraph->edgeLabel("edge", edgeId) == *l1n && (!inTheSamePartition || n > numberOfnodes))
				{
					l1o++;
				}
				else if (lineGraph->edgeLabel("edge", edgeId) == *l2n && (inTheSamePartition || n > numberOfnodes))
				{
					l2s++;
				}
				else if (lineGraph->edgeLabel("edge", edgeId) == *l2n && (!inTheSamePartition || n > numberOfnodes))
				{
					l2o++;
				}
			}

			binspack->at(v - 1)->l1n = l1n, binspack->at(v - 1)->l1s = l1s, binspack->at(v - 1)->l1o = l1o, binspack->at(v - 1)->l2n = l2n, binspack->at(v - 1)->l2s = l2s, binspack->at(v - 1)->l2o = l2o;
			binspack->at(v - 1)->pairs = pairs;
		}
	}
}

void calculateGains(htd::LabeledMultiGraph * lineGraph, const size_t &numberOfnodes, std::vector<bin *> * binspack);

void updateBinPairPartitions(std::vector<binGainPair *> * binGainPairsPartition1, std::vector<binGainPair *> * binGainPairsPartition2, std::vector<bin *> * binspack)
{
	for (int i = 0; i < 5; i++)
	{
		binGainPairsPartition1->at(i)->bin = nullptr;
		binGainPairsPartition2->at(i)->bin = nullptr;
	}

	for (int i = 0; i < binspack->size(); i++)
	{
		if (binspack->at(i)->prev == nullptr)
		{
			struct bin *temp = binspack->at(i);

			while (temp != nullptr)
			{
				if (temp->partition == 1 && !temp->locked && binGainPairsPartition1->at(temp->gain + 2)->bin == nullptr)
				{
					binGainPairsPartition1->at(temp->gain + 2)->bin = temp;
				}
				else if (temp->partition == 2 && !temp->locked && binGainPairsPartition2->at(temp->gain + 2)->bin == nullptr)
				{
					binGainPairsPartition2->at(temp->gain + 2)->bin = temp;
				}
				temp = temp->next;
			}
		}
	}
}

void updateBins(std::vector<bin *> * binspack)
{
	for (int i = 0; i < binspack->size(); i++)
	{
		struct bin *bin = binspack->at(i);
		if ((i != binspack->size() - 1) && bin->next != nullptr && bin->next->gain != bin->gain)
		{
			bin->next = nullptr;
		}
		if (bin->prev != nullptr && bin->prev->gain != bin->gain)
		{
			bin->prev = nullptr;
		}

		if (bin->next == nullptr)
		{
			for (int j = i + 1; j < binspack->size(); j++)
			{
				if (binspack->at(j)->gain == bin->gain)
				{
					bin->next = binspack->at(j);
					binspack->at(j)->prev = bin;
					binspack->at(j)->next = nullptr;
					break;
				}
			}
		}
	}
}

void calculateGainsAndUpdateBins(htd::LabeledMultiGraph * lineGraph, std::size_t numberOfnodes, std::vector<struct bin*> * binspack, std::vector<binGainPair *> * binGainPairsPartition1, std::vector<binGainPair *> * binGainPairsPartition2)
{
	calculateGains(lineGraph, numberOfnodes, binspack);

	updateBins(binspack);

	updateBinPairPartitions(binGainPairsPartition1, binGainPairsPartition2, binspack);
}

void adjustLabels(htd::LabeledMultiGraph * lineGraph, std::vector<struct bin*> * binspack, struct bin * maxGain, std::vector<htd::vertex_t> *partition1, std::vector<htd::vertex_t> *partition2)
{
	int numberOfnodes = binspack->size();

	for (auto v : lineGraph->neighbors(maxGain->node))
	{
		if (v <= numberOfnodes)
		{
			int l1s = 0, l1o = 0, l2s = 0, l2o = 0;
			const htd::ILabel *l1n = binspack->at(v - 1)->l1n;
			const htd::ILabel *l2n = binspack->at(v - 1)->l2n;

			for (auto n : lineGraph->neighbors(v))
			{
				htd::id_t edgeId = v < n ? *lineGraph->associatedEdgeIds(v, n).begin() : *lineGraph->associatedEdgeIds(n, v).begin();

				bool inTheSamePartition = binspack->at(v - 1)->partition == 1 ? (std::find(partition1->begin(), partition1->end(), n) != partition1->end()) : (std::find(partition2->begin(), partition2->end(), n) != partition2->end());

				if (lineGraph->edgeLabel("edge", edgeId) == *l1n && (inTheSamePartition || n > numberOfnodes))
				{
					l1s++;
				}
				else if (lineGraph->edgeLabel("edge", edgeId) == *l1n && (!inTheSamePartition || n > numberOfnodes))
				{
					l1o++;
				}
				else if (lineGraph->edgeLabel("edge", edgeId) == *l2n && (inTheSamePartition || n > numberOfnodes))
				{
					l2s++;
				}
				else if (lineGraph->edgeLabel("edge", edgeId) == *l2n && (!inTheSamePartition || n > numberOfnodes))
				{
					l2o++;
				}
			}

			binspack->at(v - 1)->l1s = l1s, binspack->at(v - 1)->l1o = l1o, binspack->at(v - 1)->l2s = l2s, binspack->at(v - 1)->l2o = l2o;
		}
	}
}

void calculateGains(htd::LabeledMultiGraph * lineGraph, const size_t &numberOfnodes, std::vector<bin *> * binspack)
{
	for (auto v : lineGraph->vertices())
	{
		if (v <= numberOfnodes)
		{
			struct bin *bin = binspack->at(v - 1);
			if (bin->l1s == 0 && bin->l1o > 0)
			{
				bin->gain = bin->gain + 1; //weight(lineGraph, bin->l1n);
			}
			else if (bin->l1s > 0 && bin->l1o == 0)
			{
				bin->gain = bin->gain - 1; //weight(lineGraph, bin->l1n);
			}
			else if (bin->l1s > 0 && bin->l1o > 0)
			{
				// do nothing
			}

			if (bin->l2s == 0 && bin->l2o > 0)
			{
				bin->gain = bin->gain + 1;//weight(lineGraph, bin->l2n);
			}
			else if (bin->l2s > 0 && bin->l2o == 0)
			{
				bin->gain = bin->gain - 1; //weight(lineGraph, bin->l2n);
			}
			else if (bin->l2s > 0 && bin->l2o > 0)
			{
				// do nothing
			}
		}
	}
}

bool contains(std::vector<const htd::ILabel*> *labels, const htd::ILabel* label)
{
	bool contained = false;

	for (int i = 0; i < labels->size(); i++)
	{
		if (*label == *labels->at(i))
		{
			contained = true;
			break;
		}
	}

	return contained;
}

void findDistinctLabelsInPartitions(std::vector<const htd::ILabel*> *cutLabels, std::vector<const htd::ILabel*> *part1Labels, std::vector<const htd::ILabel*> *part2Labels, std::vector<struct bin*> * binspack)
{
	for (int i = 0; i < binspack->size(); i++)
	{
		struct bin * temp = binspack->at(i);

		if (temp->l1o > 0) 
		{
			if (!contains(cutLabels, temp->l1n))
			{
				cutLabels->push_back(temp->l1n); 
			}
		}
		else if (temp->l1s > 0)
		{
			if (temp->partition == 1) 
			{ 
				if (!contains(part1Labels, temp->l1n))
				{
					part1Labels->push_back(temp->l1n); 
				}
			} 
			else 
			{ 
				if (!contains(part2Labels, temp->l1n))
				{
					part2Labels->push_back(temp->l1n); 
				}
			}
		}

		if (temp->l2o > 0)
		{
			if (!contains(cutLabels, temp->l2n))
			{
				cutLabels->push_back(temp->l2n); 
			}
		}
		else if (temp->l2s > 0) 
		{
			if (temp->partition == 1) 
			{ 
				if (!contains(part1Labels, temp->l2n))
				{
					part1Labels->push_back(temp->l2n); 			
				}
			}
			else 
			{ 
				if (!contains(part2Labels, temp->l2n))
				{
					part2Labels->push_back(temp->l2n); 
				}
			}
		}

		temp = temp->next;
	}
	return;
 }

bool sortBins(struct bin *bin1, struct bin *bin2)
{ 
	if (bin1 != nullptr && bin2 == nullptr) { return true; }
	else if (bin1 == nullptr && bin2 != nullptr) { return false; }
	else if (bin1 != nullptr && bin2 != nullptr) { return bin1->node < bin2->node; }
	else { return false; }
}

void updateW1W2(int *wSrc, int *wDst, struct bin * binToBeMoved)
{
	if (binToBeMoved->l1o > 0)
	{
		if (binToBeMoved->l1s > 0) { /* nothing changes */ }
		else { *wDst = *wDst + 1; /* weight(L1n) */ }
	}
	else if (binToBeMoved->l1o = 0) { *wSrc = *wSrc - 1; }

	if (binToBeMoved->l2o > 0)
	{
		if (binToBeMoved->l2s > 0) { /* nothing changes */ }
		else { *wDst = *wDst + 1; }
	}
	else if (binToBeMoved->l2o = 0) { *wSrc = *wSrc - 1; }

	return;
}

int findK(std::vector<struct move> * moves)
{
	int max = 0, gain = 0, k = 0;

	for (int i = 0; i < moves->size(); i++)
	{
		for (int j = 0; j < i; j++) {
			gain += moves->at(j).gain;
		}

		if (gain > max) { max = gain; k = i; }
	}

	return k;
}

void updateTheGainsOfAdjacentNodes(htd::LabeledMultiGraph * lineGraph, int w1, int w2, std::vector<struct bin*> * binspack, struct bin * binX)
{
	for (auto m : lineGraph->neighbors(binX->node)) 
	{
		if (m <= binspack->size())
		{
			htd::vertex_t v = binX->node;
			htd::id_t edgeId = v < m ? *lineGraph->associatedEdgeIds(v, m).begin() : *lineGraph->associatedEdgeIds(m, v).begin();
			const htd::ILabel *c;
			c = &lineGraph->edgeLabel("edge", edgeId);
			struct bin *binM = binspack->at(m - 1);
			if (binM->partition != binX->partition)
			{
				if (*binM->l1n == *c)
				{
					if (binM->l1s == 0 && binM->l1o == 1) {
						binX->gain = binX->gain - 2;
					}
					else if (binM->l1s == 0 && binM->l1o > 1)
					{
						binX->gain = binX->gain - 1;
					}
					else if (binM->l1s > 0 && binM->l1o == 1)
					{
						binX->gain = binX->gain - 1;
					}
					else if (binM->l1s > 0 && binM->l1o > 1) {
						// stays the same
					}
				}
				else if (*binM->l2n == *c)
				{
					if (binM->l2s == 1 && binM->l2o == 0) {
						binX->gain = binX->gain - 2;
					}
					else if (binM->l2s == 1 && binM->l2o > 0)
					{
						binX->gain = binX->gain - 1;
					}
					else if (binM->l2s > 1 && binM->l2o == 0)
					{
						binX->gain = binX->gain - 1;
					}
					else if (binM->l2s > 1 && binM->l2o > 0) {
						// stays the same
					}
				}
			}
			else
			{
				if (*binM->l1n == *c)
				{
					if (binM->l1s == 1 && binM->l1o == 0)
					{
						binX->gain = binX->gain + 2;
					}
					else if (binM->l1s > 1 && binM->l1o == 0)
					{
						binX->gain = binX->gain + 1;
					}
					else if (binM->l1s == 1 && binM->l1o > 0)
					{
						binX->gain = binX->gain + 1;
					}
					else if (binM->l1s > 1 && binM->l1o > 0)
					{
						// stays the same
					}
				}
				else
				{
					if (binM->l2s == 1 && binM->l2o == 0)
					{
						binX->gain = binX->gain + 2;
					}
					else if (binM->l2s > 1 && binM->l2o == 0)
					{
						binX->gain = binX->gain + 1;
					}
					else if (binM->l2s == 1 && binM->l2o > 0)
					{
						binX->gain = binX->gain + 1;
					}
					else if (binM->l2s > 1 && binM->l2o > 0)
					{
						// stays the same
					}
				}
			}
		}
	}
}

void copyBinspackPairs(std::vector<bin *> * binspack, std::vector<bin *> * binspackCopy)
{
	for (int i = 0; i < binspack->size(); i++)
	{
		std::vector<struct pair> * newPairs = new std::vector<struct pair>();
		for (int j = 0; j < binspack->at(i)->pairs->size(); j++)
		{
			struct pair *pair = new struct pair;
			pair->lij = binspack->at(i)->pairs->at(j).lij;
			pair->nij = binspack->at(i)->pairs->at(j).nij;
			
			newPairs->push_back(*pair);
		}

		binspackCopy->at(i)->pairs = newPairs;
	}
}

void copyBinspackBinPointers(std::vector<bin *> * binspack, std::vector<bin *> * binspackCopy)
{
	for (int i = 0; i < binspack->size(); i++)
	{
		if (i < binspack->size() - 1 && binspack->at(i)->next != nullptr)
		{
			htd::index_t index = binspack->at(i)->next->node - 1;
			binspackCopy->at(i)->next = binspackCopy->at(index);
		} 
		else
		{
			binspackCopy->at(i)->next = nullptr;
		}

		if (i > 0 && binspack->at(i)->prev != nullptr)
		{
			htd::index_t index = binspack->at(i)->prev->node - 1;
			binspackCopy->at(i)->prev = binspackCopy->at(index);
		}
		else
		{
			binspackCopy->at(i)->prev = nullptr;
		}
	}
}

void copyBinspack(std::vector<bin *> * binspack, std::vector<bin *> * binspackCopy)
{
	for (int i = 0; i < binspack->size(); i++)
	{
		struct bin * bin = binspack->at(i);

		struct bin * newBin = new struct bin;
		newBin->degree = bin->degree;
		newBin->gain = bin->gain;
		newBin->l1n = bin->l1n;
		newBin->l2n = bin->l2n;
		newBin->l1s = bin->l1s;
		newBin->l1o = bin->l1o;
		newBin->l2s = bin->l2s;
		newBin->l2o = bin->l2o;
		newBin->locked = bin->locked;
		newBin->node = bin->node;
		newBin->partition = bin->partition;

		binspackCopy->push_back(newBin);
	}

	copyBinspackPairs(binspack, binspackCopy);
	copyBinspackBinPointers(binspack, binspackCopy);
}

void copyBinGainPartition(std::vector<binGainPair *> * binGainPairsPartition, std::vector<binGainPair *> * binGainPairsPartitionCopy, std::vector<bin *> * binspackCopy)
{
	for (int i = 0; i < binGainPairsPartition->size(); i++) 
	{
		struct binGainPair * bgPair = binGainPairsPartition->at(i);
		struct binGainPair * newBGpair = new struct binGainPair;
		newBGpair->gain = bgPair->gain;
		if (bgPair->bin != nullptr)
		{
			newBGpair->bin = binspackCopy->at(bgPair->bin->node - 1);
		}
		else
		{
			newBGpair->bin = nullptr;
		}

		binGainPairsPartitionCopy->push_back(newBGpair);
	}
}

void copyBinGainPartitions(std::vector<binGainPair *> * binGainPairsPartition1, std::vector<binGainPair *> * binGainPairsPartition1Copy, std::vector<binGainPair *> * binGainPairsPartition2, std::vector<binGainPair *> * binGainPairsPartition2Copy, std::vector<bin *> * binspackCopy)
{
	copyBinGainPartition(binGainPairsPartition1, binGainPairsPartition1Copy, binspackCopy);
	copyBinGainPartition(binGainPairsPartition2, binGainPairsPartition2Copy, binspackCopy);
}

void reCalculateCutLabels(std::vector<bin *> * binspack, std::vector<const htd::ILabel*> *cutLabels)
{
	cutLabels->clear();

	for (int i = 0; i < binspack->size(); i++)
	{
		struct bin * temp = binspack->at(i);

		if (temp->l1o > 0)
		{
			if (!contains(cutLabels, temp->l1n))
			{
				cutLabels->push_back(temp->l1n);
			}
		}
	}
}

void move(htd::LabeledMultiGraph * lineGraph, struct move move, std::vector<bin *> * binspack, std::vector<const htd::ILabel*> *cutLabels, std::vector<htd::vertex_t> *partition1, std::vector<htd::vertex_t> *partition2, std::vector<binGainPair *> * binGainPairsPartition1, std::vector<binGainPair *> * binGainPairsPartition2, int *w1, int *w2)
{
	std::vector<htd::vertex_t> *srcP = new std::vector<htd::vertex_t>();
	partition1->reserve(binspack->size());

	std::vector<htd::vertex_t> *dstP = new std::vector<htd::vertex_t>();
	partition2->reserve(binspack->size());

	int wSrc = 0, wDst = 0;

	if (*w1 > *w2)
	{
		srcP = partition1, wSrc = *w1;
		dstP = partition2, wDst = *w2;
	}
	else
	{
		srcP = partition2, wSrc = *w2;
		dstP = partition1, wDst = *w1;
	}

	struct bin * bin = move.bin;

	updateW1W2(&wSrc, &wDst, bin);

	if (w1 > w2)
	{
		*w1 = wSrc, *w2 = wDst;
	}
	else
	{
		*w2 = wSrc, *w1 = wDst;
	}

	bin->locked = true;

	updateTheGainsOfAdjacentNodes(lineGraph, *w1, *w2, binspack, bin);

	updateBins(binspack);

	updateBinPairPartitions(binGainPairsPartition1, binGainPairsPartition2, binspack);

	if (srcP->size() > 0)
	{
		srcP->erase(srcP->begin());
	}

	dstP->push_back(bin->node);

	adjustLabels(lineGraph, binspack, bin, partition1, partition2);

	reCalculateCutLabels(binspack, cutLabels);
}

void getBinWithMaxGain(std::vector<binGainPair *> * binGainPairsPartition2Copy, bin * &nodeWithMaxGain)
{
	for (int i = 4; i >= 0; i--)
	{
		if (binGainPairsPartition2Copy->at(i)->bin != nullptr) { nodeWithMaxGain = binGainPairsPartition2Copy->at(i)->bin; break; }
	}
}

void lgb(htd::LabeledMultiGraph * lineGraph, float eps, std::size_t numberOfnodes, std::vector<htd::vertex_t> *partition1, std::vector<htd::vertex_t> *partition2)
{
	int LIMIT = 10; // TODO implement logic for limit
	std::vector<struct bin*> * binspack = new std::vector<struct bin*>();
	std::vector<struct move> * moves = new std::vector<struct move>();
	
	std::vector<struct binGainPair *> * binGainPairsPartition1 = new std::vector<struct binGainPair *>();
	std::vector<struct binGainPair *> * binGainPairsPartition2 = new std::vector<struct binGainPair *>();
	
	initializePartitions(lineGraph, numberOfnodes, partition1, partition2, binspack);
	std::sort(binspack->begin(), binspack->end(), sortBins);
	
	computeLs(lineGraph, numberOfnodes, partition1, partition2, binspack);

	std::vector<const htd::ILabel*> *cutLabels = new std::vector<const htd::ILabel*>();
	std::vector<const htd::ILabel*> *part1Labels = new std::vector<const htd::ILabel*>();
	std::vector<const htd::ILabel*> *part2Labels = new std::vector<const  htd::ILabel*>();
	findDistinctLabelsInPartitions(cutLabels, part1Labels, part2Labels, binspack);

	int w1 = part1Labels->size(), 
		w2 = part2Labels->size(), 
		originalW1 = w1, 
		originalW2 = w2,
		loopNo = 0, 
		loadDiff = abs(originalW1 - originalW2), 
		prevLoadDiff = 0, 
		gainMax = 0;
	do
	{
		loopNo++;
		initBins(binspack, binGainPairsPartition1, binGainPairsPartition2);

		htd::LabeledMultiGraph * tempCopy = lineGraph->clone();
		std::vector<struct bin*> * binspackCopy = new std::vector<struct bin*>();
		std::vector<struct binGainPair *> * binGainPairsPartition1Copy = new std::vector<struct binGainPair *>();
		std::vector<struct binGainPair *> * binGainPairsPartition2Copy = new std::vector<struct binGainPair *>();
		copyBinspack(binspack, binspackCopy);
		copyBinGainPartitions(binGainPairsPartition1, binGainPairsPartition1Copy, binGainPairsPartition2, binGainPairsPartition2Copy, binspackCopy);

		calculateGainsAndUpdateBins(tempCopy, numberOfnodes, binspackCopy, binGainPairsPartition1Copy, binGainPairsPartition2Copy);

		prevLoadDiff = loadDiff;
		int seqno = 0;

		std::vector<htd::vertex_t> *srcP = new std::vector<htd::vertex_t>();
		partition1->reserve(numberOfnodes);

		std::vector<htd::vertex_t> *dstP = new std::vector<htd::vertex_t>();
		partition2->reserve(numberOfnodes);

		bool found = false;
		do {
			seqno++;

			int wSrc = 0, wDst = 0;

			if (w1 > w2)
			{
				srcP = partition1, wSrc = w1;
				dstP = partition2, wDst = w2;
			}
			else
			{
				srcP = partition2, wSrc = w2;
				dstP = partition1, wDst = w1;
			}

			int gain = 0;
			struct bin * binWithMaxGain = new struct bin;
			found = false;

			if (w1 > w2)
			{
				gain = binGainPairsPartition1Copy->at(binGainPairsPartition1Copy->size() - 1)->gain;
				getBinWithMaxGain(binGainPairsPartition1Copy, binWithMaxGain);
			}
			else
			{
				gain = binGainPairsPartition2Copy->at(binGainPairsPartition2Copy->size() - 1)->gain;
				getBinWithMaxGain(binGainPairsPartition2Copy, binWithMaxGain);
			}

			if (binWithMaxGain != nullptr && !binWithMaxGain->locked)
			{ 
				found = true; 
				gainMax = binWithMaxGain->gain;
			}

			if (found)
			{
				updateW1W2(&wSrc, &wDst, binWithMaxGain);

				if (w1 > w2) 
				{ 
					w1 = wSrc, w2 = wDst; 
				}
				else 
				{ 
					w2 = wSrc, w1 = wDst; 
				}

				loadDiff = abs(wSrc - wDst);

				struct move *move = new struct move;
				move->gain = binWithMaxGain->gain;
				move->seq = seqno;
				move->loadDiff = loadDiff;
				move->bin = binWithMaxGain;
				binWithMaxGain->locked = true;
				
				updateTheGainsOfAdjacentNodes(tempCopy, w1, w2, binspackCopy, binWithMaxGain);
				updateBins(binspackCopy);

				updateBinPairPartitions(binGainPairsPartition1Copy, binGainPairsPartition2Copy, binspackCopy);
				
				if (srcP->size() > 0)
				{
					srcP->erase(srcP->begin());	
				}
				dstP->push_back(binWithMaxGain->node);

				adjustLabels(tempCopy, binspackCopy, binWithMaxGain, partition1, partition2);

				reCalculateCutLabels(binspackCopy, cutLabels);
				move->cutSize = cutLabels->size();
				moves->push_back(*move);
			}
		} while (found);

		int k = findK(moves);

		w1 = originalW1, w2 = originalW2;

		if (moves->at(k).cutSize <= cutLabels->size() || moves->at(k).loadDiff <= loadDiff)
		{
			for (int i = 0; i <= k; i++)
			{
				move(lineGraph, moves->at(i), binspack, cutLabels, partition1, partition2, binGainPairsPartition1, binGainPairsPartition2, &w1, &w2);
			}
		}

		originalW1 = w1, originalW2 = w2;

	} while ((gainMax > 0 || (gainMax = 0 && prevLoadDiff > loadDiff)) && loopNo < LIMIT); // 10
}

std::vector<htd::vertex_t> * revertLineGraph(htd::LabeledMultiGraph * graph, htd::LabeledMultiGraph * lineGraph, std::vector<htd::vertex_t> *partition1, std::vector<htd::vertex_t> *partition2)
{
	std::vector<htd::vertex_t> * separator = new std::vector<htd::vertex_t>();

	for (auto v1 : *partition1)
	{
		for (auto v2 : *partition2)
		{
			if (lineGraph->isEdge(v1, v2) || lineGraph->isEdge(v2, v1))
			{
				htd::id_t edgeId = v1 < v2 ? *lineGraph->associatedEdgeIds(v1, v2).begin() : *lineGraph->associatedEdgeIds(v2, v1).begin();
				for (htd::vertex_t u : graph->vertices())
				{
					if (graph->vertexLabel("edge", u) == lineGraph->edgeLabel("edge", edgeId) && std::find(separator->begin(), separator->end(), u) == separator->end())
					{
						separator->push_back(u);
					}
				}
			}
		}
	}

	return separator;
}

std::vector<htd::vertex_t> * htd::LgbSeparatorAlgorithm::computeSeparator(htd::LabeledMultiGraph * graph) const
{

	htd::LabeledMultiGraph * lineGraph = new htd::LabeledMultiGraph(managementInstance(), graph->edgeCount());
	
	fillLineGraph(lineGraph, graph);

	std::size_t lineGraphNoArtifCount = graph->edgeCount(); // Number of vertices in the line graph not counting artificial nodes.

	std::vector<htd::vertex_t> *partition1 = new std::vector<htd::vertex_t>();
	partition1->reserve(lineGraphNoArtifCount);

	std::vector<htd::vertex_t> *partition2 = new std::vector<htd::vertex_t>();
	partition2->reserve(lineGraphNoArtifCount);
	
	float eps = 0.5f; // TODO use eps?

	lgb(lineGraph, eps, lineGraphNoArtifCount, partition1, partition2);

	return revertLineGraph(graph, lineGraph, partition1, partition2);
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

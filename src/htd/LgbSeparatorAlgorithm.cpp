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

//struct binGainPair
//{
//	int gain;
//	struct bin * binWithThatGain;
//};
//
struct nodeGainPair
{
	int gain;
	htd::index_t index;
};

struct move
{
	int seq;
	int gain;
	int loadDiff;
};

void initBins(std::vector<struct bin*> * binspack, std::vector<nodeGainPair> * nodeGainPairsPartition1, std::vector<nodeGainPair> * nodeGainPairsPartition2)
{
	struct bin * prev = new struct bin;
	
	prev->gain = 0;
	prev->node = -10;
	binspack->at(0)->prev = prev;
	
	for (int i = 0; i < binspack->size(); i++)
	{
		if (i != binspack->size() - 1)
		{
			binspack->at(i)->next = binspack->at(i + 1);
		}
		
		struct nodeGainPair * ng = new struct nodeGainPair;
		
		ng->gain = 0;
		ng->index = i;

		if (binspack->at(i)->partition == 1)
		{
			nodeGainPairsPartition1->push_back(*ng);
		}
		else
		{
			nodeGainPairsPartition2->push_back(*ng);
		}
	}
	
	return;
}

//void deleteBin(int position)
//{
//	int i = 0;
//	if (position<1)
//	{
//		printf("Invalid Position of Deletion \n");
//		return;
//	}
//	if (head == NULL)
//	{
//		return;
//	}
//	if (position == 1)
//	{
//		head = head->next;
//		if (head != NULL)
//		{
//			head->prev = NULL;
//		}
//	}
//	else
//	{
//		struct bin *temp = head;
//		while (temp->next->next != NULL && i<position - 2)
//		{
//			i++;
//			temp = temp->next;
//		}
//		if (i<position - 2)
//		{
//			printf("Invalid Position of Deletion \n");
//			return;
//		}
//		else
//		{
//			temp->next = temp->next->next;
//			if (temp->next != NULL)
//				temp->next->prev = temp;
//		}
//	}
//}

//struct bin * get(htd::vertex_t node)
//{
//	if (head == NULL)
//	{
//		return nullptr;
//	}
//	if (head->node == node)
//	{
//		return head;
//	}
//	else
//	{
//		struct bin *temp = head;
//		while (temp->next != NULL && temp->node != node)
//		{
//			temp = temp->next;
//		}
//
//		return temp;
//	}
//
//	return nullptr;
//}

//struct bin * at(int position)
//{
//	int i = 0;
//	if (position<1)
//	{
//		printf("Invalid Position of Deletion \n");
//		return nullptr;
//	}
//	if (head == NULL)
//	{
//		return nullptr;
//	}
//	if (position == 1)
//	{
//		return head;
//	}
//	else
//	{
//		struct bin *temp = head;
//		while (temp->next->next != NULL && i<position - 2)
//		{
//			i++;
//			temp = temp->next;
//		}
//		if (i<position - 2)
//		{
//			printf("Invalid Position of Deletion \n");
//			return nullptr;
//		}
//		else
//		{
//			return temp;
//		}
//	}
//
//	return nullptr;
//}

void printBins(std::vector<struct bin*> * binspack)
{
	if (binspack->size() == 0)
	{
		printf("Empty List!! \n");
		return;
	}

	for (int i = 0; i <= binspack->size(); i++)
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

	std::cout << std::floor(numberOfnodes / 2) << std::endl;

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
			int l1s = 0;
			int l1o = 0;
			int l2s = 0;
			int l2o = 0;
		
			const htd::ILabel *l1n = nullptr;
			const htd::ILabel *l2n = nullptr;

			std::vector<struct pair> *pairs = new std::vector<struct pair>(lineGraph->neighborCount(v));

			for (auto n : lineGraph->neighbors(v))
			{
				struct pair * pair = new struct pair;
				pair->nij = binspack->at(n);

				htd::id_t edgeId = v < n ? *lineGraph->associatedEdgeIds(v, n).begin() : *lineGraph->associatedEdgeIds(n, v).begin();
				if (l1n != nullptr && !(*l1n == lineGraph->edgeLabel("edge", edgeId)))
				{
					l2n = &lineGraph->edgeLabel("edge", edgeId);
					break;
				}
				l1n = &lineGraph->edgeLabel("edge", edgeId);

				pair->lij = &lineGraph->edgeLabel("edge", edgeId);
				pairs->push_back(*pair);
			}

			for (auto n : lineGraph->neighbors(v))
			{
				htd::id_t edgeId = v < n ? *lineGraph->associatedEdgeIds(v, n).begin() : *lineGraph->associatedEdgeIds(n, v).begin();

				bool inTheSamePartition = binspack->at(v)->partition == 1 ? (std::find(partition1->begin(), partition1->end(), n) != partition1->end()) : (std::find(partition2->begin(), partition2->end(), n) != partition2->end());

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

			binspack->at(v)->l1n = l1n;
			binspack->at(v)->l1s = l1s;
			binspack->at(v)->l1o = l1o;

			binspack->at(v)->l2n = l2n;
			binspack->at(v)->l2s = l2s;
			binspack->at(v)->l2o = l2o;
			binspack->at(v)->pairs = pairs;
		}
	}
}

//int weight(htd::LabeledMultiGraph * lineGraph, const htd::ILabel * label)
//{
//	int weight = 0;
//	for (auto v : lineGraph->vertices())
//	{
//		for (auto n : lineGraph->neighbors(v))
//		{
//			if (v < n)
//			{
//				htd::id_t edgeId = *lineGraph->associatedEdgeIds(v, n).begin();
//				if (*label == lineGraph->edgeLabel("edge", edgeId)) 
//				{
//					weight++;
//				}
//			}
//		}
//	}
//
//	return weight;
//}

void calculateGains(htd::LabeledMultiGraph * lineGraph, std::size_t numberOfnodes, std::vector<struct bin*> * binspack)
{
	for (auto v : lineGraph->vertices())
	{
		struct bin *bin = binspack->at(v);
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

void findDistinctLabelsInPartitions(std::vector<const htd::ILabel*> *cutLabels, std::vector<const htd::ILabel*> *part1Labels, std::vector<const htd::ILabel*> *part2Labels, std::vector<struct bin*> * binspack)
{
	for (int i = 0; i <= binspack->size(); i++)
	{
		struct bin * temp = binspack->at(i);

		if (temp->l1o > 0) 
		{
			cutLabels->push_back(temp->l1n);
		}
		else if (temp->l1s > 0)
		{
			if (temp->partition == 1) 
			{
				part1Labels->push_back(temp->l1n);
			} 
			else 
			{
				part2Labels->push_back(temp->l1n);
			}
		}

		if (temp->l2o > 0)
		{
			cutLabels->push_back(temp->l2n);
		}
		else if (temp->l2s > 0)
		{
			if (temp->partition == 1)
			{
				part1Labels->push_back(temp->l2n);
			}
			else
			{
				part2Labels->push_back(temp->l2n);
			}
		}

		temp = temp->next;
	}
	return;
 }

bool sortFunction(struct bin *bin1, struct bin *bin2)
{ 
	if (bin1 != nullptr && bin2 == nullptr)
	{
		return true;
	}
	else if (bin1 == nullptr && bin2 != nullptr)
	{
		return false;
	}
	else if (bin1 != nullptr && bin2 != nullptr)
	{
		return bin1->node < bin2->node; 
	}
	else
	{
		return false;
	}
}

bool sortNodeGainPairs(struct nodeGainPair ng1, struct nodeGainPair ng2) { return ng1.gain <= ng2.gain; }

void updateW1W2(int *wSrc, int *wDst, struct bin * binToBeMoved)
{
	if (binToBeMoved->l1o > 0)
	{
		if (binToBeMoved->l1s > 0)
		{
			// nothing changes
		}
		else
		{
			*wDst = *wDst + 1; // weight(L1n)
		}
	}
	else if (binToBeMoved->l1o = 0)
	{
		*wSrc = *wSrc - 1; // weight(L1n)
	}

	if (binToBeMoved->l2o > 0)
	{
		if (binToBeMoved->l2s > 0)
		{
			// nothing changes
		}
		else
		{
			*wDst = *wDst + 1; // weight(L2n)
		}
	}
	else if (binToBeMoved->l2o = 0)
	{
		*wSrc = *wSrc - 1; // weight(L2n)
	}

	return;
}

void moveInTemp(std::vector<struct bin*> * binspack, struct bin * nodeToBeMoved)
{

	return;
}

int findK(std::vector<struct move> * moves)
{
	int k = 0;
	int gain = 0;
	
	for (int i = 0; i < moves->size(); i++)
	{
		for (int j = 0; j < i; j++) {
			gain += moves->at(j).gain;
		}

		if (gain > k) 
		{
			k = gain;
		}
	}

	return k;
}

void updateTheGainsOfAdjacentNodes(htd::LabeledMultiGraph * lineGraph, int w1, int w2, std::vector<struct bin*> * binspack, struct bin * node)
{

	return;
}

void lgb(htd::LabeledMultiGraph * lineGraph, float eps, std::size_t numberOfnodes, std::vector<htd::vertex_t> *partition1, std::vector<htd::vertex_t> *partition2)
{
	int LIMIT = 10; // TODO implement logic for limit
	// Structures
	std::vector<struct bin*> * binspack = new std::vector<struct bin*>(numberOfnodes);
	
	std::vector<struct nodeGainPair> * nodeGainPairsPartition1 = new std::vector<struct nodeGainPair>();
	nodeGainPairsPartition1->reserve(numberOfnodes);
	
	std::vector<struct nodeGainPair> * nodeGainPairsPartition2 = new std::vector<struct nodeGainPair>();
	nodeGainPairsPartition2->reserve(numberOfnodes);


	std::vector<struct move> * moves = new std::vector<struct move>();
	// 1, 7
	initializePartitions(lineGraph, numberOfnodes, partition1, partition2, binspack);
	std::sort(binspack->begin(), binspack->end(), sortFunction);
	
	// 2
	computeLs(lineGraph, numberOfnodes, partition1, partition2, binspack);

	std::vector<const htd::ILabel*> *cutLabels = new std::vector<const htd::ILabel*>();
	std::vector<const htd::ILabel*> *part1Labels = new std::vector<const htd::ILabel*>();
	std::vector<const htd::ILabel*> *part2Labels = new std::vector<const  htd::ILabel*>();
	//3,4
	findDistinctLabelsInPartitions(cutLabels, part1Labels, part2Labels, binspack);

	// 5 - The weight are 1 atm.
	int w1 = part1Labels->size();
	int w2 = part2Labels->size();

	// 6
	int originalW1 = w1;
	int originalW2 = w2;

	// 7 is done already together with 1
	
	// 8
	int loopno = 0;
	int loadDiff = abs(originalW1 - originalW2);

	int prevLoadDiff = 0;
	int gainMax = 0;
	
	// 9
	do 
	{
		loopno++;
		
		// 9.1
		htd::LabeledMultiGraph * tempCopy = lineGraph->clone();
		
		// 9.2
		initBins(binspack, nodeGainPairsPartition1, nodeGainPairsPartition2);
		
		// Copies
		std::vector<struct bin*> * binspackCopy(binspack);
		std::vector<struct nodeGainPair> * nodeGainPairsPartition1Copy(nodeGainPairsPartition1);
		std::vector<struct nodeGainPair> * nodeGainPairsPartition2Copy(nodeGainPairsPartition2);

		// 9.3
		calculateGains(tempCopy, numberOfnodes, binspackCopy);
		
		// 9.4
		prevLoadDiff = loadDiff;
		int seqno = 0;
		
		// 9.5
		bool found = false;
		do {
			seqno++;

			std::vector<htd::vertex_t> *srcP = new std::vector<htd::vertex_t>();
			partition1->reserve(numberOfnodes);

			std::vector<htd::vertex_t> *dstP = new std::vector<htd::vertex_t>();
			partition2->reserve(numberOfnodes);

			int wSrc = 0;
			int wDst = 0;
			// 9.5.1
			if (w1 > w2)
			{
				srcP = partition1;
				wSrc = w1;

				dstP = partition2;
				wDst = w2;
			}
			else 
			{
				srcP = partition2;
				wSrc = w2;

				dstP = partition1;
				wDst = w1;
			}
			
			// 9.5.2
			std::sort(nodeGainPairsPartition1Copy->begin(), nodeGainPairsPartition1Copy->end(), sortNodeGainPairs);
			
			std::sort(nodeGainPairsPartition2Copy->begin(), nodeGainPairsPartition2Copy->end(), sortNodeGainPairs);
			
			std::sort(nodeGainPairsPartition1->begin(), nodeGainPairsPartition1->end(), sortNodeGainPairs);

			std::sort(nodeGainPairsPartition2->begin(), nodeGainPairsPartition2->end(), sortNodeGainPairs);

			int gain = 0;
			htd::index_t maxGainIndex = 0;
			if (w1 > w2)
			{
				//gainMax = nodeGainPairsPartition1->at(0).gain;
				gain = nodeGainPairsPartition1Copy->at(nodeGainPairsPartition1Copy->size() - 1).gain;
				maxGainIndex = nodeGainPairsPartition1Copy->at(0).index;
			}
			else
			{
				//gainMax = nodeGainPairsPartition2->at(0).gain;
				gain = nodeGainPairsPartition2Copy->at(nodeGainPairsPartition2Copy->size() - 1).gain;
				maxGainIndex = nodeGainPairsPartition2Copy->at(0).index;
			}

			struct bin * nodeWithMaxGain = binspackCopy->at(maxGainIndex);
			// 9.5.3
			if (nodeWithMaxGain->node != NULL && !nodeWithMaxGain->locked)
			{
				found = true;
			}

			// 9.5.3.1
			if (found)
			{
				// 9.5.3.2
				updateW1W2(&wSrc, &wDst, nodeWithMaxGain);
				if (w1 > w2)
				{
					w1 = wSrc;
					w2 = wDst;
				}
				else
				{
					w2 = wSrc;
					w1 = wDst;
				}
				loadDiff = abs(wSrc - wDst);

				// 9.5.3.3
				struct move *move = new struct move;
				move->gain = nodeWithMaxGain->gain;
				move->seq = seqno;
				move->loadDiff = loadDiff;
				moves->push_back(*move);

				updateTheGainsOfAdjacentNodes(tempCopy, w1, w2, binspackCopy, nodeWithMaxGain);
				moveInTemp(binspackCopy, nodeWithMaxGain); // TODO
			}
		} while (found);

		int k = findK(moves);

		w1 = originalW1;
		w2 = originalW2;

		// TODO 9.9

		originalW1 = w1;
		originalW2 = w2;
	} while ((gainMax > 0 || (gainMax = 0 && prevLoadDiff > loadDiff)) && loopno < LIMIT); // 10

	printBins(binspack);
	
	return;
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
				for (htd::vertex_t u : graph->vertices)
				{
					if (graph->vertexLabel("edge", u) == lineGraph->edgeLabel("edge", edgeId))
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
	std::cout << "LINE GRAPH VERTICES:  " << lineGraph->vertexCount() << std::endl;
	for (auto v1 : lineGraph->vertices())
	{
		for (auto v2 : lineGraph->neighbors(v1))
		{
			if (lineGraph->isEdge(v1, v2)) {
				std::cout << "" <<std::endl;
				std::cout << v1 << " -  " << v2 << std::endl;
				std::ostream stream(nullptr);
				stream.rdbuf(std::cout.rdbuf());
				std::cout << "Label ";
				lineGraph->edgeLabel("edge", *lineGraph->associatedEdgeIds(v1, v2).begin()).print(stream);
				std::cout << "" << std::endl;

			}
		}
	}

	std::size_t lineGraphNoArtifCount = graph->edgeCount(); // Number of vertices in the line graph not counting artificial nodes.

	std::vector<htd::vertex_t> *partition1 = new std::vector<htd::vertex_t>();
	partition1->reserve(lineGraphNoArtifCount);

	std::vector<htd::vertex_t> *partition2 = new std::vector<htd::vertex_t>();
	partition2->reserve(lineGraphNoArtifCount);
	
	float eps = 0.5f;

	lgb(lineGraph, eps, lineGraphNoArtifCount, partition1, partition2);
	system("pause");
	return revertLineGraph(graph, lineGraph, partition1, partition2);;
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

/*
 * File:   TreeDecompositionViaSeparatorAlgorithmTest.cpp
 *
 * Author: ABSEHER Michael (abseher@dbai.tuwien.ac.at)
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

#include <gtest/gtest.h>

#include <htd/main.hpp>

#include <vector>

class TreeDecompositionViaSeparatorAlgorithmTest : public ::testing::Test
{
    public:
		TreeDecompositionViaSeparatorAlgorithmTest(void)
        {

        }

        virtual ~TreeDecompositionViaSeparatorAlgorithmTest()
        {

        }

        void SetUp()
        {

        }

        void TearDown()
        {

        }
};

bool isValidTreeDecomposition(const htd::IMultiHypergraph & graph, const htd::ITreeDecomposition & decomposition)
{
	bool ret = true;
	
	for (htd::vertex_t vertex : decomposition.vertices())
	{
		const std::vector<htd::vertex_t> & bag = decomposition.bagContent(vertex);
		const htd::FilteredHyperedgeCollection & inducedHyperedges = decomposition.inducedHyperedges(vertex);

		std::unordered_set<htd::vertex_t> missingVertices(bag.begin(), bag.end());
		for (auto it = inducedHyperedges.begin(); it != inducedHyperedges.end() && !missingVertices.empty(); ++it)
		{
			const htd::Hyperedge & retrievedEdge = graph.hyperedge(it->id());

			ret = ret && retrievedEdge.elements() == it->elements();
			ret = ret && retrievedEdge.sortedElements() == it->sortedElements();
			std::cout << "22*********HERE***********22" << std::endl;
			for (htd::vertex_t vertex2 : it->sortedElements())
			{
				missingVertices.erase(vertex2);
			}
		}		
		while (!missingVertices.empty())
		{
			htd::vertex_t missingVertex = *(missingVertices.begin());

			ret = ret && graph.edgeCount(missingVertex) == 0 && graph.isIsolatedVertex(missingVertex);
			
			missingVertices.erase(missingVertex);
		}
		
	}
	std::cout << "PROSLO" << std::endl;
	return ret;
}


TEST(TreeDecompositionViaSeparatorTest, CheckDecomposition)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

	htd::MultiGraph & graph = * new htd::MultiGraph(libraryInstance, 11);

	graph.addEdge(1, 2);
	graph.addEdge(1, 3);
	graph.addEdge(2, 3);
	graph.addEdge(3, 4);
	graph.addEdge(4, 5);
	graph.addEdge(3, 5);
	graph.addEdge(5, 6);
	graph.addEdge(2, 6);
	graph.addEdge(1, 4);
	graph.addEdge(6, 7);
	graph.addEdge(7, 8);
	graph.addEdge(7, 9);
	graph.addEdge(6, 9);
	graph.addEdge(8, 9);
	graph.addEdge(6, 10);
	graph.addEdge(8, 10);
	graph.addEdge(1, 8);
	graph.addEdge(1, 11);
	graph.addEdge(2, 11);

	htd::TreeDecompositionViaSeparatorAlgorithm algorithm(libraryInstance);

	htd::ITreeDecomposition & decomposition = *algorithm.computeDecomposition(graph);

	htd::TreeDecompositionVerifier *tdv = new htd::TreeDecompositionVerifier();

	ASSERT_TRUE(tdv->verify(graph, decomposition));
 
    delete libraryInstance;
}

TEST(TreeDecompositionViaSeparatorTest, CheckDecomposition2)
{
	htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

	htd::MultiGraph & graph = *new htd::MultiGraph(libraryInstance, 25);

	graph.addEdge(1, 2);
	graph.addEdge(1, 7);
	graph.addEdge(2, 7);
	graph.addEdge(3, 4);
	graph.addEdge(3, 9);
	graph.addEdge(4, 5);
	graph.addEdge(4, 10);
	graph.addEdge(5, 10);
	graph.addEdge(6, 11);
	graph.addEdge(6, 12);
	graph.addEdge(7, 8);
	graph.addEdge(7, 12);
	graph.addEdge(7, 13);
	graph.addEdge(8, 9);
	graph.addEdge(8, 13);
	graph.addEdge(8, 14);
	graph.addEdge(9, 14);
	graph.addEdge(9, 15);
	graph.addEdge(11, 17);
	graph.addEdge(11, 16);
	graph.addEdge(12, 13);
	graph.addEdge(12, 18);
	graph.addEdge(13, 14);
	graph.addEdge(13, 18);
	graph.addEdge(14, 15);
	graph.addEdge(15, 20);
	graph.addEdge(16, 21);
	graph.addEdge(17, 18);
	graph.addEdge(17, 22);
	graph.addEdge(18, 19);
	graph.addEdge(18, 23);		
	graph.addEdge(19, 24);
	graph.addEdge(19, 25);
	graph.addEdge(20, 25);
	graph.addEdge(21, 22);
	graph.addEdge(22, 23);
	graph.addEdge(23, 24);

	htd::TreeDecompositionViaSeparatorAlgorithm algorithm(libraryInstance);
	algorithm.setAlgorithmType(1);
	algorithm.setCriteriaType(3);
	htd::ITreeDecomposition & decomposition = *algorithm.computeDecomposition(graph);

	htd::TreeDecompositionVerifier *tdv = new htd::TreeDecompositionVerifier();

	ASSERT_TRUE(tdv->verifyConnectednessCriterion(graph, decomposition));

	delete libraryInstance;
}

TEST(TreeDecompositionViaSeparatorTest, CheckDecompositio3)
{
	htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

	htd::MultiGraph & graph = *new htd::MultiGraph(libraryInstance, 25);

	graph.addEdge(1, 2);
	graph.addEdge(2, 3);
	graph.addEdge(2, 7);
	graph.addEdge(3, 4);
	graph.addEdge(4, 5);
	graph.addEdge(4, 23);
	graph.addEdge(5, 6);
	graph.addEdge(5, 20);
	graph.addEdge(6, 7);
	graph.addEdge(6, 13);
	graph.addEdge(7, 8);
	graph.addEdge(7, 12);
	graph.addEdge(8, 9);
	graph.addEdge(9, 10);
	graph.addEdge(10, 11);
	graph.addEdge(11, 12);
	graph.addEdge(13, 14);
	graph.addEdge(13, 18);
	graph.addEdge(14, 15);
	graph.addEdge(15, 16);
	graph.addEdge(16, 17);
	graph.addEdge(17, 18);
	graph.addEdge(19, 20);
	graph.addEdge(20, 21);
	graph.addEdge(20, 25);
	graph.addEdge(21, 22);
	graph.addEdge(22, 23);
	graph.addEdge(23, 24);
	graph.addEdge(24, 25);

	htd::TreeDecompositionViaSeparatorAlgorithm algorithm(libraryInstance);
	algorithm.setAlgorithmType(2);
	algorithm.setCriteriaType(3);
	htd::ITreeDecomposition & decomposition = *algorithm.computeDecomposition(graph);

	htd::TreeDecompositionVerifier *tdv = new htd::TreeDecompositionVerifier();

	ASSERT_TRUE(tdv->verify(graph, decomposition));

	delete libraryInstance;
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}

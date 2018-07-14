/*
 * File:   SeparatorAlgorithmTest.cpp
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

class SeparatorAlgorithmTest : public ::testing::Test
{
    public:
		SeparatorAlgorithmTest(void)
        {

        }

        virtual ~SeparatorAlgorithmTest()
        {

        }

        void SetUp()
        {

        }

        void TearDown()
        {

        }
};

TEST(SeparatorAlgorithmTest, CompareAlgorithms)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

	std::vector<htd::vertex_t> * hdSeparator = nullptr;
	std::vector<htd::vertex_t> * lgbSeparator = nullptr;

	htd::MultiGraph * graph = new htd::MultiGraph(libraryInstance, 9);
	graph->addEdge(1, 2);
	graph->addEdge(2, 3);
	graph->addEdge(3, 4);
	graph->addEdge(3, 6);
	graph->addEdge(3, 7);
	graph->addEdge(4, 5);
	graph->addEdge(6, 9);
	graph->addEdge(7, 8);

	htd::HighestDegreeSeparatorAlgorithm * separatorAlgorithm = new htd::HighestDegreeSeparatorAlgorithm(libraryInstance);
	hdSeparator = separatorAlgorithm->computeSeparator(graph);

	htd::LgbSeparatorAlgorithm * separatorAlgorithm = new htd::LgbSeparatorAlgorithm(libraryInstance);
	lgbSeparator = separatorAlgorithm->computeSeparator(graph);

	ASSERT_TRUE(hdSeparator > lgbSeparator);
	ASSERT_EQ((std::size_t)1, hdSeparator->size());
	ASSERT_EQ((std::size_t)1, lgbSeparator->size());
	
	delete hdSeparator;
	delete lgbSeparator;
	delete separatorAlgorithm;
    delete libraryInstance;
}

TEST(SeparatorAlgorithmTest, CompareAlgorithms2)
{
	htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

	std::vector<htd::vertex_t> * hdSeparator = nullptr;
	std::vector<htd::vertex_t> * lgbSeparator = nullptr;

	htd::MultiGraph * graph = new htd::MultiGraph(libraryInstance, 9);
	graph->addEdge(1, 2);
	graph->addEdge(1, 4);
	graph->addEdge(3, 5);
	graph->addEdge(3, 6);
	graph->addEdge(4, 5);
	graph->addEdge(4, 7);
	graph->addEdge(5, 7);
	graph->addEdge(5, 8);
	graph->addEdge(6, 8);
	graph->addEdge(6, 9);

	htd::HighestDegreeSeparatorAlgorithm * separatorAlgorithm = new htd::HighestDegreeSeparatorAlgorithm(libraryInstance);
	hdSeparator = separatorAlgorithm->computeSeparator(graph);

	htd::LgbSeparatorAlgorithm * separatorAlgorithm = new htd::LgbSeparatorAlgorithm(libraryInstance);
	lgbSeparator = separatorAlgorithm->computeSeparator(graph);

	ASSERT_TRUE(hdSeparator > lgbSeparator);
	ASSERT_EQ((std::size_t)1, hdSeparator->size());
	ASSERT_EQ((std::size_t)1, lgbSeparator->size());

	delete hdSeparator;
	delete lgbSeparator;
	delete separatorAlgorithm;
	delete libraryInstance;
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}

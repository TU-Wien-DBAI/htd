#include <gtest/gtest.h>

#include <htd/main.hpp>
#include <htd/FitnessMinimizingTreeDecompositionAlgorithm.hpp>

class FitnessMinimizingTreeDecompositionAlgorithmTest : public ::testing::Test
{
    public:
        FitnessMinimizingTreeDecompositionAlgorithmTest(void)
        {

        }

        virtual ~FitnessMinimizingTreeDecompositionAlgorithmTest()
        {

        }

        void SetUp()
        {

        }

        void TearDown()
        {

        }
};

class WidthFitnessFunction : public htd::ITreeDecompositionFitnessFunction
{
    public:
        WidthFitnessFunction(void)
        {

        }

        ~WidthFitnessFunction()
        {

        }

        htd::FitnessEvaluation * fitness(const htd::IMultiHypergraph & graph, const htd::ITreeDecomposition & decomposition) const
        {
            HTD_UNUSED(graph)

            return new htd::FitnessEvaluation(1, -(double)(decomposition.maximumBagSize()));
        }

        WidthFitnessFunction * clone(void) const
        {
            return new WidthFitnessFunction();
        }
};

TEST(FitnessMinimizingTreeDecompositionAlgorithmTest,)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

    htd::MultiHypergraph graph(libraryInstance);

    htd::vertex_t vertex1 = graph.addVertex();
    htd::vertex_t vertex2 = graph.addVertex();
    htd::vertex_t vertex3 = graph.addVertex();

    graph.addEdge(vertex1, vertex2);
    graph.addEdge(vertex2, vertex3);

    std::string decompString1 = std::string("s td 2 2 3\n"
                                            "b 1 2 3 \n"
                                            "b 2 1 2 \n"
                                            "1 2");

    htd::FileTreeDecompositionAlgorithm fileAlgorithm1(libraryInstance, decompString1);

    std::string decompString2 = std::string("s td 1 3 3\n"
                                            "b 1 1 2 3");

    htd::FileTreeDecompositionAlgorithm fileAlgorithm2(libraryInstance, decompString2);

    std::vector<htd::ITreeDecompositionAlgorithm *> detAlgorithms;
    detAlgorithms.push_back(&fileAlgorithm1);
    detAlgorithms.push_back(&fileAlgorithm2);

    htd::FitnessMinimizingTreeDecompositionAlgorithm algorithm(libraryInstance, new WidthFitnessFunction(), detAlgorithms, std::vector<htd::ITreeDecompositionAlgorithm *>());

    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph);
    ASSERT_EQ(decomposition->maximumBagSize(), (std::size_t)2);
}

int main(int argc, char ** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}

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

TEST(EmptyAlgorithm,)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

    htd::MultiHypergraph graph(libraryInstance);

    htd::vertex_t vertex1 = graph.addVertex();
    htd::vertex_t vertex2 = graph.addVertex();
    htd::vertex_t vertex3 = graph.addVertex();

    graph.addEdge(vertex1, vertex2);
    graph.addEdge(vertex2, vertex3);

    htd::FitnessMinimizingTreeDecompositionAlgorithm algorithm(libraryInstance, new WidthFitnessFunction(), std::vector<htd::ITreeDecompositionAlgorithm *>(), std::vector<htd::ITreeDecompositionAlgorithm *>());

    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph);
    ASSERT_EQ(decomposition, nullptr);
}

TEST(FitnessMinimizingTreeDecompositionAlgorithmTestStatic,)
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

TEST(FitnessMinimizingTreeDecompositionAlgorithmTestDynamic,)
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

    std::vector<htd::ITreeDecompositionAlgorithm *> dynAlgorithms;
    dynAlgorithms.push_back(&fileAlgorithm1);
    dynAlgorithms.push_back(&fileAlgorithm2);

    htd::FitnessMinimizingTreeDecompositionAlgorithm algorithm(libraryInstance, new WidthFitnessFunction(), std::vector<htd::ITreeDecompositionAlgorithm *>(), dynAlgorithms);

    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph);
    ASSERT_EQ(decomposition->maximumBagSize(), (std::size_t)2);
}

TEST(FitnessMinimizingTreeDecompositionAlgorithmTestDynamicAndStatic,)
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

    std::vector<htd::ITreeDecompositionAlgorithm *> dynAlgorithms;
    dynAlgorithms.push_back(&fileAlgorithm1);
    std::vector<htd::ITreeDecompositionAlgorithm *> staticAlgorithms;
    staticAlgorithms.push_back(&fileAlgorithm2);

    htd::FitnessMinimizingTreeDecompositionAlgorithm algorithm(libraryInstance, new WidthFitnessFunction(), staticAlgorithms, dynAlgorithms);

    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph);
    ASSERT_EQ(decomposition->maximumBagSize(), (std::size_t)2);
}

namespace htd
{
    int id_DynamicTestAlgorithm = 0;

    class DynamicTestAlgorithm : public htd::ITreeDecompositionAlgorithm
    {
        public:

            std::vector<htd::ITreeDecompositionAlgorithm *> algorithms;

            HTD_API DynamicTestAlgorithm()
            {

            }

            HTD_API DynamicTestAlgorithm(const htd::DynamicTestAlgorithm & original)
            {
                algorithms = original.algorithms;
            }

            HTD_API ~DynamicTestAlgorithm()
            {};

            HTD_API htd::ITreeDecomposition * computeDecomposition(const htd::IMultiHypergraph & graph) const HTD_OVERRIDE
            {
                ITreeDecomposition * decomp = algorithms[id_DynamicTestAlgorithm]->computeDecomposition(graph);
                id_DynamicTestAlgorithm++;
                return decomp;
            }

            HTD_API htd::ITreeDecomposition * computeDecomposition(const htd::IMultiHypergraph & graph, const std::vector<htd::IDecompositionManipulationOperation *> &) const HTD_OVERRIDE
            {
                return computeDecomposition(graph);
            }

            HTD_API htd::ITreeDecomposition * computeDecomposition(const htd::IMultiHypergraph & graph, int, ...) const
            {
                return computeDecomposition(graph);
            }

            HTD_API htd::ITreeDecomposition * computeDecomposition(const htd::IMultiHypergraph & graph, std::size_t, std::size_t) const
            {
                return computeDecomposition(graph);
            }

            HTD_API htd::ITreeDecomposition * computeDecomposition(const htd::IMultiHypergraph & graph, const std::vector<htd::IDecompositionManipulationOperation *> &, std::size_t, std::size_t) const
            {
                return computeDecomposition(graph);
            }

            HTD_API htd::ITreeDecomposition * computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph &) const HTD_OVERRIDE
            {
                return computeDecomposition(graph);
            }

            HTD_API htd::ITreeDecomposition * computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph &, const std::vector<htd::IDecompositionManipulationOperation *> &) const HTD_OVERRIDE
            {
                return computeDecomposition(graph);
            }

            HTD_API htd::ITreeDecomposition * computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph &, int, ...) const
            {
                return computeDecomposition(graph);
            }

            HTD_API htd::ITreeDecomposition * computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph &, std::size_t, std::size_t) const
            {
                return computeDecomposition(graph);
            }

            HTD_API htd::ITreeDecomposition * computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph &, const std::vector<htd::IDecompositionManipulationOperation *> &, std::size_t, std::size_t) const
            {
                return computeDecomposition(graph);
            }

            HTD_API void setManipulationOperations(const std::vector<htd::IDecompositionManipulationOperation *> &) HTD_OVERRIDE
            {

            }

            HTD_API void addManipulationOperation(htd::IDecompositionManipulationOperation *) HTD_OVERRIDE
            {

            }

            HTD_API void addManipulationOperations(const std::vector<htd::IDecompositionManipulationOperation *> &) HTD_OVERRIDE
            {

            }

            HTD_API bool isSafelyInterruptible(void) const HTD_OVERRIDE
            {
                return false;
            }

            HTD_API bool isComputeInducedEdgesEnabled(void) const HTD_OVERRIDE
            {
                return false;
            }

            HTD_API void setComputeInducedEdgesEnabled(bool) HTD_OVERRIDE
            {

            }

            HTD_API const htd::LibraryInstance * managementInstance(void) const HTD_NOEXCEPT HTD_OVERRIDE
            {
                return nullptr;
            }

            HTD_API void setManagementInstance(const htd::LibraryInstance * const) HTD_OVERRIDE
            {

            }

            HTD_API DynamicTestAlgorithm * clone(void) const HTD_OVERRIDE
            {
                return new DynamicTestAlgorithm(*this);
            }
    };
}

TEST(TestDynamicAlgorithm,)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

    htd::MultiHypergraph graph(libraryInstance);

    htd::vertex_t vertex1 = graph.addVertex();
    htd::vertex_t vertex2 = graph.addVertex();
    htd::vertex_t vertex3 = graph.addVertex();
    htd::vertex_t vertex4 = graph.addVertex();
    htd::vertex_t vertex5 = graph.addVertex();

    graph.addEdge(vertex1, vertex2);
    graph.addEdge(vertex2, vertex3);
    graph.addEdge(vertex3, vertex4);
    graph.addEdge(vertex4, vertex5);

    htd::FileTreeDecompositionAlgorithm fileAlgorithm1(libraryInstance, std::string("s td 1 5 5\n"
                                                                                    "b 1 1 2 3 4 5"));

    std::vector<htd::ITreeDecompositionAlgorithm *> staticAlgorithms;
    staticAlgorithms.push_back(&fileAlgorithm1);

    htd::FileTreeDecompositionAlgorithm fileAlgorithm2(libraryInstance, std::string("s td 3 3 5\n"
                                                                                    "b 1 1 2 3\n"
                                                                                    "b 2 3 4\n"
                                                                                    "b 3 4 5\n"));

    htd::FileTreeDecompositionAlgorithm fileAlgorithm3(libraryInstance, std::string("s td 4 2 5\n"
                                                                                    "b 1 1 2\n"
                                                                                    "b 2 2 3\n"
                                                                                    "b 3 3 4\n"
                                                                                    "b 4 4 5\n"));

    htd::FileTreeDecompositionAlgorithm fileAlgorithm4(libraryInstance, std::string("s td 1 3 3\n"
                                                                                    "b 1 1 2 3"));

    htd::FileTreeDecompositionAlgorithm fileAlgorithm5(libraryInstance, std::string("s td 3 3 5\n"
                                                                                    "b 1 1 2 3\n"
                                                                                    "b 2 2 3\n"
                                                                                    "b 3 2 3\n"));

    htd::DynamicTestAlgorithm dynAlgorithm;
    dynAlgorithm.algorithms.push_back(&fileAlgorithm2);
    dynAlgorithm.algorithms.push_back(&fileAlgorithm3);
    dynAlgorithm.algorithms.push_back(&fileAlgorithm4);
    dynAlgorithm.algorithms.push_back(&fileAlgorithm5);

    std::vector<htd::ITreeDecompositionAlgorithm *> dynAlgorithms;
    dynAlgorithms.push_back(&dynAlgorithm);

    htd::FitnessMinimizingTreeDecompositionAlgorithm algorithm(libraryInstance, new WidthFitnessFunction(), staticAlgorithms, dynAlgorithms);

    htd::id_DynamicTestAlgorithm = 0;
    algorithm.setIterationCount(0);
    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph);
    ASSERT_EQ(decomposition->maximumBagSize(), (std::size_t)5);
    delete decomposition;

    htd::id_DynamicTestAlgorithm = 0;
    algorithm.setIterationCount(1);
    decomposition = algorithm.computeDecomposition(graph);
    ASSERT_EQ(decomposition->maximumBagSize(), (std::size_t)3);
    delete decomposition;

    htd::id_DynamicTestAlgorithm = 0;
    algorithm.setIterationCount(2);
    decomposition = algorithm.computeDecomposition(graph);
    ASSERT_EQ(decomposition->maximumBagSize(), (std::size_t)2);
    delete decomposition;

    htd::id_DynamicTestAlgorithm = 0;
    algorithm.setIterationCount(3);
    decomposition = algorithm.computeDecomposition(graph);
    ASSERT_EQ(decomposition->maximumBagSize(), (std::size_t)2);
    delete decomposition;

    htd::id_DynamicTestAlgorithm = 0;
    algorithm.setIterationCount(4);
    decomposition = algorithm.computeDecomposition(graph);
    ASSERT_EQ(decomposition->maximumBagSize(), (std::size_t)2);
    delete decomposition;
}

int main(int argc, char ** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}

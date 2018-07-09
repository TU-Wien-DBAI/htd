#include <gtest/gtest.h>

#include <htd/main.hpp>
#include <htd_io/TdFormatExporter.hpp>

class ExternalPipeTreeDecompositionAlgorithmTest : public ::testing::Test
{
    public:
        ExternalPipeTreeDecompositionAlgorithmTest(void)
        {

        }

        virtual ~ExternalPipeTreeDecompositionAlgorithmTest()
        {

        }

        void SetUp()
        {

        }

        void TearDown()
        {

        }
};

#ifndef BUILD_EXTERNAL_SOLVERS
#ifdef MRPRAJESH_PATH

TEST(Mrprajesh, TestMrprajesh)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

    htd::MultiHypergraph graph(libraryInstance);

    std::string path(MRPRAJESH_PATH);

    htd::ExternalPipeTreeDecompositionAlgorithm algorithm(libraryInstance, path + "/tw-heuristic", 0);

    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph);

    ASSERT_NE(decomposition, nullptr);

    ASSERT_EQ(decomposition->vertexCount(), (std::size_t)1);

    EXPECT_EQ(decomposition->edgeCount(), (std::size_t)0);

    ASSERT_EQ(decomposition->root(), (htd::vertex_t)1);

    const std::vector<htd::vertex_t> & bag = decomposition->bagContent(1);

    EXPECT_EQ(bag.size(), (std::size_t)0);

    delete decomposition;

    delete libraryInstance;
}

TEST(SimpleGraph_Mrprajesh, CheckResultSimpleGraph_Mrprajesh)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

    htd::MultiHypergraph graph(libraryInstance);

    std::string path(MRPRAJESH_PATH);

    htd::vertex_t vertex1 = graph.addVertex();
    htd::vertex_t vertex2 = graph.addVertex();
    htd::vertex_t vertex3 = graph.addVertex();

    graph.addEdge(vertex1, vertex2);
    graph.addEdge(vertex2, vertex3);

    htd::ExternalPipeTreeDecompositionAlgorithm algorithm(libraryInstance, path + "/tw-heuristic", 0);

    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph);

    ASSERT_NE(decomposition, nullptr);

    ASSERT_GE(decomposition->vertexCount(), (std::size_t)1);

    EXPECT_EQ(decomposition->edgeCount(), decomposition->vertexCount() - 1);

    ASSERT_EQ(decomposition->root(), (htd::vertex_t)1);

    ASSERT_LE(decomposition->minimumBagSize(), decomposition->maximumBagSize());

    htd::TreeDecompositionVerifier verifier;

    ASSERT_TRUE(verifier.verify(graph, *decomposition));

    delete decomposition;

    delete libraryInstance;
}

#endif
#ifdef BZTREEWIDTH_PATH

TEST(BZTreewidth, TestBZTreewidth)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

    htd::MultiHypergraph graph(libraryInstance);

    std::string path(BZTREEWIDTH_PATH);

    htd::ExternalPipeTreeDecompositionAlgorithm algorithm(libraryInstance, path + "/bin/BZTreewidth-DFS.exe", 0);

    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph);

    ASSERT_NE(decomposition, nullptr);

    ASSERT_EQ(decomposition->vertexCount(), (std::size_t)1);

    EXPECT_EQ(decomposition->edgeCount(), (std::size_t)0);

    ASSERT_EQ(decomposition->root(), (htd::vertex_t)1);

    const std::vector<htd::vertex_t> & bag = decomposition->bagContent(1);

    EXPECT_EQ(bag.size(), (std::size_t)0);

    delete decomposition;

    delete libraryInstance;
}

TEST(SimpleGraph_BZTreewidth, CheckResultSimpleGraph_BZTreewidth)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

    htd::MultiHypergraph graph(libraryInstance);

    std::string path(BZTREEWIDTH_PATH);

    htd::vertex_t vertex1 = graph.addVertex();
    htd::vertex_t vertex2 = graph.addVertex();
    htd::vertex_t vertex3 = graph.addVertex();

    graph.addEdge(vertex1, vertex2);
    graph.addEdge(vertex2, vertex3);

    htd::ExternalPipeTreeDecompositionAlgorithm algorithm(libraryInstance, path + "/bin/BZTreewidth-DFS.exe", 0);

    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph);

    ASSERT_NE(decomposition, nullptr);

    ASSERT_GE(decomposition->vertexCount(), (std::size_t)1);

    EXPECT_EQ(decomposition->edgeCount(), decomposition->vertexCount() - 1);

    ASSERT_EQ(decomposition->root(), (htd::vertex_t)1);

    ASSERT_LE(decomposition->minimumBagSize(), decomposition->maximumBagSize());

    htd::TreeDecompositionVerifier verifier;

    ASSERT_TRUE(verifier.verify(graph, *decomposition));

    delete decomposition;

    delete libraryInstance;
}

#endif
#ifdef JDRASIL_PATH

TEST(Jdrasil, TestJdrasil)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

    htd::MultiHypergraph graph(libraryInstance);

    std::string path(JDRASIL_PATH);

    htd::ExternalPipeTreeDecompositionAlgorithm algorithm(libraryInstance, path + "/tw-heuristic", 0);

    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph);

    ASSERT_NE(decomposition, nullptr);

    ASSERT_EQ(decomposition->vertexCount(), (std::size_t)1);

    EXPECT_EQ(decomposition->edgeCount(), (std::size_t)0);

    ASSERT_EQ(decomposition->root(), (htd::vertex_t)1);

    const std::vector<htd::vertex_t> & bag = decomposition->bagContent(1);

    EXPECT_EQ(bag.size(), (std::size_t)0);

    delete decomposition;

    delete libraryInstance;
}

TEST(SimpleGraph_Jdrasil, CheckResultSimpleGraph_Jdrasil)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

    htd::MultiHypergraph graph(libraryInstance);

    std::string path(JDRASIL_PATH);

    htd::vertex_t vertex1 = graph.addVertex();
    htd::vertex_t vertex2 = graph.addVertex();
    htd::vertex_t vertex3 = graph.addVertex();

    graph.addEdge(vertex1, vertex2);
    graph.addEdge(vertex2, vertex3);

    htd::ExternalPipeTreeDecompositionAlgorithm algorithm(libraryInstance, path + "/tw-heuristic", 0);

    algorithm.setDirectory(path);

    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph);

    ASSERT_NE(decomposition, nullptr);

    ASSERT_GE(decomposition->vertexCount(), (std::size_t)1);

    EXPECT_EQ(decomposition->edgeCount(), decomposition->vertexCount() - 1);

    ASSERT_EQ(decomposition->root(), (htd::vertex_t)1);

    ASSERT_LE(decomposition->minimumBagSize(), decomposition->maximumBagSize());

    htd::TreeDecompositionVerifier verifier;

    ASSERT_TRUE(verifier.verify(graph, *decomposition));

    delete decomposition;

    delete libraryInstance;
}

#endif
#ifdef MFJONES_PATH

TEST(mfjones, Testmfjones)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

    htd::MultiHypergraph graph(libraryInstance);

    std::string path(MFJONES_PATH);

    htd::ExternalPipeTreeDecompositionAlgorithm algorithm(libraryInstance, path + "/tw-heuristic", 0);

    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph);

    ASSERT_NE(decomposition, nullptr);

    ASSERT_EQ(decomposition->vertexCount(), (std::size_t)1);

    EXPECT_EQ(decomposition->edgeCount(), (std::size_t)0);

    ASSERT_EQ(decomposition->root(), (htd::vertex_t)1);

    const std::vector<htd::vertex_t> & bag = decomposition->bagContent(1);

    EXPECT_EQ(bag.size(), (std::size_t)0);

    delete decomposition;

    delete libraryInstance;
}

TEST(SimpleGraph_mfjones, CheckResultSimpleGraph_mfjones)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

    htd::MultiHypergraph graph(libraryInstance);

    std::string path(MFJONES_PATH);

    htd::vertex_t vertex1 = graph.addVertex();
    htd::vertex_t vertex2 = graph.addVertex();
    htd::vertex_t vertex3 = graph.addVertex();

    graph.addEdge(vertex1, vertex2);
    graph.addEdge(vertex2, vertex3);

    htd::ExternalPipeTreeDecompositionAlgorithm algorithm(libraryInstance, path + "/tw-heuristic", 0);

    algorithm.setDirectory(path);

    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph);

    ASSERT_NE(decomposition, nullptr);

    ASSERT_GE(decomposition->vertexCount(), (std::size_t)1);

    EXPECT_EQ(decomposition->edgeCount(), decomposition->vertexCount() - 1);

    ASSERT_EQ(decomposition->root(), (htd::vertex_t)1);

    ASSERT_LE(decomposition->minimumBagSize(), decomposition->maximumBagSize());

    htd::TreeDecompositionVerifier verifier;

    ASSERT_TRUE(verifier.verify(graph, *decomposition));

    delete decomposition;

    delete libraryInstance;
}

#endif
#ifdef MINFILL_MRS_PATH

TEST(minfill_mrs, Testminfill_mrs)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

    htd::MultiHypergraph graph(libraryInstance);

    std::string path(MINFILL_MRS_PATH);

    htd::ExternalPipeTreeDecompositionAlgorithm algorithm(libraryInstance, path + "/minfill_mrs", 0);

    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph);

    ASSERT_NE(decomposition, nullptr);

    ASSERT_EQ(decomposition->vertexCount(), (std::size_t)1);

    EXPECT_EQ(decomposition->edgeCount(), (std::size_t)0);

    ASSERT_EQ(decomposition->root(), (htd::vertex_t)1);

    const std::vector<htd::vertex_t> & bag = decomposition->bagContent(1);

    EXPECT_EQ(bag.size(), (std::size_t)0);

    delete decomposition;

    delete libraryInstance;
}

TEST(SimpleGraph_minfill_mrs, CheckResultSimpleGraph_minfill_mrs)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

    htd::MultiHypergraph graph(libraryInstance);

    std::string path(MINFILL_MRS_PATH);

    htd::vertex_t vertex1 = graph.addVertex();
    htd::vertex_t vertex2 = graph.addVertex();
    htd::vertex_t vertex3 = graph.addVertex();

    graph.addEdge(vertex1, vertex2);
    graph.addEdge(vertex2, vertex3);

    htd::ExternalPipeTreeDecompositionAlgorithm algorithm(libraryInstance, path + "/minfill_mrs", 0);

    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph);

    ASSERT_NE(decomposition, nullptr);

    ASSERT_GE(decomposition->vertexCount(), (std::size_t)1);

    EXPECT_EQ(decomposition->edgeCount(), decomposition->vertexCount() - 1);

    ASSERT_EQ(decomposition->root(), (htd::vertex_t)1);

    ASSERT_LE(decomposition->minimumBagSize(), decomposition->maximumBagSize());

    htd::TreeDecompositionVerifier verifier;

    ASSERT_TRUE(verifier.verify(graph, *decomposition));

    delete decomposition;

    delete libraryInstance;
}

#endif
#endif

int main(int argc, char ** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
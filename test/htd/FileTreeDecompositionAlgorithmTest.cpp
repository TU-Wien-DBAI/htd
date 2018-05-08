#include <gtest/gtest.h>

#include <htd/main.hpp>

class FileTreeDecompositionAlgorithmTest : public ::testing::Test
{
    public:
        FileTreeDecompositionAlgorithmTest(void)
        {

        }

        virtual ~FileTreeDecompositionAlgorithmTest()
        {

        }

        void SetUp()
        {

        }

        void TearDown()
        {

        }
};

std::string graph_1 = std::string("s td 1 0 0\n"
                                  "b 1 ");
std::string graph_2 = std::string("s td 3 1 3\n"
                                  "b 1 1 \n"
                                  "b 2 3 \n"
                                  "b 3 2 \n"
                                  "1 2\n"
                                  "2 3");
std::string graph_3 = std::string("s td 2 2 3\n"
                                  "b 1 2 3 \n"
                                  "b 2 1 2 \n"
                                  "1 2");

TEST(FileTreeDecompositionAlgorithmTest, CheckResultEmptyGraph)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

    htd::MultiHypergraph graph(libraryInstance);

    htd::FileTreeDecompositionAlgorithm algorithm(libraryInstance, graph_1, false);

    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph);

    ASSERT_NE(decomposition, nullptr);

    ASSERT_EQ(decomposition->vertexCount(), (std::size_t) 1);

    EXPECT_EQ(decomposition->edgeCount(), (std::size_t) 0);

    ASSERT_EQ(decomposition->root(), (htd::vertex_t) 1);

    const std::vector<htd::vertex_t> & bag = decomposition->bagContent(1);

    EXPECT_EQ(bag.size(), (std::size_t) 0);

    delete decomposition;

    delete libraryInstance;
}

TEST(FileTreeDecompositionAlgorithmTest, CheckResultDisconnectedGraph)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

    htd::MultiHypergraph graph(libraryInstance);

    graph.addVertex();
    graph.addVertex();
    graph.addVertex();

    htd::FileTreeDecompositionAlgorithm algorithm(libraryInstance, graph_2, false);

    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph);

    ASSERT_NE(decomposition, nullptr);

    ASSERT_GE(decomposition->vertexCount(), (std::size_t) 1);

    EXPECT_EQ(decomposition->edgeCount(), decomposition->vertexCount() - 1);

    ASSERT_EQ(decomposition->root(), (htd::vertex_t) 1);

    ASSERT_LE(decomposition->minimumBagSize(), decomposition->maximumBagSize());

    htd::TreeDecompositionVerifier verifier;

    ASSERT_TRUE(verifier.verify(graph, *decomposition));

    delete decomposition;

    delete libraryInstance;
}

TEST(FileTreeDecompositionAlgorithmTest, CheckResultSimpleGraph)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

    htd::MultiHypergraph graph(libraryInstance);

    htd::vertex_t vertex1 = graph.addVertex();
    htd::vertex_t vertex2 = graph.addVertex();
    htd::vertex_t vertex3 = graph.addVertex();

    graph.addEdge(vertex1, vertex2);
    graph.addEdge(vertex2, vertex3);

    htd::FileTreeDecompositionAlgorithm algorithm(libraryInstance, graph_3, false);

    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph);

    ASSERT_NE(decomposition, nullptr);

    ASSERT_GE(decomposition->vertexCount(), (std::size_t) 1);

    EXPECT_EQ(decomposition->edgeCount(), decomposition->vertexCount() - 1);

    ASSERT_EQ(decomposition->root(), (htd::vertex_t) 1);

    ASSERT_LE(decomposition->minimumBagSize(), decomposition->maximumBagSize());

    htd::TreeDecompositionVerifier verifier;

    ASSERT_TRUE(verifier.verify(graph, *decomposition));

    delete decomposition;

    delete libraryInstance;
}

class BagSizeLabelingFunction : public htd::ILabelingFunction
{
    public:
        BagSizeLabelingFunction(const htd::LibraryInstance * const manager) : managementInstance_(manager)
        {

        }

        virtual ~BagSizeLabelingFunction()
        {

        }

        std::string name() const HTD_OVERRIDE
        {
            return "BAG_SIZE";
        }

        htd::ILabel * computeLabel(const std::vector<htd::vertex_t> & vertices, const htd::ILabelCollection & labels) const HTD_OVERRIDE
        {
            HTD_UNUSED(labels)

            return new htd::Label<std::size_t>(vertices.size());
        }

        htd::ILabel * computeLabel(const htd::ConstCollection<htd::vertex_t> & vertices, const htd::ILabelCollection & labels) const HTD_OVERRIDE
        {
            HTD_UNUSED(labels)

            return new htd::Label<std::size_t>(vertices.size());
        }

        const htd::LibraryInstance * managementInstance(void) const HTD_NOEXCEPT HTD_OVERRIDE
        {
            return managementInstance_;
        }

        void setManagementInstance(const htd::LibraryInstance * const manager) HTD_OVERRIDE
        {
            HTD_ASSERT(manager != nullptr)

            managementInstance_ = manager;
        }

#ifndef HTD_USE_VISUAL_STUDIO_COMPATIBILITY_MODE

        BagSizeLabelingFunction * clone(void) const HTD_OVERRIDE
        {
            return new BagSizeLabelingFunction(managementInstance_);
        }

#else
        BagSizeLabelingFunction * clone(void) const
        {
            return new BagSizeLabelingFunction(managementInstance_);
        }

        htd::ILabelingFunction * cloneLabelingFunction(void) const HTD_OVERRIDE
        {
            return new BagSizeLabelingFunction(managementInstance_);
        }

        htd::IDecompositionManipulationOperation * cloneDecompositionManipulationOperation(void) const HTD_OVERRIDE
        {
            return new BagSizeLabelingFunction(managementInstance_);
        }
#endif

    private:
        const htd::LibraryInstance * managementInstance_;
};

TEST(FileTreeDecompositionAlgorithmTest, CheckResultSimpleGraphWithLabelingFunction)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

    htd::MultiHypergraph graph(libraryInstance);

    htd::vertex_t vertex1 = graph.addVertex();
    htd::vertex_t vertex2 = graph.addVertex();
    htd::vertex_t vertex3 = graph.addVertex();

    graph.addEdge(vertex1, vertex2);
    graph.addEdge(vertex2, vertex3);

    htd::FileTreeDecompositionAlgorithm algorithm(libraryInstance, graph_3, false);

    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph, 1, new BagSizeLabelingFunction(libraryInstance));

    ASSERT_NE(decomposition, nullptr);

    ASSERT_GE(decomposition->vertexCount(), (std::size_t) 1);

    EXPECT_EQ(decomposition->edgeCount(), decomposition->vertexCount() - 1);

    ASSERT_EQ(decomposition->root(), (htd::vertex_t) 1);

    ASSERT_LE(decomposition->minimumBagSize(), decomposition->maximumBagSize());

    htd::TreeDecompositionVerifier verifier;

    ASSERT_TRUE(verifier.verify(graph, *decomposition));

    for (htd::vertex_t vertex : decomposition->vertices())
    {
        ASSERT_EQ(decomposition->bagSize(vertex), htd::accessLabel<std::size_t>(decomposition->vertexLabel("BAG_SIZE", vertex)));
    }

    delete decomposition;

    delete libraryInstance;
}

class BagSizeLabelingFunction2 : public htd::ILabelingFunction
{
    public:
        BagSizeLabelingFunction2(const htd::LibraryInstance * const manager) : managementInstance_(manager)
        {

        }

        virtual ~BagSizeLabelingFunction2()
        {

        }

        std::string name() const HTD_OVERRIDE
        {
            return "BAG_SIZE_TIMES_2";
        }

        htd::ILabel * computeLabel(const std::vector<htd::vertex_t> & vertices, const htd::ILabelCollection & labels) const HTD_OVERRIDE
        {
            HTD_UNUSED(labels)

            return new htd::Label<std::size_t>(vertices.size() + htd::accessLabel<std::size_t>(labels.label("BAG_SIZE")));
        }

        htd::ILabel * computeLabel(const htd::ConstCollection<htd::vertex_t> & vertices, const htd::ILabelCollection & labels) const HTD_OVERRIDE
        {
            HTD_UNUSED(labels)

            return new htd::Label<std::size_t>(vertices.size() + htd::accessLabel<std::size_t>(labels.label("BAG_SIZE")));
        }

        const htd::LibraryInstance * managementInstance(void) const HTD_NOEXCEPT HTD_OVERRIDE
        {
            return managementInstance_;
        }

        void setManagementInstance(const htd::LibraryInstance * const manager) HTD_OVERRIDE
        {
            HTD_ASSERT(manager != nullptr)

            managementInstance_ = manager;
        }

#ifndef HTD_USE_VISUAL_STUDIO_COMPATIBILITY_MODE

        BagSizeLabelingFunction2 * clone(void) const HTD_OVERRIDE
        {
            return new BagSizeLabelingFunction2(managementInstance_);
        }

#else
        BagSizeLabelingFunction2 * clone(void) const
        {
            return new BagSizeLabelingFunction2(managementInstance_);
        }

        htd::ILabelingFunction * cloneLabelingFunction(void) const HTD_OVERRIDE
        {
            return new BagSizeLabelingFunction2(managementInstance_);
        }

        htd::IDecompositionManipulationOperation * cloneDecompositionManipulationOperation(void) const HTD_OVERRIDE
        {
            return new BagSizeLabelingFunction2(managementInstance_);
        }
#endif

    private:
        const htd::LibraryInstance * managementInstance_;
};

TEST(FileTreeDecompositionAlgorithmTest, CheckResultSimpleGraphWithLabelingFunctionVector)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

    htd::MultiHypergraph graph(libraryInstance);

    htd::vertex_t vertex1 = graph.addVertex();
    htd::vertex_t vertex2 = graph.addVertex();
    htd::vertex_t vertex3 = graph.addVertex();

    graph.addEdge(vertex1, vertex2);
    graph.addEdge(vertex2, vertex3);

    htd::FileTreeDecompositionAlgorithm algorithm(libraryInstance, graph_3, false);

    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph, {new BagSizeLabelingFunction(libraryInstance),
                                                                                     new BagSizeLabelingFunction2(libraryInstance)});

    ASSERT_NE(decomposition, nullptr);

    ASSERT_GE(decomposition->vertexCount(), (std::size_t) 1);

    EXPECT_EQ(decomposition->edgeCount(), decomposition->vertexCount() - 1);

    ASSERT_EQ(decomposition->root(), (htd::vertex_t) 1);

    ASSERT_LE(decomposition->minimumBagSize(), decomposition->maximumBagSize());

    htd::TreeDecompositionVerifier verifier;

    ASSERT_TRUE(verifier.verify(graph, *decomposition));

    for (htd::vertex_t vertex : decomposition->vertices())
    {
        ASSERT_EQ(decomposition->bagSize(vertex), htd::accessLabel<std::size_t>(decomposition->vertexLabel("BAG_SIZE", vertex)));
        ASSERT_EQ(decomposition->bagSize(vertex) * 2, htd::accessLabel<std::size_t>(decomposition->vertexLabel("BAG_SIZE_TIMES_2", vertex)));
    }

    delete decomposition;

    delete libraryInstance;
}

TEST(FileTreeDecompositionAlgorithmTest, CheckResultSimpleGraphWithLabelingFunctionVectorAndManipulationOperation)
{
    htd::LibraryInstance * libraryInstance = htd::createManagementInstance(htd::Id::FIRST);

    htd::MultiHypergraph graph(libraryInstance);

    htd::vertex_t vertex1 = graph.addVertex();
    htd::vertex_t vertex2 = graph.addVertex();
    htd::vertex_t vertex3 = graph.addVertex();

    graph.addEdge(vertex1, vertex2);
    graph.addEdge(vertex2, vertex3);

    htd::FileTreeDecompositionAlgorithm algorithm(libraryInstance, {new BagSizeLabelingFunction(libraryInstance),
                                                                    new htd::JoinNodeReplacementOperation(libraryInstance)}, graph_3, false);

    htd::ITreeDecomposition * decomposition = algorithm.computeDecomposition(graph, {new BagSizeLabelingFunction2(libraryInstance)});

    ASSERT_NE(decomposition, nullptr);

    ASSERT_GE(decomposition->vertexCount(), (std::size_t) 1);

    EXPECT_EQ(decomposition->edgeCount(), decomposition->vertexCount() - 1);

    ASSERT_EQ(decomposition->root(), (htd::vertex_t) 1);

    ASSERT_LE(decomposition->minimumBagSize(), decomposition->maximumBagSize());

    htd::TreeDecompositionVerifier verifier;

    ASSERT_TRUE(verifier.verify(graph, *decomposition));

    for (htd::vertex_t vertex : decomposition->vertices())
    {
        ASSERT_EQ(decomposition->bagSize(vertex), htd::accessLabel<std::size_t>(decomposition->vertexLabel("BAG_SIZE", vertex)));
        ASSERT_EQ(decomposition->bagSize(vertex) * 2, htd::accessLabel<std::size_t>(decomposition->vertexLabel("BAG_SIZE_TIMES_2", vertex)));
    }

    delete decomposition;

    delete libraryInstance;
}


int main(int argc, char ** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}

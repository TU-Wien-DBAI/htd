#ifndef HTD_FILEGRAPHDECOMPOSITIONALGORITHM_CPP
#define HTD_FILEGRAPHDECOMPOSITIONALGORITHM_CPP

#include <htd/Globals.hpp>
#include <htd/Helpers.hpp>
#include <htd/FileGraphDecompositionAlgorithm.hpp>
#include <htd/OrderingAlgorithmFactory.hpp>
#include <htd/GraphDecompositionFactory.hpp>
#include <htd/IWidthLimitableOrderingAlgorithm.hpp>
#include <htd/GraphPreprocessorFactory.hpp>

#include <cstdarg>
#include <stack>
#include <unordered_map>
#include <sstream>
#include <fstream>

/**
 *  Private implementation details of class htd::FileGraphDecompostionAlgorithm.
 */
struct htd::FileGraphDecompositionAlgorithm::Implementation
{
    /**
     *  Constructor for the implementation details structure.
     *
     *  @param[in] manager   The management instance to which the current object instance belongs.
     */
    Implementation(const htd::LibraryInstance * const manager, const std::string & decomposition, const bool & isPath) : managementInstance_(manager), orderingAlgorithm_(manager->orderingAlgorithmFactory().createInstance()), labelingFunctions_(), postProcessingOperations_(), compressionEnabled_(true), computeInducedEdges_(true)
    {
        if (isPath)
        {
            std::string inputLine;

            std::ifstream fileIn(decomposition);

            std::stringbuf treeD;

            while (getline(fileIn, inputLine))
            {
                treeD.sputn(inputLine.c_str(), inputLine.size());

                treeD.sputn("\n", 1);
            }
            this->decomposition = std::string(treeD.str());
        }
        else
        {
            this->decomposition = std::string(decomposition);
        }
    }

    virtual ~Implementation()
    {
        delete orderingAlgorithm_;

        for (auto & labelingFunction : labelingFunctions_)
        {
            delete labelingFunction;
        }

        for (auto & postProcessingOperation : postProcessingOperations_)
        {
            delete postProcessingOperation;
        }
    }

    /**
     *  The path to the file containing the graph decompostion.
     */
    std::string decomposition;

    /**
     *  The management instance to which the current object instance belongs.
     */
    const htd::LibraryInstance * managementInstance_;

    /**
     *  The ordering algorithm which shall be used to compute the vertex elimination ordering.
     */
    htd::IOrderingAlgorithm * orderingAlgorithm_;

    /**
     *  The labeling functions which are applied after the decomposition was computed.
     */
    std::vector<htd::ILabelingFunction *> labelingFunctions_;

    /**
     *  The manipuation operations which are applied after the decomposition was computed.
     */
    std::vector<htd::IGraphDecompositionManipulationOperation *> postProcessingOperations_;

    /**
     *  A boolean flag indicating whether the computed decompositions shall contain only subset-maximal bags.
     */
    bool compressionEnabled_;

    /**
     *  A boolean flag indicating whether the hyperedges induced by a respective bag shall be computed.
     */
    bool computeInducedEdges_;

    /**
     *  Compute a new mutable graph decompostion of the given graph.
     *
     *  @param[in] graph    The graph which shall be decomposed.
     *  @param[in] ordering The vertex ordering which shall be used to compute the decomposition.
     *
     *  @return A mutable graph decompostion of the given graph based on the provided vertex ordering.
     */
    htd::IMutableGraphDecomposition * computeMutableDecomposition(const htd::IMultiHypergraph & graph, const std::vector<htd::vertex_t> & ordering) const;

    /**
     *  Compute a new mutable graph decompostion of the given graph.
     *
     *  @param[in] graph                The graph which shall be decomposed.
     *  @param[in] preprocessedGraph    The input graph in preprocessed format.
     *  @param[in] maxBagSize           The upper bound for the maximum bag size of the decomposition.
     *  @param[in] maxIterationCount    The maximum number of iterations resulting in a higher maximum bag size than maxBagSize after which a null-pointer is returned.
     *
     *  @return A pair consisting of a mutable graph decompostion of the given graph or a null-pointer in case that no decomposition with a appropriate maximum bag size could be found after maxIterationCount iterations and the number of iterations actually needed to find the decomposition at hand.
     */
    std::pair<htd::IMutableGraphDecomposition *, std::size_t> computeMutableDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, std::size_t maxBagSize, std::size_t maxIterationCount) const;

    /**
     *  Get the vertex which is ranked first in the vertex elimination ordering.
     *
     *  @param[in] vertices         The set of vertices which shall be investigated.
     *  @param[in] ordering         The vertex elimination ordering.
     *  @param[in] vertexIndices    The indices of the vertices in the vertex elimination ordering.
     *
     *  @return The vertex which is ranked first in the vertex elimination ordering.
     */
    htd::vertex_t getMinimumVertex(const std::vector<htd::vertex_t> & vertices, const std::vector<htd::vertex_t> & ordering, const std::vector<htd::index_t> & vertexIndices) const;

    /**
     *  Get the vertex which is ranked first in the vertex elimination ordering.
     *
     *  @param[in] vertices         The set of vertices which shall be investigated.
     *  @param[in] ordering         The vertex elimination ordering.
     *  @param[in] vertexIndices    The indices of the vertices in the vertex elimination ordering.
     *  @param[in] excludedVertex   The vertex which shall be ignored during the determination of the first vertex in the ordering.
     *
     *  @return The vertex which is ranked first in the vertex elimination ordering.
     */
    htd::vertex_t getMinimumVertex(const std::vector<htd::vertex_t> & vertices, const std::vector<htd::vertex_t> & ordering, const std::vector<htd::index_t> & vertexIndices, htd::vertex_t excludedVertex) const;

    /**
     *  Compress the given decomposition by retaining only subset-maximal bags.
     *
     *  @param[in] startingVertex       The starting vertex.
     *  @param[in] neighbors            The neighborhood relation which shall be updated.
     *  @param[in] bagContent           The bag contents which might be swapped during the traversal.
     *  @param[in] unvisitedVertices    The set of unvisited vertices which is updated during the traversal.
     *  @param[in] relevantVertices     The set of relevant vertices which will be updated where relevance refers to subset-maximality.
     *  @param[in] inducedEdges         A vector holding the indices of the edges which are induced by the bag content associated with a vertex.
     *  @param[in] edgeTarget           A vector holding the first target node for each edge.
     */
    void compressDecomposition(htd::vertex_t startingVertex,
                               std::vector<std::vector<htd::vertex_t>> & neighbors,
                               std::vector<std::vector<htd::vertex_t>> & bagContent,
                               std::unordered_set<htd::vertex_t> & unvisitedVertices,
                               std::vector<htd::vertex_t> & relevantVertices,
                               std::vector<std::vector<htd::index_t>> & inducedEdges,
                               std::vector<htd::index_t> & edgeTarget) const;

    /**
     *  Compress the given decomposition by retaining only subset-maximal bags.
     *
     *  If the bag of the vertex is a subset of the bag of its parent, the vertex is removed. If the bag of
     *  the vertex is a superset of the bag of its parent, the bag contents are swapped and the vertex is
     *  removed. Otherwise, the decomposition is left unchanged.
     *
     *  @param[in] vertex           The vertex.
     *  @param[in] parent           The parent of the vertex during the traversal.
     *  @param[in] neighbors        The neighborhood relation which shall be updated.
     *  @param[in] bagContent       The bag contents which might be swapped during the traversal.
     *  @param[in] relevantVertices The set of relevant vertices which will be updated where relevance refers to subset-maximality.
     *  @param[in] inducedEdges     A vector holding the indices of the edges which are induced by the bag content associated with a vertex.
     *  @param[in] edgeTarget       A vector holding the first target node for each edge.
     */
    void compressDecomposition(htd::vertex_t vertex, htd::vertex_t parent,
                               std::vector<std::vector<htd::vertex_t>> & neighbors,
                               std::vector<std::vector<htd::vertex_t>> & bagContent,
                               std::vector<htd::vertex_t> & relevantVertices,
                               std::vector<std::vector<htd::index_t>> & inducedEdges,
                               std::vector<htd::index_t> & edgeTarget) const;

    /**
     *  Update the given decomposition by performing pre-order traversal.
     *
     *  @param[in] graph                    The graph from which the decomposition was computed.
     *  @param[in] decomposition            The decomposition which shall be updated.
     *  @param[in] startingVertex           The starting vertex.
     *  @param[in] neighbors                The neighborhood relation which shall be used.
     *  @param[in] bagContent               The bag contents.
     *  @param[in] unvisitedVertices        The set of unvisited vertices which is updated during the traversal.
     *  @param[in] inducedEdges             A vector holding the indices of the edges which are induced by the bag content associated with a vertex.
     *  @param[in] decompositionVertices    A mapping between the vertices and their counterparts in the decomposition.
     *
     *  @note The bag contents and the induced edges are moved into the decomposition during this operation.
     */
    void updateDecomposition(const htd::IMultiHypergraph & graph,
                             htd::IMutableGraphDecomposition & decomposition,
                             htd::vertex_t startingVertex,
                             const std::vector<std::vector<htd::vertex_t>> & neighbors,
                             std::vector<std::vector<htd::vertex_t>> & bagContent,
                             std::vector<std::vector<htd::index_t>> & inducedEdges,
                             std::unordered_set<htd::vertex_t> & unvisitedVertices,
                             std::unordered_map<htd::vertex_t, htd::vertex_t> & decompositionVertices) const;

    /**
     *  Check whether two sets are subset-maximal with respect to the other set.
     *
     *  @param[in] set1 The first set.
     *  @param[in] set2 The second set.
     *
     *  @return This function returns -1 if the first set is a superset of or identical to the second set.
     *  If the second set is a proper superset of the first set, the return value is 1. Otherwise, the
     *  return value is 0.
     */
    int is_maximal(const std::vector<htd::vertex_t> & set1, const std::vector<htd::vertex_t> & set2) const;

    /**
     *  Distribute a given edge, identified by its index, in the decomposition so that the information about induced edges is updated.
     *
     *  @param[in] edgeIndex        The index of the edge which shall be distributed.
     *  @param[in] edge             The sorted elements of the edge which shall be distributed.
     *  @param[in] startBucket      The identifier of the node from which the process shall start.
     *  @param[in] buckets          The available buckets.
     *  @param[in] neighbors        The neighbors of the buckets.
     *  @param[in] inducedEdges     The set of edge indices induced by a bucket.
     *  @param[in] lastAssignedEdge The identifier of the last edge which was assigned to a bucket.
     *  @param[in] originStack      The stack instance used for backtracking.
     */
    void distributeEdge(htd::index_t edgeIndex,
                        const std::vector<htd::vertex_t> & edge,
                        htd::vertex_t startBucket,
                        const std::vector<std::vector<htd::vertex_t>> & buckets,
                        const std::vector<std::vector<htd::vertex_t>> & neighbors,
                        std::vector<std::vector<htd::index_t>> & inducedEdges,
                        std::vector<htd::id_t> & lastAssignedEdge,
                        std::stack<htd::vertex_t> & originStack) const;

    /**
     *  Distribute a given edge, identified by its index, in the decomposition so that the information about induced edges is updated.
     *
     *  @param[in] edgeIndex        The index of the edge which shall be distributed.
     *  @param[in] vertex1          The first vertex (i.e., the one with lower ID) of the edge which shall be distributed.
     *  @param[in] vertex2          The second vertex (i.e., the one with higher ID) of the edge which shall be distributed.
     *  @param[in] startBucket      The identifier of the node from which the process shall start.
     *  @param[in] buckets          The available buckets.
     *  @param[in] neighbors        The neighbors of the buckets.
     *  @param[in] inducedEdges     The set of edge indices induced by a bucket.
     *  @param[in] lastAssignedEdge The identifier of the last edge which was assigned to a bucket.
     *  @param[in] originStack      The stack instance used for backtracking.
     */
    void distributeEdge(htd::index_t edgeIndex,
                        htd::vertex_t vertex1,
                        htd::vertex_t vertex2,
                        htd::vertex_t startBucket,
                        const std::vector<std::vector<htd::vertex_t>> & buckets,
                        const std::vector<std::vector<htd::vertex_t>> & neighbors,
                        std::vector<std::vector<htd::index_t>> & inducedEdges,
                        std::vector<htd::id_t> & lastAssignedEdge,
                        std::stack<htd::vertex_t> & originStack) const;

    /**
     *  Compute the set union of two sets and store the result in the first set.
     *
     *  @param[in,out] set1         The first set
     *  @param[in] set2             The second set.
     *  @param[in] ignoredVertex    The vertex which shall be ignored if it occurs in the second set.
     */
    void set_union(std::vector<htd::vertex_t> & set1,
                   const std::vector<htd::vertex_t> & set2,
                   htd::vertex_t ignoredVertex) const
    {
        std::vector<htd::vertex_t> tmp;
        tmp.reserve(set2.size());

        auto first1 = set1.begin();
        auto first2 = set2.begin();

        auto last1 = set1.end();
        auto last2 = set2.end();

        while (first1 != last1 && first2 != last2)
        {
            if (*first1 < *first2)
            {
                ++first1;
            }
            else if (*first2 < *first1)
            {
                if (*first2 != ignoredVertex)
                {
                    tmp.push_back(*first2);
                }

                ++first2;
            }
            else
            {
                ++first1;

                //Skip common value in set 2.
                ++first2;
            }
        }

        std::copy_if(first2, last2, std::back_inserter(tmp), [&](const htd::vertex_t vertex) { return vertex != ignoredVertex; });

        if (!tmp.empty())
        {
            htd::inplace_merge(set1, tmp);
        }
    }

    void parseEdgeLine(std::string & item, IMutableGraphDecomposition * ret) const;

    void parseBagLine(std::string item, IMutableGraphDecomposition * ret, const IMultiHypergraph & graph) const;

    void getInducedEdges(std::vector<vertex_t> & bucket, const IMultiHypergraph & graph, std::vector<index_t> & inducedEdges) const;
};

htd::FileGraphDecompositionAlgorithm::FileGraphDecompositionAlgorithm(const htd::LibraryInstance * const manager, const std::string & decompPath, const bool & isPath) : implementation_(new Implementation(manager, decompPath, isPath))
{

}

htd::FileGraphDecompositionAlgorithm::FileGraphDecompositionAlgorithm(const htd::LibraryInstance * const manager, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations, const std::string & decompPath, const bool & isPath) : implementation_(new Implementation(manager, decompPath, isPath))
{
    setManipulationOperations(manipulationOperations);
}

htd::FileGraphDecompositionAlgorithm::~FileGraphDecompositionAlgorithm()
{

}

htd::IGraphDecomposition * htd::FileGraphDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph) const
{
    return computeDecomposition(graph, std::vector<htd::IDecompositionManipulationOperation *>());
}

htd::IGraphDecomposition * htd::FileGraphDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations) const
{
    return computeDecomposition(graph, manipulationOperations, (std::size_t) -1, 1).first;
}

std::pair<htd::IGraphDecomposition *, std::size_t> htd::FileGraphDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, std::size_t maxBagSize, std::size_t maxIterationCount) const
{
    return computeDecomposition(graph, std::vector<htd::IDecompositionManipulationOperation *>(), maxBagSize, maxIterationCount);
}

std::pair<htd::IGraphDecomposition *, std::size_t> htd::FileGraphDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations, std::size_t maxBagSize, std::size_t maxIterationCount) const
{
    htd::IGraphPreprocessor * preprocessor = implementation_->managementInstance_->graphPreprocessorFactory().createInstance();

    htd::IPreprocessedGraph * preprocessedGraph = preprocessor->prepare(graph);

    std::pair<htd::IGraphDecomposition *, std::size_t> ret =
        computeDecomposition(graph, *preprocessedGraph, manipulationOperations, maxBagSize, maxIterationCount);

    delete preprocessedGraph;
    delete preprocessor;

    return ret;
}

htd::IGraphDecomposition * htd::FileGraphDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph) const
{
    return computeDecomposition(graph, preprocessedGraph, std::vector<htd::IDecompositionManipulationOperation *>());
}

htd::IGraphDecomposition * htd::FileGraphDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations) const
{
    return computeDecomposition(graph, preprocessedGraph, manipulationOperations, (std::size_t) -1, 1).first;
}

std::pair<htd::IGraphDecomposition *, std::size_t> htd::FileGraphDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, std::size_t maxBagSize, std::size_t maxIterationCount) const
{
    return computeDecomposition(graph, preprocessedGraph, std::vector<htd::IDecompositionManipulationOperation *>(), maxBagSize, maxIterationCount);
}

std::pair<htd::IGraphDecomposition *, std::size_t> htd::FileGraphDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations, std::size_t maxBagSize, std::size_t maxIterationCount) const
{
    std::pair<htd::IMutableGraphDecomposition *, std::size_t> ret = implementation_->computeMutableDecomposition(graph, preprocessedGraph, maxBagSize, maxIterationCount);

    htd::IMutableGraphDecomposition * decomposition = ret.first;

    if (decomposition != nullptr)
    {
        std::vector<htd::ILabelingFunction *> labelingFunctions;

        std::vector<htd::IGraphDecompositionManipulationOperation *> postProcessingOperations;

        for (htd::IDecompositionManipulationOperation * operation : manipulationOperations)
        {
            htd::ILabelingFunction * labelingFunction = dynamic_cast<htd::ILabelingFunction *>(operation);

            if (labelingFunction != nullptr)
            {
                labelingFunctions.push_back(labelingFunction);
            }

            htd::IGraphDecompositionManipulationOperation * manipulationOperation = dynamic_cast<htd::IGraphDecompositionManipulationOperation *>(operation);

            if (manipulationOperation != nullptr)
            {
                postProcessingOperations.push_back(manipulationOperation);
            }
        }

        for (const auto & operation : implementation_->postProcessingOperations_)
        {
            operation->apply(graph, *decomposition);
        }

        for (const auto & operation : postProcessingOperations)
        {
            operation->apply(graph, *decomposition);
        }

        for (const auto & labelingFunction : implementation_->labelingFunctions_)
        {
            for (htd::vertex_t vertex : decomposition->vertices())
            {
                htd::ILabelCollection * labelCollection = decomposition->labelings().exportVertexLabelCollection(vertex);

                htd::ILabel * newLabel = labelingFunction->computeLabel(decomposition->bagContent(vertex), *labelCollection);

                delete labelCollection;

                decomposition->setVertexLabel(labelingFunction->name(), vertex, newLabel);
            }
        }

        for (const auto & labelingFunction : labelingFunctions)
        {
            for (htd::vertex_t vertex : decomposition->vertices())
            {
                htd::ILabelCollection * labelCollection = decomposition->labelings().exportVertexLabelCollection(vertex);

                htd::ILabel * newLabel = labelingFunction->computeLabel(decomposition->bagContent(vertex), *labelCollection);

                delete labelCollection;

                decomposition->setVertexLabel(labelingFunction->name(), vertex, newLabel);
            }
        }

        for (auto & labelingFunction : labelingFunctions)
        {
            delete labelingFunction;
        }

        for (auto & postProcessingOperation : postProcessingOperations)
        {
            delete postProcessingOperation;
        }
    }
    else
    {
        for (htd::IDecompositionManipulationOperation * operation : manipulationOperations)
        {
            delete operation;
        }
    }

    return ret;
}

htd::IGraphDecomposition * htd::FileGraphDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, int manipulationOperationCount, ...) const
{
    va_list arguments;

    va_start(arguments, manipulationOperationCount);

    std::vector<htd::IDecompositionManipulationOperation *> manipulationOperations;

    for (int manipulationOperationIndex = 0; manipulationOperationIndex < manipulationOperationCount; manipulationOperationIndex++)
    {
        manipulationOperations.push_back(va_arg(arguments, htd::IDecompositionManipulationOperation *));
    }

    va_end(arguments);

    return computeDecomposition(graph, manipulationOperations);
}

htd::IGraphDecomposition * htd::FileGraphDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, int manipulationOperationCount, ...) const
{
    va_list arguments;

    va_start(arguments, manipulationOperationCount);

    std::vector<htd::IDecompositionManipulationOperation *> manipulationOperations;

    for (int manipulationOperationIndex = 0; manipulationOperationIndex < manipulationOperationCount; manipulationOperationIndex++)
    {
        manipulationOperations.push_back(va_arg(arguments, htd::IDecompositionManipulationOperation *));
    }

    va_end(arguments);

    return computeDecomposition(graph, preprocessedGraph, manipulationOperations);
}

void htd::FileGraphDecompositionAlgorithm::setOrderingAlgorithm(htd::IOrderingAlgorithm * algorithm)
{
    HTD_ASSERT(algorithm != nullptr)

    delete implementation_->orderingAlgorithm_;

    implementation_->orderingAlgorithm_ = algorithm;
}

void htd::FileGraphDecompositionAlgorithm::setManipulationOperations(const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations)
{
    for (auto & labelingFunction : implementation_->labelingFunctions_)
    {
        delete labelingFunction;
    }

    for (auto & postProcessingOperation : implementation_->postProcessingOperations_)
    {
        delete postProcessingOperation;
    }

    implementation_->labelingFunctions_.clear();

    implementation_->postProcessingOperations_.clear();

    addManipulationOperations(manipulationOperations);
}

void htd::FileGraphDecompositionAlgorithm::addManipulationOperation(htd::IDecompositionManipulationOperation * manipulationOperation)
{
    bool assigned = false;

    htd::ILabelingFunction * labelingFunction = dynamic_cast<htd::ILabelingFunction *>(manipulationOperation);

    if (labelingFunction != nullptr)
    {
        implementation_->labelingFunctions_.emplace_back(labelingFunction);

        assigned = true;
    }

    htd::IGraphDecompositionManipulationOperation * newManipulationOperation = dynamic_cast<htd::IGraphDecompositionManipulationOperation *>(manipulationOperation);

    if (newManipulationOperation != nullptr)
    {
        implementation_->postProcessingOperations_.emplace_back(newManipulationOperation);

        assigned = true;
    }

    if (!assigned)
    {
        delete manipulationOperation;
    }
}

void htd::FileGraphDecompositionAlgorithm::addManipulationOperations(const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations)
{
    for (htd::IDecompositionManipulationOperation * operation : manipulationOperations)
    {
        addManipulationOperation(operation);
    }
}

bool htd::FileGraphDecompositionAlgorithm::isSafelyInterruptible(void) const
{
    return false;
}

const htd::LibraryInstance * htd::FileGraphDecompositionAlgorithm::managementInstance(void) const HTD_NOEXCEPT
{
    return implementation_->managementInstance_;
}

void htd::FileGraphDecompositionAlgorithm::setManagementInstance(const htd::LibraryInstance * const manager)
{
    HTD_ASSERT(manager != nullptr)

    implementation_->managementInstance_ = manager;
}

bool htd::FileGraphDecompositionAlgorithm::isCompressionEnabled(void) const
{
    return implementation_->compressionEnabled_;
}

void htd::FileGraphDecompositionAlgorithm::setCompressionEnabled(bool compressionEnabled)
{
    implementation_->compressionEnabled_ = compressionEnabled;
}

bool htd::FileGraphDecompositionAlgorithm::isComputeInducedEdgesEnabled(void) const
{
    return implementation_->computeInducedEdges_;
}

void htd::FileGraphDecompositionAlgorithm::setComputeInducedEdgesEnabled(bool computeInducedEdgesEnabled)
{
    implementation_->computeInducedEdges_ = computeInducedEdgesEnabled;
}

htd::FileGraphDecompositionAlgorithm * htd::FileGraphDecompositionAlgorithm::clone(void) const
{
    htd::FileGraphDecompositionAlgorithm * ret = new htd::FileGraphDecompositionAlgorithm(implementation_->managementInstance_, implementation_->decomposition, false);

    ret->setCompressionEnabled(implementation_->compressionEnabled_);
    ret->setComputeInducedEdgesEnabled(implementation_->computeInducedEdges_);

    for (const auto & labelingFunction : implementation_->labelingFunctions_)
    {
#ifndef HTD_USE_VISUAL_STUDIO_COMPATIBILITY_MODE
        ret->addManipulationOperation(labelingFunction->clone());
#else
        ret->addManipulationOperation(labelingFunction->cloneLabelingFunction());
#endif
    }

    for (const auto & postProcessingOperation : implementation_->postProcessingOperations_)
    {
#ifndef HTD_USE_VISUAL_STUDIO_COMPATIBILITY_MODE
        ret->addManipulationOperation(postProcessingOperation->clone());
#else
        ret->addManipulationOperation(postProcessingOperation->cloneGraphDecompositionManipulationOperation());
#endif
    }

#ifndef HTD_USE_VISUAL_STUDIO_COMPATIBILITY_MODE
    ret->setOrderingAlgorithm(implementation_->orderingAlgorithm_->clone());
#else
    ret->setOrderingAlgorithm(implementation_->orderingAlgorithm_->cloneOrderingAlgorithm());
#endif

    return ret;
}

std::pair<htd::IMutableGraphDecomposition *, std::size_t> htd::FileGraphDecompositionAlgorithm::Implementation::computeMutableDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, std::size_t maxBagSize, std::size_t maxIterationCount) const
{
    const htd::LibraryInstance & managementInstance = *managementInstance_;

    htd::IMutableGraphDecomposition * ret = managementInstance.graphDecompositionFactory().createInstance();

    std::size_t iterations = 0;

    std::size_t size = graph.vertexCount();

    if (size > 0)
    {
        htd::vertex_t lastVertex = graph.vertexAtPosition(size - 1);

        std::vector<std::vector<htd::vertex_t>> buckets(lastVertex + 1);

        std::vector<std::vector<htd::index_t>> inducedEdges(lastVertex + 1);

        std::stringstream ss(decomposition);

        std::string item;

        while (getline(ss, item))
        {
            //ignore empty line
            if (item.length() > 0)
            {
                char type = item.at(0);
                if (type == 'c')
                {
                    //comment line (ignore)
                }
                else if (type == 's')
                {
                    //start line (ignore)
                }
                else if (type == 'b')
                {
                    //bag line
                    parseBagLine(item, ret, graph);
                }
                else
                {
                    //edge line
                    parseEdgeLine(item, ret);
                }
            }
        }

        ConstIterator<vertex_t> iter;

        // check if all vertices are in a bag
        for (iter = graph.vertices().begin(); iter != graph.vertices().end(); ++iter)
        {
            bool found = false;

            for (ConstIterator<vertex_t> iter_decomp = ret->vertices().begin(); iter_decomp != ret->vertices().end(); ++iter_decomp)
            {
                std::vector<vertex_t> & bagContent = ret->mutableBagContent(*iter_decomp);

                if (std::find(bagContent.begin(), bagContent.end(), *iter) != bagContent.end())
                {
                    found = true;

                    break;
                }
            }
            if (!found)
            {
                delete ret;

                ret = nullptr;

                return std::make_pair(ret, iterations);
            }
        }

        std::size_t edgeCount = graph.edgeCount();

        const htd::ConstCollection<htd::Hyperedge> & hyperedges = graph.hyperedges();

        auto hyperedgePosition = hyperedges.begin();

        // check for every edge if it is in a bag
        for (htd::index_t index = 0; index < edgeCount; ++index)
        {
            const std::vector<htd::vertex_t> & elements = hyperedgePosition->sortedElements();

            bool found = false;

            for (ConstIterator<vertex_t> iter_decomp = ret->vertices().begin(); iter_decomp != ret->vertices().end() & !found; ++iter_decomp)
            {
                std::vector<vertex_t> & bagContent = ret->mutableBagContent(*iter_decomp);

                for (unsigned long i = 0; i <= elements.size(); i++)
                {
                    if (i == elements.size())
                    {
                        found = true;

                        break;
                    }

                    if (std::find(bagContent.begin(), bagContent.end(), elements[i]) == bagContent.end())
                    {
                        break;
                    }
                }
            }
            if (!found)
            {
                delete ret;

                ret = nullptr;

                return std::make_pair(ret, iterations);
            }
            ++hyperedgePosition;
        }
    }
    else
    {
        if (!managementInstance.isTerminated())
        {
            ret->addVertex();
        }
    }
    return std::make_pair(ret, iterations);
}

void htd::FileGraphDecompositionAlgorithm::Implementation::parseBagLine(std::string item, IMutableGraphDecomposition * ret, const IMultiHypergraph & graph) const
{
    std::vector<index_t> inducedEdges;

    std::vector<vertex_t> buckets;

    std::stringstream sline(item);

    std::string i;

    std::getline(sline, i, ' '); //b

    std::getline(sline, i, ' '); //bag number

    while (getline(sline, i, ' ')) //vertices
    {
        if (i[0] != '\r')
        {
            buckets.push_back(stoul(i));
        }
    }

    getInducedEdges(buckets, graph, inducedEdges);

    ret->addVertex(std::vector<htd::vertex_t>(buckets), graph.hyperedgesAtPositions(inducedEdges));
}

void htd::FileGraphDecompositionAlgorithm::Implementation::getInducedEdges(std::vector<vertex_t> & bucket, const IMultiHypergraph & graph, std::vector<index_t> & inducedEdges) const
{

    std::size_t edgeCount = graph.edgeCount();

    const htd::ConstCollection<htd::Hyperedge> & hyperedges = graph.hyperedges();

    auto hyperedgePosition = hyperedges.begin();

    for (htd::index_t index = 0; index < edgeCount; ++index)
    {
        const std::vector<htd::vertex_t> & elements = hyperedgePosition->sortedElements();

        switch (elements.size())
        {
            case 1:
            {
                htd::vertex_t vertex = elements[0];

                if (std::find(bucket.begin(), bucket.end(), vertex) != bucket.end())
                {
                    inducedEdges.push_back(index);
                }

                break;
            }
            case 2:
            {
                htd::vertex_t vertex1 = elements[0];

                htd::vertex_t vertex2 = elements[1];

                if (std::find(bucket.begin(), bucket.end(), vertex1) != bucket.end() && std::find(bucket.begin(), bucket.end(), vertex2) != bucket.end())
                {
                    inducedEdges.push_back(index);
                }

                break;
            }
            default:
            {
                for (unsigned long i = 0; i <= elements.size(); i++)
                {
                    if (i == elements.size())
                    {
                        inducedEdges.push_back(index);

                        break;
                    }
                    if (std::find(bucket.begin(), bucket.end(), elements[i]) == bucket.end())
                    {
                        break;
                    }
                }

                break;
            }
        }

        ++hyperedgePosition;
    }
}


void htd::FileGraphDecompositionAlgorithm::Implementation::parseEdgeLine(std::string & item, IMutableGraphDecomposition * ret) const
{
    std::stringstream sline(item);

    std::string i;

    std::vector<vertex_t> edgeNodes;

    while (getline(sline, i, ' ')) //vertices
    {
        if (i[0] != '\r')
        {
            edgeNodes.push_back(stoul(i));
        }
    }
    ret->addEdge(edgeNodes);
}

#endif /* HTD_FILEGRAPHDECOMPOSITIONALGORITHM_CPP */

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
     *  @param[in] manager          The management instance to which the current object instance belongs.
     *  @param[in] decompostion     String containing the tree decomposition or the path to the file containing the tree decomposition.
     */
    Implementation(const htd::LibraryInstance * const manager, const std::string & decomposition) : managementInstance_(manager), labelingFunctions_(), postProcessingOperations_(), decomposition_(decomposition), computeInducedEdges_(true)
    {
    }

    virtual ~Implementation()
    {
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
     *  The management instance to which the current object instance belongs.
     */
    const htd::LibraryInstance * managementInstance_;

    /**
     *  The labeling functions which are applied after the decomposition was computed.
     */
    std::vector<htd::ILabelingFunction *> labelingFunctions_;

    /**
     *  The manipuation operations which are applied after the decomposition was computed.
     */
    std::vector<htd::IGraphDecompositionManipulationOperation *> postProcessingOperations_;

    /**
     *  The string containing the decomposition in td format.
     */
    std::string decomposition_;

    /**
     *  A boolean flag indicating whether the hyperedges induced by a respective bag shall be computed.
     */
    bool computeInducedEdges_;

    /**
     *  Compute a new mutable graph decompostion of the given graph.
     *
     *  @param[in] graph                The graph which shall be decomposed.
     *  @param[in] preprocessedGraph    The input graph in preprocessed format.
     *  @param[in] maxBagSize           The upper bound for the maximum bag size of the decomposition.
     *
     *  @return A pair consisting of a mutable graph decompostion of the given graph or a null-pointer in case that the decomposition does not have a appropriate maximum bag size or the decomposition is not a valid decomposition of the graph.
     */
    htd::IMutableGraphDecomposition * computeMutableDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, std::size_t maxBagSize) const;

    /**
     * Parses the given line as edge line of the decomposition.
     *
     * @param[in] line          the edge line
     * @param[in,out] decomp    the decomposition
     */
    void parseEdgeLine(std::string & line, IMutableGraphDecomposition * decomp) const;

    /**
     * Parses the given line as bag line of the decomposition.
     *
     * @param[in] line              the bag line
     * @param[in,out] decomp        the decomposition
     * @param[in] graph             the base graph of the decomposition
     * @param[in] unfoundEdges     the base graph of the decomposition
     */
    void parseBagLine(std::string line, IMutableGraphDecomposition * decomp, const IMultiHypergraph & graph, std::unordered_set<std::vector<vertex_t>> & unfoundEdges, std::unordered_set<vertex_t> & unfoundVertices) const;

    /**
     * Computes the edges induced by the given bag.
     *
     * @param[in] bag                   the bag of the decomposition
     * @param[in] graph                 the base graph of the decomposition
     * @param[in,out] inducedEdges      the edges induced by the bag
     * @param[in,out] notfoundEdges     the edges induced by the bag
     */
    void getInducedEdges(std::vector<vertex_t> & bag, const IMultiHypergraph & graph, std::vector<index_t> & inducedEdges, std::unordered_set<std::vector<vertex_t>> & notfoundEdges) const;
};

htd::FileGraphDecompositionAlgorithm::FileGraphDecompositionAlgorithm(const htd::LibraryInstance * const manager, const std::string & decomposition) : implementation_(new Implementation(manager, decomposition))
{

}

htd::FileGraphDecompositionAlgorithm::FileGraphDecompositionAlgorithm(const htd::LibraryInstance * const manager, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations, const std::string & decomposition) : implementation_(new Implementation(manager, decomposition))
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
    return computeDecomposition(graph, manipulationOperations, (std::size_t)-1);
}

htd::IGraphDecomposition * htd::FileGraphDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, std::size_t maxBagSize) const
{
    return computeDecomposition(graph, std::vector<htd::IDecompositionManipulationOperation *>(), maxBagSize);
}

htd::IGraphDecomposition * htd::FileGraphDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations, std::size_t maxBagSize) const
{
    htd::IGraphPreprocessor * preprocessor = implementation_->managementInstance_->graphPreprocessorFactory().createInstance();

    htd::IPreprocessedGraph * preprocessedGraph = preprocessor->prepare(graph);

    htd::IGraphDecomposition * ret = computeDecomposition(graph, *preprocessedGraph, manipulationOperations, maxBagSize);

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
    return computeDecomposition(graph, preprocessedGraph, manipulationOperations, (std::size_t)-1);
}

htd::IGraphDecomposition * htd::FileGraphDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, std::size_t maxBagSize) const
{
    return computeDecomposition(graph, preprocessedGraph, std::vector<htd::IDecompositionManipulationOperation *>(), maxBagSize);
}

htd::IGraphDecomposition * htd::FileGraphDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations, std::size_t maxBagSize) const
{
    htd::IMutableGraphDecomposition * ret = implementation_->computeMutableDecomposition(graph, preprocessedGraph, maxBagSize);

    htd::IMutableGraphDecomposition * decomposition = ret;

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
    htd::FileGraphDecompositionAlgorithm * ret = new htd::FileGraphDecompositionAlgorithm(implementation_->managementInstance_, implementation_->decomposition_);

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

    return ret;
}

htd::IMutableGraphDecomposition * htd::FileGraphDecompositionAlgorithm::Implementation::computeMutableDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph &, std::size_t maxBagSize) const
{
    const htd::LibraryInstance & managementInstance = *managementInstance_;

    htd::IMutableGraphDecomposition * decomposition = managementInstance.graphDecompositionFactory().createInstance();

    std::size_t size = graph.vertexCount();

    if (size > 0)
    {
        htd::vertex_t lastVertex = graph.vertexAtPosition(size - 1);

        std::vector<std::vector<htd::vertex_t>> buckets(lastVertex + 1);

        std::vector<std::vector<htd::index_t>> inducedEdges(lastVertex + 1);

        std::string decompositionString;

        std::ifstream test(decomposition_);
        if (!test)
        {
            decompositionString = std::string(decomposition_);
        }
        else
        {
            test.close();

            std::string inputLine;

            std::ifstream fileIn(decomposition_);

            std::stringbuf treeD;

            while (getline(fileIn, inputLine))
            {
                treeD.sputn(inputLine.c_str(), inputLine.size());

                treeD.sputn("\n", 1);
            }
            decompositionString = std::string(treeD.str());
        }

        std::stringstream ss(decompositionString);

        std::string line;

        std::unordered_set<std::vector<vertex_t>> unfoundEdges;

        std::unordered_set<vertex_t> unfoundVertices;

        for (const htd::Hyperedge & edge: graph.hyperedges())
        {
            unfoundEdges.insert(edge.sortedElements());
        }

        for (const htd::vertex_t & vertex: graph.vertices())
        {
            unfoundVertices.insert(vertex);
        }

        while (getline(ss, line))
        {
            //ignore empty line
            if (line.length() > 1)
            {
                char type = line.at(0);
                switch (type)
                {
                    case 'c': // comment line - nothing to do
                        break;

                    case 's': // start line - nothing to do
                        break;

                    case 'b': // bag line
                        parseBagLine(line, decomposition, graph, unfoundEdges, unfoundVertices);
                        break;

                    default: // edge line
                        parseEdgeLine(line, decomposition);
                        break;
                }
            }
        }

        // vertices missing
        if (!unfoundVertices.empty())
        {
            delete decomposition;

            decomposition = nullptr;

            return decomposition;
        }

        // check edges
        if (!unfoundEdges.empty())
        {
            delete decomposition;

            decomposition = nullptr;

            return decomposition;
        };

        if (decomposition->maximumBagSize() > maxBagSize)
        {
            delete decomposition;

            decomposition = nullptr;

            return decomposition;
        }
    }
    else
    {
        if (!managementInstance.isTerminated())
        {
            decomposition->addVertex();
        }
    }
    return decomposition;
}

void htd::FileGraphDecompositionAlgorithm::Implementation::parseBagLine(std::string line, IMutableGraphDecomposition * decomp, const IMultiHypergraph & graph, std::unordered_set<std::vector<vertex_t>> & unfoundEdges, std::unordered_set<vertex_t> & unfoundVertices) const
{
    std::vector<index_t> inducedEdges;

    std::vector<vertex_t> buckets;

    std::stringstream sline(line);

    std::string i;

    std::getline(sline, i, ' '); //b

    std::getline(sline, i, ' '); //bag number

    while (getline(sline, i, ' ')) //vertices
    {
        if (i[0] != '\r')
        {
            buckets.push_back(stoul(i));
            unfoundVertices.erase(stoul(i));
        }
    }

    std::sort(buckets.begin(), buckets.end());

    getInducedEdges(buckets, graph, inducedEdges, unfoundEdges);

    decomp->addVertex(std::vector<htd::vertex_t>(buckets), graph.hyperedgesAtPositions(inducedEdges));
}

void htd::FileGraphDecompositionAlgorithm::Implementation::getInducedEdges(std::vector<vertex_t> & bag, const IMultiHypergraph & graph, std::vector<index_t> & inducedEdges, std::unordered_set<std::vector<vertex_t>> & notfoundEdges) const
{

    const std::unordered_set<std::vector<vertex_t>> notfoundEdgesCopy = notfoundEdges;
    for (const std::vector<vertex_t> & edge:notfoundEdgesCopy)
    {
        if (std::includes(bag.begin(), bag.end(), edge.begin(), edge.end()))
        {
            notfoundEdges.erase(edge);
        }
    }

    if (computeInducedEdges_)
    {
        for (std::size_t i = 0; i < graph.hyperedges().size(); i++)
        {
            const htd::Hyperedge & edge = graph.hyperedges()[i];
            if (std::includes(bag.begin(), bag.end(), edge.sortedElements().begin(), edge.sortedElements().end()))
            {
                inducedEdges.push_back(i);
            }
        }
    }
}


void htd::FileGraphDecompositionAlgorithm::Implementation::parseEdgeLine(std::string & line, IMutableGraphDecomposition * decomp) const
{
    std::stringstream sline(line);

    std::string i;

    std::vector<vertex_t> edgeNodes;

    while (getline(sline, i, ' ')) //vertices
    {
        if (i[0] != '\r')
        {
            edgeNodes.push_back(stoul(i));
        }
    }
    decomp->addEdge(edgeNodes);
}

#endif /* HTD_FILEGRAPHDECOMPOSITIONALGORITHM_CPP */

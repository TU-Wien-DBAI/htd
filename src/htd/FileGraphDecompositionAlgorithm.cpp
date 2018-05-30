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
    std::pair<htd::IMutableGraphDecomposition *, std::size_t> computeMutableDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, std::size_t maxBagSize, std::size_t) const;

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

        std::copy_if(first2, last2, std::back_inserter(tmp), [&](const htd::vertex_t vertex){ return vertex != ignoredVertex; });

        if (!tmp.empty())
        {
            htd::inplace_merge(set1, tmp);
        }
    }

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
     * @param[in] line          the bag line
     * @param[in,out] decomp    the decomposition
     * @param[in] graph         the base graph of the decomposition
     * @param[in] notfoundEdges         the base graph of the decomposition
     */
    void parseBagLine(std::string line, IMutableGraphDecomposition * decomp, const IMultiHypergraph & graph, std::unordered_set<std::vector<vertex_t>> & notfoundEdges) const;

    /**
     * Computes the edges induced by the given bag.
     *
     * @param[in] bag               the bag of the decomposition
     * @param[in] graph             the base graph of the decomposition
     * @param[in,out] inducedEdges  the edges induced by the bag
     * @param[in,out] notfoundEdges  the edges induced by the bag
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
    return computeDecomposition(graph, manipulationOperations, (std::size_t)-1, 1).first;
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
    return computeDecomposition(graph, preprocessedGraph, manipulationOperations, (std::size_t)-1, 1).first;
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

std::pair<htd::IMutableGraphDecomposition *, std::size_t> htd::FileGraphDecompositionAlgorithm::Implementation::computeMutableDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, std::size_t maxBagSize, std::size_t) const
{
    const htd::LibraryInstance & managementInstance = *managementInstance_;

    htd::IMutableGraphDecomposition * decomposition = managementInstance.graphDecompositionFactory().createInstance();

    std::size_t iterations = 0;

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

        std::size_t edgeCount = graph.hyperedges().size();

        const htd::ConstCollection<htd::Hyperedge> & hyperedges = graph.hyperedges();

        auto hyperedgePosition = hyperedges.begin();

        std::unordered_set<std::vector<vertex_t>> notfoundEdges;
        for (htd::index_t index = 0; index < edgeCount; ++index)
        {
            const std::vector<htd::vertex_t> & elements = hyperedgePosition->sortedElements();
            notfoundEdges.insert(elements);
            ++hyperedgePosition;
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
                        parseBagLine(line, decomposition, graph, notfoundEdges);
                        break;

                    default: // edge line
                        parseEdgeLine(line, decomposition);
                        break;
                }
            }
        }

        ConstIterator<vertex_t> iter;

        // check if all vertices are in a bag
        for (iter = graph.vertices().begin(); iter != graph.vertices().end(); ++iter)
        {
            bool found = false;

            for (ConstIterator<vertex_t> iter_decomp = decomposition->vertices().begin(); iter_decomp != decomposition->vertices().end(); ++iter_decomp)
            {
                std::vector<vertex_t> & bagContent = decomposition->mutableBagContent(*iter_decomp);

                if (std::find(bagContent.begin(), bagContent.end(), *iter) != bagContent.end())
                {
                    found = true;

                    break;
                }
            }
            if (!found)
            {
                delete decomposition;

                decomposition = nullptr;

                return std::make_pair(decomposition, iterations);
            }
        }

        // check edges
        if (notfoundEdges.size() > 0)
        {
            delete decomposition;

            decomposition = nullptr;

            return std::make_pair(decomposition, iterations);
        };

        // check bag size
        for (ConstIterator<vertex_t> iter_decomp = decomposition->vertices().begin(); iter_decomp != decomposition->vertices().end(); ++iter_decomp)
        {
            if (decomposition->mutableBagContent(*iter_decomp).size() > maxBagSize)
            {
                delete decomposition;

                decomposition = nullptr;

                return std::make_pair(decomposition, iterations);
            };
        }
    }
    else
    {
        if (!managementInstance.isTerminated())
        {
            decomposition->addVertex();
        }
    }
    return std::make_pair(decomposition, iterations);
}

void htd::FileGraphDecompositionAlgorithm::Implementation::parseBagLine(std::string line, IMutableGraphDecomposition * decomp, const IMultiHypergraph & graph, std::unordered_set<std::vector<vertex_t>> & notfoundEdges) const
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
        }
    }
    getInducedEdges(buckets, graph, inducedEdges, notfoundEdges);

    decomp->addVertex(std::vector<htd::vertex_t>(buckets), graph.hyperedgesAtPositions(inducedEdges));
}

void htd::FileGraphDecompositionAlgorithm::Implementation::getInducedEdges(std::vector<vertex_t> & bag, const IMultiHypergraph & graph, std::vector<index_t> & inducedEdges, std::unordered_set<std::vector<vertex_t>> & notfoundEdges) const
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

                if (std::find(bag.begin(), bag.end(), vertex) != bag.end())
                {
                    inducedEdges.push_back(index);

                    notfoundEdges.erase(elements);
                }

                break;
            }
            case 2:
            {
                htd::vertex_t vertex1 = elements[0];

                htd::vertex_t vertex2 = elements[1];

                if (std::find(bag.begin(), bag.end(), vertex1) != bag.end() && std::find(bag.begin(), bag.end(), vertex2) != bag.end())
                {
                    inducedEdges.push_back(index);

                    notfoundEdges.erase(elements);
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

                        notfoundEdges.erase(elements);

                        break;
                    }
                    if (std::find(bag.begin(), bag.end(), elements[i]) == bag.end())
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

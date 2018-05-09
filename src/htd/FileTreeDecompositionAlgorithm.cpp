#ifndef HTD_FILETREEDECOMPOSITIONALGORITHM_CPP
#define HTD_FILETREEDECOMPOSITIONALGORITHM_CPP

#include <htd/Globals.hpp>
#include <htd/Helpers.hpp>

#include <htd/ConnectedComponentAlgorithmFactory.hpp>
#include <htd/GraphDecompositionFactory.hpp>
#include <htd/TreeDecompositionFactory.hpp>
#include <htd/GraphLabeling.hpp>
#include <htd/ILabelingFunction.hpp>
#include <htd/CompressionOperation.hpp>
#include <htd/BreadthFirstGraphTraversal.hpp>
#include <htd/GraphPreprocessorFactory.hpp>
#include <htd/FileTreeDecompositionAlgorithm.hpp>

#include <cstdarg>
#include <sstream>
#include <fstream>
#include <htd/FileGraphDecompositionAlgorithm.hpp>

/**
 *  Private implementation details of class htd::FileTreeDecompostionAlgorithm.
 */
struct htd::FileTreeDecompositionAlgorithm::Implementation
{
    /**
     *  Constructor for the implementation details structure.
     *
     *  @param[in] manager  The management instance to which the current object instance belongs.
     */
    Implementation(const htd::LibraryInstance * const manager, const std::string & decomposition, const bool & isPath) : managementInstance_(manager), labelingFunctions_(), postProcessingOperations_()
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
        baseAlgorithm_ = new htd::FileGraphDecompositionAlgorithm(manager, this->decomposition, false);
    }

    /**
     *  Copy constructor for the implementation details structure.
     *
     *  @param[in] original The original implementation details structure.
     */
    Implementation(const Implementation & original) : decomposition(original.decomposition), managementInstance_(original.managementInstance_), baseAlgorithm_(original.baseAlgorithm_->clone()), labelingFunctions_(), postProcessingOperations_()
    {
        for (htd::ILabelingFunction * labelingFunction : original.labelingFunctions_)
        {
#ifndef HTD_USE_VISUAL_STUDIO_COMPATIBILITY_MODE
            labelingFunctions_.push_back(labelingFunction->clone());
#else
            labelingFunctions_.push_back(labelingFunction->cloneLabelingFunction());
#endif
        }

        for (htd::ITreeDecompositionManipulationOperation * postProcessingOperation : original.postProcessingOperations_)
        {
#ifndef HTD_USE_VISUAL_STUDIO_COMPATIBILITY_MODE
            postProcessingOperations_.push_back(postProcessingOperation->clone());
#else
            postProcessingOperations_.push_back(postProcessingOperation->cloneTreeDecompositionManipulationOperation());
#endif
        }
    }

    virtual ~Implementation()
    {
        delete baseAlgorithm_;

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
     *  The path to the file containing the tree decompostion.
     */
    std::string decomposition;

    /**
     *  The management instance to which the current object instance belongs.
     */
    const htd::LibraryInstance * managementInstance_;

    /**
     *  The underlying graph decomposition algorithm.
     */
    htd::FileGraphDecompositionAlgorithm * baseAlgorithm_;

    /**
     *  The labeling functions which are applied after the decomposition was computed.
     */
    std::vector<htd::ILabelingFunction *> labelingFunctions_;

    /**
     *  The manipuation operations which are applied after the decomposition was computed.
     */
    std::vector<htd::ITreeDecompositionManipulationOperation *> postProcessingOperations_;

    /**
     *  Compute a new mutable tree decompostion of the given graph.
     *
     *  @param[in] graph                The graph which shall be decomposed.
     *  @param[in] preprocessedGraph    The input graph in preprocessed format.
     *  @param[in] maxBagSize           The upper bound for the maximum bag size of the decomposition.
     *  @param[in] maxIterationCount    The maximum number of iterations resulting in a higher maximum bag size than maxBagSize after which a null-pointer is returned.
     *
     *  @return A pair consisting of a mutable tree decompostion of the given graph or a null-pointer in case that no decomposition with a appropriate maximum bag size could be found after maxIterationCount iterations and the number of iterations actually needed to find the decomposition at hand.
     */
    std::pair<htd::IMutableTreeDecomposition *, std::size_t> computeMutableDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, std::size_t maxBagSize, std::size_t maxIterationCount) const;
};

htd::FileTreeDecompositionAlgorithm::FileTreeDecompositionAlgorithm(const htd::LibraryInstance * const manager, const std::string & decomposition, const bool & isPath) : implementation_(new Implementation(manager, decomposition, isPath))
{

}

htd::FileTreeDecompositionAlgorithm::FileTreeDecompositionAlgorithm(const htd::LibraryInstance * const manager, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations, const std::string & decomposition, const bool & isPath) : implementation_(new Implementation(manager, decomposition, isPath))
{
    setManipulationOperations(manipulationOperations);
}

htd::FileTreeDecompositionAlgorithm::FileTreeDecompositionAlgorithm(const htd::FileTreeDecompositionAlgorithm & original) : implementation_(new Implementation(*(original.implementation_)))
{

}

htd::FileTreeDecompositionAlgorithm::~FileTreeDecompositionAlgorithm()
{

}

htd::ITreeDecomposition * htd::FileTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph) const
{
    return computeDecomposition(graph, std::vector<htd::IDecompositionManipulationOperation *>());
}

htd::ITreeDecomposition * htd::FileTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations) const
{
    return computeDecomposition(graph, manipulationOperations, (std::size_t) -1, 1).first;
}

std::pair<htd::ITreeDecomposition *, std::size_t> htd::FileTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, std::size_t maxBagSize, std::size_t maxIterationCount) const
{
    return computeDecomposition(graph, std::vector<htd::IDecompositionManipulationOperation *>(), maxBagSize, maxIterationCount);
}

std::pair<htd::ITreeDecomposition *, std::size_t> htd::FileTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations, std::size_t maxBagSize, std::size_t maxIterationCount) const
{
    htd::IGraphPreprocessor * preprocessor = implementation_->managementInstance_->graphPreprocessorFactory().createInstance();

    htd::IPreprocessedGraph * preprocessedGraph = preprocessor->prepare(graph);

    std::pair<htd::ITreeDecomposition *, std::size_t> ret =
        computeDecomposition(graph, *preprocessedGraph, manipulationOperations, maxBagSize, maxIterationCount);

    delete preprocessedGraph;
    delete preprocessor;

    return ret;
}

htd::ITreeDecomposition * htd::FileTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph) const
{
    return computeDecomposition(graph, preprocessedGraph, std::vector<htd::IDecompositionManipulationOperation *>());
}

htd::ITreeDecomposition * htd::FileTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations) const
{
    return computeDecomposition(graph, preprocessedGraph, manipulationOperations, (std::size_t) -1, 1).first;
}

std::pair<htd::ITreeDecomposition *, std::size_t> htd::FileTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, std::size_t maxBagSize, std::size_t maxIterationCount) const
{
    return computeDecomposition(graph, preprocessedGraph, std::vector<htd::IDecompositionManipulationOperation *>(), maxBagSize, maxIterationCount);
}

std::pair<htd::ITreeDecomposition *, std::size_t> htd::FileTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations, std::size_t maxBagSize, std::size_t maxIterationCount) const
{
    std::pair<htd::IMutableTreeDecomposition *, std::size_t> ret = implementation_->computeMutableDecomposition(graph, preprocessedGraph, maxBagSize, maxIterationCount);

    htd::IMutableTreeDecomposition * decomposition = ret.first;

    if (decomposition != nullptr)
    {
        std::vector<htd::ILabelingFunction *> labelingFunctions;

        std::vector<htd::ITreeDecompositionManipulationOperation *> postProcessingOperations;

        for (htd::IDecompositionManipulationOperation * operation : manipulationOperations)
        {
            htd::ILabelingFunction * labelingFunction = dynamic_cast<htd::ILabelingFunction *>(operation);

            if (labelingFunction != nullptr)
            {
                labelingFunctions.push_back(labelingFunction);
            }

            htd::ITreeDecompositionManipulationOperation * manipulationOperation = dynamic_cast<htd::ITreeDecompositionManipulationOperation *>(operation);

            if (manipulationOperation != nullptr)
            {
                postProcessingOperations.push_back(manipulationOperation);
            }
        }

        for (const htd::ITreeDecompositionManipulationOperation * operation : implementation_->postProcessingOperations_)
        {
            operation->apply(graph, *decomposition);
        }

        for (htd::ITreeDecompositionManipulationOperation * operation : postProcessingOperations)
        {
            operation->apply(graph, *decomposition);
        }

        for (const htd::ILabelingFunction * labelingFunction : implementation_->labelingFunctions_)
        {
            for (htd::vertex_t vertex : decomposition->vertices())
            {
                htd::ILabelCollection * labelCollection = decomposition->labelings().exportVertexLabelCollection(vertex);

                htd::ILabel * newLabel = labelingFunction->computeLabel(decomposition->bagContent(vertex), *labelCollection);

                delete labelCollection;

                decomposition->setVertexLabel(labelingFunction->name(), vertex, newLabel);
            }
        }

        for (htd::ILabelingFunction * labelingFunction : labelingFunctions)
        {
            for (htd::vertex_t vertex : decomposition->vertices())
            {
                htd::ILabelCollection * labelCollection = decomposition->labelings().exportVertexLabelCollection(vertex);

                htd::ILabel * newLabel = labelingFunction->computeLabel(decomposition->bagContent(vertex), *labelCollection);

                delete labelCollection;

                decomposition->setVertexLabel(labelingFunction->name(), vertex, newLabel);
            }
        }
    }

    for (htd::IDecompositionManipulationOperation * operation : manipulationOperations)
    {
        delete operation;
    }

    return ret;
}

htd::ITreeDecomposition * htd::FileTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, int manipulationOperationCount, ...) const
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

htd::ITreeDecomposition * htd::FileTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, int manipulationOperationCount, ...) const
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

void htd::FileTreeDecompositionAlgorithm::setOrderingAlgorithm(htd::IOrderingAlgorithm * algorithm)
{
    implementation_->baseAlgorithm_->setOrderingAlgorithm(algorithm);
}

void htd::FileTreeDecompositionAlgorithm::setManipulationOperations(const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations)
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

void htd::FileTreeDecompositionAlgorithm::addManipulationOperation(htd::IDecompositionManipulationOperation * manipulationOperation)
{
    bool assigned = false;

    htd::ILabelingFunction * labelingFunction = dynamic_cast<htd::ILabelingFunction *>(manipulationOperation);

    if (labelingFunction != nullptr)
    {
        implementation_->labelingFunctions_.emplace_back(labelingFunction);

        assigned = true;
    }

    htd::ITreeDecompositionManipulationOperation * newManipulationOperation = dynamic_cast<htd::ITreeDecompositionManipulationOperation *>(manipulationOperation);

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

void htd::FileTreeDecompositionAlgorithm::addManipulationOperations(const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations)
{
    for (htd::IDecompositionManipulationOperation * operation : manipulationOperations)
    {
        addManipulationOperation(operation);
    }
}

bool htd::FileTreeDecompositionAlgorithm::isSafelyInterruptible(void) const
{
    return false;
}

const htd::LibraryInstance * htd::FileTreeDecompositionAlgorithm::managementInstance(void) const HTD_NOEXCEPT
{
    return implementation_->managementInstance_;
}

void htd::FileTreeDecompositionAlgorithm::setManagementInstance(const htd::LibraryInstance * const manager)
{
    HTD_ASSERT(manager != nullptr)

    implementation_->managementInstance_ = manager;
}

bool htd::FileTreeDecompositionAlgorithm::isCompressionEnabled(void) const
{
    return implementation_->baseAlgorithm_->isCompressionEnabled();
}

void htd::FileTreeDecompositionAlgorithm::setCompressionEnabled(bool compressionEnabled)
{
    implementation_->baseAlgorithm_->setCompressionEnabled(compressionEnabled);
}

bool htd::FileTreeDecompositionAlgorithm::isComputeInducedEdgesEnabled(void) const
{
    return implementation_->baseAlgorithm_->isComputeInducedEdgesEnabled();
}

void htd::FileTreeDecompositionAlgorithm::setComputeInducedEdgesEnabled(bool computeInducedEdgesEnabled)
{
    implementation_->baseAlgorithm_->setComputeInducedEdgesEnabled(computeInducedEdgesEnabled);
}

htd::FileTreeDecompositionAlgorithm * htd::FileTreeDecompositionAlgorithm::clone(void) const
{
    return new htd::FileTreeDecompositionAlgorithm(*this);
}

std::pair<htd::IMutableTreeDecomposition *, std::size_t> htd::FileTreeDecompositionAlgorithm::Implementation::computeMutableDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, std::size_t maxBagSize, std::size_t maxIterationCount) const
{
    htd::IMutableTreeDecomposition * ret = managementInstance_->treeDecompositionFactory().createInstance();

    std::size_t iterations = 1;

    if (graph.vertexCount() > 0)
    {
        std::pair<htd::IGraphDecomposition *, std::size_t> graphDecomposition = baseAlgorithm_->computeDecomposition(graph, preprocessedGraph, maxBagSize, maxIterationCount);

        if (graphDecomposition.first != nullptr)
        {
            htd::IMutableGraphDecomposition & mutableGraphDecomposition = managementInstance_->graphDecompositionFactory().accessMutableInstance(*(graphDecomposition.first));

            if (!managementInstance_->isTerminated())
            {
                if (mutableGraphDecomposition.edgeCount() + 1 != mutableGraphDecomposition.vertexCount() || mutableGraphDecomposition.isolatedVertexCount() > 0)
                {
                    htd::IConnectedComponentAlgorithm * connectedComponentAlgorithm = managementInstance_->connectedComponentAlgorithmFactory().createInstance();

                    HTD_ASSERT(connectedComponentAlgorithm != nullptr)

                    std::vector<std::vector<htd::vertex_t>> components;

                    connectedComponentAlgorithm->determineComponents(*(graphDecomposition.first), components);

                    delete connectedComponentAlgorithm;

                    std::size_t componentCount = components.size();

                    if (componentCount > 1)
                    {
                        for (htd::index_t index = 0; index < componentCount - 1; ++index)
                        {
                            const std::vector<htd::vertex_t> & component1 = components[index];

                            const std::vector<htd::vertex_t> & component2 = components[index + 1];

                            /* Coverity complains about std::rand() being not safe for security related operations. We are happy with a pseudo-random number here. */
                            // coverity[dont_call]
                            htd::vertex_t vertex1 = component1[std::rand() % component1.size()];

                            /* Coverity complains about std::rand() being not safe for security related operations. We are happy with a pseudo-random number here. */
                            // coverity[dont_call]
                            htd::vertex_t vertex2 = component2[std::rand() % component2.size()];

                            mutableGraphDecomposition.addEdge(vertex1, vertex2);
                        }
                    }
                }

                htd::vertex_t node = htd::Vertex::UNKNOWN;

                std::unordered_map<htd::vertex_t, htd::vertex_t> vertexMapping;

                htd::BreadthFirstGraphTraversal graphTraversal(managementInstance_);

                graphTraversal.traverse(*(graphDecomposition.first), graphDecomposition.first->vertexAtPosition(0), [&](htd::vertex_t vertex, htd::vertex_t predecessor, std::size_t distanceFromStartingVertex)
                {
                    HTD_UNUSED(distanceFromStartingVertex)

                    if (predecessor == htd::Vertex::UNKNOWN)
                    {
                        node = ret->insertRoot(std::move(mutableGraphDecomposition.mutableBagContent(vertex)), std::move(mutableGraphDecomposition.mutableInducedHyperedges(vertex)));
                    }
                    else
                    {
                        node = ret->addChild(vertexMapping.at(predecessor), std::move(mutableGraphDecomposition.mutableBagContent(vertex)), std::move(mutableGraphDecomposition.mutableInducedHyperedges(vertex)));
                    }

                    vertexMapping.emplace(vertex, node);
                });
            }
            else
            {
                delete ret;

                ret = nullptr;
            }

            delete graphDecomposition.first;
        }
        else
        {
            delete ret;

            ret = nullptr;
        }

        iterations = graphDecomposition.second;
    }
    else
    {
        ret->insertRoot();
    }

    return std::make_pair(ret, iterations);
}


#endif /* HTD_FILETREEDECOMPOSITIONALGORITHM_CPP */

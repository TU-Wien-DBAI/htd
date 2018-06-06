#ifndef HTD_FITNESSMINIMIZINGTREEDECOMPOSITIONALGORITHM_CPP
#define HTD_FITNESSMINIMIZINGTREEDECOMPOSITIONALGORITHM_CPP

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
#include <htd/FitnessMinimizingTreeDecompositionAlgorithm.hpp>

#include <cstdarg>
#include <sstream>
#include <fstream>
#include <htd/FileTreeDecompositionAlgorithm.hpp>
#include <htd_io/GrFormatExporter.hpp>

#include <unistd.h>
#include <sys/wait.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <poll.h>
#include <fcntl.h>
#include <wordexp.h>
#include <htd/ITreeDecompositionFitnessFunction.hpp>

/**
 *  Private implementation details of class htd::FileTreeDecompostionAlgorithm.
 */
struct htd::FitnessMinimizingTreeDecompositionAlgorithm::Implementation
{
    /**
     *  Constructor for the implementation details structure.
     *
     *  @param[in] manager          The management instance to which the current object instance belongs.
     *  @param[in] decompostion     String containing the tree decomposition or the path to the file containing the tree decomposition.
     */
    Implementation(const htd::LibraryInstance * const manager, ITreeDecompositionFitnessFunction * fitnessFunction, std::vector<ITreeDecompositionAlgorithm *> deterministicAlgorithms, std::vector<ITreeDecompositionAlgorithm *> nonDeterministicAlgorithms) : managementInstance_(manager), labelingFunctions_(), postProcessingOperations_(), deterministicAlgorithms_(deterministicAlgorithms), nonDeterministicAlgorithms_(nonDeterministicAlgorithms), fitnessFunction_(fitnessFunction)
    {
    }

    /**
     *  Copy constructor for the implementation details structure.
     *
     *  @param[in] original The original implementation details structure.
     */
    Implementation(const Implementation & original) : managementInstance_(original.managementInstance_), labelingFunctions_(), postProcessingOperations_(), deterministicAlgorithms_(original.deterministicAlgorithms_), nonDeterministicAlgorithms_(original.nonDeterministicAlgorithms_)
    {
        fitnessFunction_ = original.fitnessFunction_;
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
    std::vector<htd::ITreeDecompositionManipulationOperation *> postProcessingOperations_;

    /**
     *  A boolean flag indicating whether the hyperedges induced by a respective bag shall be computed.
     */
    bool computeInducedEdges_;

    /**
     * TODO
     */
    std::vector<ITreeDecompositionAlgorithm *> deterministicAlgorithms_;

    /**
     * TODO
     */
    std::vector<ITreeDecompositionAlgorithm *> nonDeterministicAlgorithms_;

    /**
     * TODO
     */
    ITreeDecompositionFitnessFunction * fitnessFunction_;

    /**
     *  Compute a new mutable tree decompostion of the given graph.
     *
     *  @param[in] graph                The graph which shall be decomposed.
     *  @param[in] preprocessedGraph    The input graph in preprocessed format.
     *  @param[in] maxBagSize           The upper bound for the maximum bag size of the decomposition.
     *  @param[in] maxIterationCount    The maximum number of iterations resulting in a higher maximum bag size than maxBagSize after which a null-pointer is returned.
     *
     *  @return A pair consisting of a mutable tree decompostion of the given graph or a null-pointer in case that the decomposition does not have a appropriate maximum bag size or the decomposition is not a valid decomposition of the graph.
     */
    std::pair<htd::IMutableTreeDecomposition *, std::size_t> computeMutableDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, std::size_t maxBagSize, std::size_t maxIterationCount) const;

    std::string convert(const IMultiHypergraph & graph) const;
};

htd::FitnessMinimizingTreeDecompositionAlgorithm::FitnessMinimizingTreeDecompositionAlgorithm(const htd::LibraryInstance * const manager, ITreeDecompositionFitnessFunction * fitnessFunction, std::vector<ITreeDecompositionAlgorithm *> deterministicAlgorithms, std::vector<ITreeDecompositionAlgorithm *> nonDeterministicAlgorithms) : implementation_(new Implementation(manager, fitnessFunction, deterministicAlgorithms, nonDeterministicAlgorithms))
{

}

htd::FitnessMinimizingTreeDecompositionAlgorithm::FitnessMinimizingTreeDecompositionAlgorithm(const htd::LibraryInstance * const manager, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations, ITreeDecompositionFitnessFunction * fitnessFunction, std::vector<ITreeDecompositionAlgorithm *> deterministicAlgorithms, std::vector<ITreeDecompositionAlgorithm *> nonDeterministicAlgorithms) : implementation_(new Implementation(manager, fitnessFunction, deterministicAlgorithms, nonDeterministicAlgorithms))
{
    setManipulationOperations(manipulationOperations);
}

htd::FitnessMinimizingTreeDecompositionAlgorithm::FitnessMinimizingTreeDecompositionAlgorithm(const htd::FitnessMinimizingTreeDecompositionAlgorithm & original) : implementation_(new Implementation(*(original.implementation_)))
{

}

htd::FitnessMinimizingTreeDecompositionAlgorithm::~FitnessMinimizingTreeDecompositionAlgorithm()
{

}

htd::ITreeDecomposition * htd::FitnessMinimizingTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph) const
{
    return computeDecomposition(graph, std::vector<htd::IDecompositionManipulationOperation *>());
}

htd::ITreeDecomposition * htd::FitnessMinimizingTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations) const
{
    return computeDecomposition(graph, manipulationOperations, (std::size_t)-1, 1).first;
}

std::pair<htd::ITreeDecomposition *, std::size_t> htd::FitnessMinimizingTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, std::size_t maxBagSize, std::size_t maxIterationCount) const
{
    return computeDecomposition(graph, std::vector<htd::IDecompositionManipulationOperation *>(), maxBagSize, maxIterationCount);
}

std::pair<htd::ITreeDecomposition *, std::size_t> htd::FitnessMinimizingTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations, std::size_t maxBagSize, std::size_t maxIterationCount) const
{
    htd::IGraphPreprocessor * preprocessor = implementation_->managementInstance_->graphPreprocessorFactory().createInstance();

    htd::IPreprocessedGraph * preprocessedGraph = preprocessor->prepare(graph);

    std::pair<htd::ITreeDecomposition *, std::size_t> ret =
        computeDecomposition(graph, *preprocessedGraph, manipulationOperations, maxBagSize, maxIterationCount);

    delete preprocessedGraph;
    delete preprocessor;

    return ret;
}

htd::ITreeDecomposition * htd::FitnessMinimizingTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph) const
{
    return computeDecomposition(graph, preprocessedGraph, std::vector<htd::IDecompositionManipulationOperation *>());
}

htd::ITreeDecomposition * htd::FitnessMinimizingTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations) const
{
    return computeDecomposition(graph, preprocessedGraph, manipulationOperations, (std::size_t)-1, 1).first;
}

std::pair<htd::ITreeDecomposition *, std::size_t> htd::FitnessMinimizingTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, std::size_t maxBagSize, std::size_t maxIterationCount) const
{
    return computeDecomposition(graph, preprocessedGraph, std::vector<htd::IDecompositionManipulationOperation *>(), maxBagSize, maxIterationCount);
}

std::pair<htd::ITreeDecomposition *, std::size_t> htd::FitnessMinimizingTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations, std::size_t maxBagSize, std::size_t maxIterationCount) const
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

htd::ITreeDecomposition * htd::FitnessMinimizingTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, int manipulationOperationCount, ...) const
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

htd::ITreeDecomposition * htd::FitnessMinimizingTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, int manipulationOperationCount, ...) const
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

void htd::FitnessMinimizingTreeDecompositionAlgorithm::setManipulationOperations(const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations)
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

void htd::FitnessMinimizingTreeDecompositionAlgorithm::addManipulationOperation(htd::IDecompositionManipulationOperation * manipulationOperation)
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

void htd::FitnessMinimizingTreeDecompositionAlgorithm::addManipulationOperations(const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations)
{
    for (htd::IDecompositionManipulationOperation * operation : manipulationOperations)
    {
        addManipulationOperation(operation);
    }
}

bool htd::FitnessMinimizingTreeDecompositionAlgorithm::isSafelyInterruptible(void) const
{
    return false;
}

const htd::LibraryInstance * htd::FitnessMinimizingTreeDecompositionAlgorithm::managementInstance(void) const HTD_NOEXCEPT
{
    return implementation_->managementInstance_;
}

void htd::FitnessMinimizingTreeDecompositionAlgorithm::setManagementInstance(const htd::LibraryInstance * const manager)
{
    HTD_ASSERT(manager != nullptr)

    implementation_->managementInstance_ = manager;
}

bool htd::FitnessMinimizingTreeDecompositionAlgorithm::isComputeInducedEdgesEnabled(void) const
{
    return implementation_->computeInducedEdges_;
}

void htd::FitnessMinimizingTreeDecompositionAlgorithm::setComputeInducedEdgesEnabled(bool computeInducedEdgesEnabled)
{
    implementation_->computeInducedEdges_ = computeInducedEdgesEnabled;
}

htd::FitnessMinimizingTreeDecompositionAlgorithm * htd::FitnessMinimizingTreeDecompositionAlgorithm::clone(void) const
{
    return new htd::FitnessMinimizingTreeDecompositionAlgorithm(*this);
}

std::pair<htd::IMutableTreeDecomposition *, std::size_t> htd::FitnessMinimizingTreeDecompositionAlgorithm::Implementation::computeMutableDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph &, std::size_t maxBagSize, std::size_t maxIterationCount) const
{
    ITreeDecomposition * decomp = nullptr;

    for (const ITreeDecompositionAlgorithm * deterministicAlgorithm : deterministicAlgorithms_)
    {
        ITreeDecomposition * tmpDecomp = deterministicAlgorithm->computeDecomposition(graph);

        if (decomp == 0 || fitnessFunction_->fitness(graph, *tmpDecomp) < fitnessFunction_->fitness(graph, *decomp))
        {
            delete decomp;

            decomp = tmpDecomp;
        }
    }

    for (std::size_t i = 0; i < maxIterationCount; i++)
    {
        for (const ITreeDecompositionAlgorithm * nonDeterministicAlgorithm : nonDeterministicAlgorithms_)
        {
            ITreeDecomposition * tmpDecomp = nonDeterministicAlgorithm->computeDecomposition(graph);

            if (decomp == 0 || fitnessFunction_->fitness(graph, *tmpDecomp) < fitnessFunction_->fitness(graph, *decomp))
            {
                delete decomp;

                decomp = tmpDecomp;
            }
        }
    }

    htd::IMutableTreeDecomposition & mutableTreeDecomposition = managementInstance_->treeDecompositionFactory().accessMutableInstance(*decomp);

    return std::pair<htd::IMutableTreeDecomposition *, std::size_t>(&mutableTreeDecomposition, 0);
}


#endif /* HTD_FITNESSMINIMIZINGTREEDECOMPOSITIONALGORITHM_CPP */

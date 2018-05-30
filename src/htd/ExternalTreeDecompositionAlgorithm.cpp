#ifndef HTD_TRELLISTREEDECOMPOSITIONALGORITHM_CPP
#define HTD_TRELLISTREEDECOMPOSITIONALGORITHM_CPP

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
#include <htd/ExternalTreeDecompositionAlgorithm.hpp>

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

/**
 *  Private implementation details of class htd::FileTreeDecompostionAlgorithm.
 */
struct htd::ExternalTreeDecompositionAlgorithm::Implementation
{
    /**
     *  Constructor for the implementation details structure.
     *
     *  @param[in] manager          The management instance to which the current object instance belongs.
     *  @param[in] decompostion     String containing the tree decomposition or the path to the file containing the tree decomposition.
     */
    Implementation(const htd::LibraryInstance * const manager, char ** cmd, __useconds_t timeout) : managementInstance_(manager), labelingFunctions_(), postProcessingOperations_(), cmd_(cmd), timeout_(timeout)
    {

    }

    /**
     *  Copy constructor for the implementation details structure.
     *
     *  @param[in] original The original implementation details structure.
     */
    Implementation(const Implementation & original) : managementInstance_(original.managementInstance_), labelingFunctions_(), postProcessingOperations_(), cmd_(original.cmd_), timeout_(original.timeout_)
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
     *  The string containing the decompostion in td format.
     */
    char ** cmd_;

    __useconds_t timeout_;

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

htd::ExternalTreeDecompositionAlgorithm::ExternalTreeDecompositionAlgorithm(const htd::LibraryInstance * const manager, char ** cmd, __useconds_t timeout) : implementation_(new Implementation(manager, cmd, timeout))
{

}

htd::ExternalTreeDecompositionAlgorithm::ExternalTreeDecompositionAlgorithm(const htd::LibraryInstance * const manager, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations, char ** cmd, __useconds_t timeout) : implementation_(new Implementation(manager, cmd, timeout))
{
    setManipulationOperations(manipulationOperations);
}

htd::ExternalTreeDecompositionAlgorithm::ExternalTreeDecompositionAlgorithm(const htd::ExternalTreeDecompositionAlgorithm & original) : implementation_(new Implementation(*(original.implementation_)))
{

}

htd::ExternalTreeDecompositionAlgorithm::~ExternalTreeDecompositionAlgorithm()
{

}

htd::ITreeDecomposition * htd::ExternalTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph) const
{
    return computeDecomposition(graph, std::vector<htd::IDecompositionManipulationOperation *>());
}

htd::ITreeDecomposition * htd::ExternalTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations) const
{
    return computeDecomposition(graph, manipulationOperations, (std::size_t)-1, 1).first;
}

std::pair<htd::ITreeDecomposition *, std::size_t> htd::ExternalTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, std::size_t maxBagSize, std::size_t maxIterationCount) const
{
    return computeDecomposition(graph, std::vector<htd::IDecompositionManipulationOperation *>(), maxBagSize, maxIterationCount);
}

std::pair<htd::ITreeDecomposition *, std::size_t> htd::ExternalTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations, std::size_t maxBagSize, std::size_t maxIterationCount) const
{
    htd::IGraphPreprocessor * preprocessor = implementation_->managementInstance_->graphPreprocessorFactory().createInstance();

    htd::IPreprocessedGraph * preprocessedGraph = preprocessor->prepare(graph);

    std::pair<htd::ITreeDecomposition *, std::size_t> ret =
        computeDecomposition(graph, *preprocessedGraph, manipulationOperations, maxBagSize, maxIterationCount);

    delete preprocessedGraph;
    delete preprocessor;

    return ret;
}

htd::ITreeDecomposition * htd::ExternalTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph) const
{
    return computeDecomposition(graph, preprocessedGraph, std::vector<htd::IDecompositionManipulationOperation *>());
}

htd::ITreeDecomposition * htd::ExternalTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations) const
{
    return computeDecomposition(graph, preprocessedGraph, manipulationOperations, (std::size_t)-1, 1).first;
}

std::pair<htd::ITreeDecomposition *, std::size_t> htd::ExternalTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, std::size_t maxBagSize, std::size_t maxIterationCount) const
{
    return computeDecomposition(graph, preprocessedGraph, std::vector<htd::IDecompositionManipulationOperation *>(), maxBagSize, maxIterationCount);
}

std::pair<htd::ITreeDecomposition *, std::size_t> htd::ExternalTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations, std::size_t maxBagSize, std::size_t maxIterationCount) const
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

htd::ITreeDecomposition * htd::ExternalTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, int manipulationOperationCount, ...) const
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

htd::ITreeDecomposition * htd::ExternalTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, int manipulationOperationCount, ...) const
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

void htd::ExternalTreeDecompositionAlgorithm::setManipulationOperations(const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations)
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

void htd::ExternalTreeDecompositionAlgorithm::addManipulationOperation(htd::IDecompositionManipulationOperation * manipulationOperation)
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

void htd::ExternalTreeDecompositionAlgorithm::addManipulationOperations(const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations)
{
    for (htd::IDecompositionManipulationOperation * operation : manipulationOperations)
    {
        addManipulationOperation(operation);
    }
}

bool htd::ExternalTreeDecompositionAlgorithm::isSafelyInterruptible(void) const
{
    return false;
}

const htd::LibraryInstance * htd::ExternalTreeDecompositionAlgorithm::managementInstance(void) const HTD_NOEXCEPT
{
    return implementation_->managementInstance_;
}

void htd::ExternalTreeDecompositionAlgorithm::setManagementInstance(const htd::LibraryInstance * const manager)
{
    HTD_ASSERT(manager != nullptr)

    implementation_->managementInstance_ = manager;
}

bool htd::ExternalTreeDecompositionAlgorithm::isComputeInducedEdgesEnabled(void) const
{
    return implementation_->computeInducedEdges_;
}

void htd::ExternalTreeDecompositionAlgorithm::setComputeInducedEdgesEnabled(bool computeInducedEdgesEnabled)
{
    implementation_->computeInducedEdges_ = computeInducedEdgesEnabled;
}

htd::ExternalTreeDecompositionAlgorithm * htd::ExternalTreeDecompositionAlgorithm::clone(void) const
{
    return new htd::ExternalTreeDecompositionAlgorithm(*this);
}

std::string htd::ExternalTreeDecompositionAlgorithm::Implementation::convert(const htd::IMultiHypergraph & graph) const
{
    std::unordered_map<htd::vertex_t, std::size_t> indices;
    std::string graphString;

    graphString += "p tw " + std::to_string(graph.vertexCount()) + " " + std::to_string(graph.edgeCount()) + "\n";

    if (graph.vertexCount() > 0)
    {
        std::stringstream tmpStream;

        const htd::ConstCollection<htd::Hyperedge> & hyperedgeCollection = graph.hyperedges();

        std::size_t edgeCount = graph.edgeCount();

        auto it = hyperedgeCollection.begin();

        for (htd::index_t index = 0; index < edgeCount; ++index)
        {
            htd::vertex_t vertex1 = (*it)[0];
            htd::vertex_t vertex2 = (*it)[1];

            graphString += std::to_string(vertex1) + " " + std::to_string(vertex2) + "\n";

            ++it;
        }
    }
    return graphString;
}

std::pair<htd::IMutableTreeDecomposition *, std::size_t> htd::ExternalTreeDecompositionAlgorithm::Implementation::computeMutableDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, std::size_t maxBagSize, std::size_t maxIterationCount) const
{
    const char * graphString = convert(graph).c_str();

    std::string decompString;

    pid_t pid = 0;
    int inpipefd[2];
    int outpipefd[2];
    char buf[256];

    pipe2(inpipefd, O_NONBLOCK);
    pipe2(outpipefd, O_NONBLOCK);
    pid = fork();

    if (pid == 0)
    {
        dup2(outpipefd[0], STDIN_FILENO);
        close(outpipefd[1]);
        dup2(inpipefd[1], STDOUT_FILENO);
        close(inpipefd[0]);

        execvp(cmd_[0], cmd_);
        _Exit(0);
    }
    //close unused pipe ends
    close(outpipefd[0]);
    close(inpipefd[1]);
    write(outpipefd[1], graphString, strlen(graphString));
    close(outpipefd[1]);

    ssize_t outputLen = -1;

    long starttime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    while (true)
    {
        if(outputLen<0){
            usleep(200000);
        }
        outputLen = read(inpipefd[0], buf, 256);
        if (outputLen > 0)
        {
            decompString.append(buf, outputLen);
        }
        else if (outputLen == 0)
        { break; }
        else if (timeout_ > 0 && (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() - starttime) > timeout_)
        {
            kill(pid, SIGTERM);
        }
    }
    close(inpipefd[0]);

    // call File Tree Decomposition Algorithm
    FileTreeDecompositionAlgorithm algorithm(managementInstance_, decompString);

    ITreeDecomposition * treeDecomposition = algorithm.computeDecomposition(graph);

    htd::IMutableTreeDecomposition & mutableTreeDecomposition = managementInstance_->treeDecompositionFactory().accessMutableInstance(*treeDecomposition);

    return std::make_pair(&mutableTreeDecomposition, (size_t)0);
}


#endif /* HTD_TRELLISTREEDECOMPOSITIONALGORITHM_CPP */

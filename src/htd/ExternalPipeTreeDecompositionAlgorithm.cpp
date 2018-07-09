#ifndef HTD_EXTERNALPIPETREEDECOMPOSITIONALGORITHM_CPP
#define HTD_EXTERNALPIPETREEDECOMPOSITIONALGORITHM_CPP

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
#include <htd/ExternalPipeTreeDecompositionAlgorithm.hpp>

#include <cstdarg>
#include <sstream>
#include <fstream>
#include <htd/FileTreeDecompositionAlgorithm.hpp>
#include <htd_io/GrFormatExporter.hpp>

#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <thread>

/**
*  Private implementation details of class htd::FileTreeDecompostionAlgorithm.
*/
struct htd::ExternalPipeTreeDecompositionAlgorithm::Implementation
{
    /**
    *  Constructor for the implementation details structure.
    *
    *  @param[in] manager          The management instance to which the current object instance belongs.
    *  @param[in] cmd              The command used to call the external solver.
    *  @param[in] timeout          The timeout for the external solver in milliseconds.
    */
    Implementation(const htd::LibraryInstance * const manager, const std::string & cmd, unsigned int timeout) : managementInstance_(manager), labelingFunctions_(), postProcessingOperations_(), cmd_(cmd), timeout_(timeout)
    {
    }

    /**
    *  Copy constructor for the implementation details structure.
    *
    *  @param[in] original The original implementation details structure.
    */
    Implementation(const Implementation & original) : managementInstance_(original.managementInstance_), labelingFunctions_(), postProcessingOperations_(), cmd_(original.cmd_), timeout_(original.timeout_), dir_(original.dir_)
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
    *  The command to call the external solver with.
    */
    std::string cmd_;

    /**
    * The timeout for the command in milliseconds
    */
    unsigned int timeout_;

    /**
    * The directory in which the solver gets executed
    */
    std::string dir_;

    /**
    *  Compute a new mutable tree decompostion of the given graph.
    *
    *  @param[in] graph                The graph which shall be decomposed.
    *  @param[in] preprocessedGraph    The input graph in preprocessed format.
    *  @param[in] maxBagSize           The upper bound for the maximum bag size of the decomposition.
    *
    *  @return A pair consisting of a mutable tree decompostion of the given graph or a null-pointer in case that the decomposition does not have a appropriate maximum bag size or the decomposition is not a valid decomposition of the graph.
    */
    htd::IMutableTreeDecomposition * computeMutableDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, std::size_t maxBagSize) const;

    std::string convert(const IMultiHypergraph & graph) const;

    std::string getDecomp(const htd::IMultiHypergraph & graph) const;

    void setDirectory(const std::string dir)
    {
        dir_ = dir;
    }

    std::string getDirectory()
    {
        return dir_;
    }
};

void htd::ExternalPipeTreeDecompositionAlgorithm::setDirectory(const std::string dir)
{
    implementation_->setDirectory(dir);
}

std::string htd::ExternalPipeTreeDecompositionAlgorithm::getDirectory()
{
    return implementation_->getDirectory();
}

htd::ExternalPipeTreeDecompositionAlgorithm::ExternalPipeTreeDecompositionAlgorithm(const htd::LibraryInstance * const manager, std::string cmd, unsigned int timeout) : implementation_(new Implementation(manager, cmd, timeout))
{

}

htd::ExternalPipeTreeDecompositionAlgorithm::ExternalPipeTreeDecompositionAlgorithm(const htd::LibraryInstance * const manager, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations, std::string cmd, unsigned int timeout) : implementation_(new Implementation(manager, cmd, timeout))
{
    setManipulationOperations(manipulationOperations);
}

htd::ExternalPipeTreeDecompositionAlgorithm::ExternalPipeTreeDecompositionAlgorithm(const htd::ExternalPipeTreeDecompositionAlgorithm & original) : implementation_(new Implementation(*(original.implementation_)))
{

}

htd::ExternalPipeTreeDecompositionAlgorithm::~ExternalPipeTreeDecompositionAlgorithm()
{

}

htd::ITreeDecomposition * htd::ExternalPipeTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph) const
{
    return computeDecomposition(graph, std::vector<htd::IDecompositionManipulationOperation *>());
}

htd::ITreeDecomposition * htd::ExternalPipeTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations) const
{
    return computeDecomposition(graph, manipulationOperations, (std::size_t)-1);
}

htd::ITreeDecomposition * htd::ExternalPipeTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, std::size_t maxBagSize) const
{
    return computeDecomposition(graph, std::vector<htd::IDecompositionManipulationOperation *>(), maxBagSize);
}

htd::ITreeDecomposition * htd::ExternalPipeTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations, std::size_t maxBagSize) const
{
    htd::IGraphPreprocessor * preprocessor = implementation_->managementInstance_->graphPreprocessorFactory().createInstance();

    htd::IPreprocessedGraph * preprocessedGraph = preprocessor->prepare(graph);

    htd::ITreeDecomposition * ret = computeDecomposition(graph, *preprocessedGraph, manipulationOperations, maxBagSize);

    delete preprocessedGraph;
    delete preprocessor;

    return ret;
}

htd::ITreeDecomposition * htd::ExternalPipeTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph) const
{
    return computeDecomposition(graph, preprocessedGraph, std::vector<htd::IDecompositionManipulationOperation *>());
}

htd::ITreeDecomposition * htd::ExternalPipeTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations) const
{
    return computeDecomposition(graph, preprocessedGraph, manipulationOperations, (std::size_t)-1);
}

htd::ITreeDecomposition * htd::ExternalPipeTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, std::size_t maxBagSize) const
{
    return computeDecomposition(graph, preprocessedGraph, std::vector<htd::IDecompositionManipulationOperation *>(), maxBagSize);
}

htd::ITreeDecomposition * htd::ExternalPipeTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations, std::size_t maxBagSize) const
{
    htd::IMutableTreeDecomposition * decomposition = implementation_->computeMutableDecomposition(graph, preprocessedGraph, maxBagSize);

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

    return decomposition;
}

htd::ITreeDecomposition * htd::ExternalPipeTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, int manipulationOperationCount, ...) const
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

htd::ITreeDecomposition * htd::ExternalPipeTreeDecompositionAlgorithm::computeDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph & preprocessedGraph, int manipulationOperationCount, ...) const
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

void htd::ExternalPipeTreeDecompositionAlgorithm::setManipulationOperations(const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations)
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

void htd::ExternalPipeTreeDecompositionAlgorithm::addManipulationOperation(htd::IDecompositionManipulationOperation * manipulationOperation)
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

void htd::ExternalPipeTreeDecompositionAlgorithm::addManipulationOperations(const std::vector<htd::IDecompositionManipulationOperation *> & manipulationOperations)
{
    for (htd::IDecompositionManipulationOperation * operation : manipulationOperations)
    {
        addManipulationOperation(operation);
    }
}

bool htd::ExternalPipeTreeDecompositionAlgorithm::isSafelyInterruptible(void) const
{
    return false;
}

const htd::LibraryInstance * htd::ExternalPipeTreeDecompositionAlgorithm::managementInstance(void) const HTD_NOEXCEPT
{
    return implementation_->managementInstance_;
}

void htd::ExternalPipeTreeDecompositionAlgorithm::setManagementInstance(const htd::LibraryInstance * const manager)
{
    HTD_ASSERT(manager != nullptr)

    implementation_->managementInstance_ = manager;
}

bool htd::ExternalPipeTreeDecompositionAlgorithm::isComputeInducedEdgesEnabled(void) const
{
    return implementation_->computeInducedEdges_;
}

void htd::ExternalPipeTreeDecompositionAlgorithm::setComputeInducedEdgesEnabled(bool computeInducedEdgesEnabled)
{
    implementation_->computeInducedEdges_ = computeInducedEdgesEnabled;
}

htd::ExternalPipeTreeDecompositionAlgorithm * htd::ExternalPipeTreeDecompositionAlgorithm::clone(void) const
{
    return new htd::ExternalPipeTreeDecompositionAlgorithm(*this);
}

void htd::ExternalPipeTreeDecompositionAlgorithm::setCommand(const std::string & cmd)
{
    implementation_->cmd_ = cmd;
}

std::string htd::ExternalPipeTreeDecompositionAlgorithm::getCommand()
{
    return implementation_->cmd_;
}

void htd::ExternalPipeTreeDecompositionAlgorithm::setTimeout(const unsigned int & timeout)
{
    implementation_->timeout_ = timeout;
}

unsigned int htd::ExternalPipeTreeDecompositionAlgorithm::getTimeout()
{
    return implementation_->timeout_;
}

std::string htd::ExternalPipeTreeDecompositionAlgorithm::Implementation::convert(const htd::IMultiHypergraph & graph) const
{
    std::unordered_map<htd::vertex_t, std::size_t> indices;

    std::string graphString = "p tw " + std::to_string(graph.vertexCount()) + " " + std::to_string(graph.edgeCount()) + "\n";

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

htd::IMutableTreeDecomposition * htd::ExternalPipeTreeDecompositionAlgorithm::Implementation::computeMutableDecomposition(const htd::IMultiHypergraph & graph, const htd::IPreprocessedGraph &, std::size_t maxBagSize) const
{
    std::string decompString = getDecomp(graph);

    // call File Tree Decomposition Algorithm
    FileTreeDecompositionAlgorithm algorithm(managementInstance_, decompString);

    ITreeDecomposition * treeDecomposition = algorithm.computeDecomposition(graph, maxBagSize);

    if (treeDecomposition != nullptr)
    {
        htd::IMutableTreeDecomposition & mutableTreeDecomposition = managementInstance_->treeDecompositionFactory().accessMutableInstance(*treeDecomposition);

        return &mutableTreeDecomposition;
    }

    return nullptr;
}

#ifdef _MSC_VER

#include <windows.h>
#include <tchar.h>
#include <stdio.h> 
#include <strsafe.h>
#pragma comment(lib, "User32.lib")

std::string htd::ExternalPipeTreeDecompositionAlgorithm::Implementation::getDecomp(const htd::IMultiHypergraph & graph) const
{
    std::string decompString = convert(graph);

    HANDLE g_hChildStd_IN_Rd = NULL;

    HANDLE g_hChildStd_IN_Wr = NULL;

    HANDLE g_hChildStd_OUT_Rd = NULL;

    HANDLE g_hChildStd_OUT_Wr = NULL;

    SECURITY_ATTRIBUTES saAttr;

    saAttr.nLength = sizeof(SECURITY_ATTRIBUTES);

    saAttr.bInheritHandle = TRUE;

    saAttr.lpSecurityDescriptor = NULL;

    CreatePipe(&g_hChildStd_OUT_Rd, &g_hChildStd_OUT_Wr, &saAttr, 0);

    SetHandleInformation(g_hChildStd_OUT_Rd, HANDLE_FLAG_INHERIT, 0);

    CreatePipe(&g_hChildStd_IN_Rd, &g_hChildStd_IN_Wr, &saAttr, 0);

    SetHandleInformation(g_hChildStd_IN_Wr, HANDLE_FLAG_INHERIT, 0);

    TCHAR *szCmdline;

    szCmdline = (TCHAR*) malloc(sizeof(char)*(cmd_.length() + 1));

    strcpy_s(szCmdline, cmd_.length(), cmd_.c_str());

    PROCESS_INFORMATION piProcInfo;

    STARTUPINFO siStartInfo;

    BOOL bSuccess = FALSE;

    ZeroMemory(&piProcInfo, sizeof(PROCESS_INFORMATION));

    ZeroMemory(&siStartInfo, sizeof(STARTUPINFO));

    siStartInfo.cb = sizeof(STARTUPINFO);

    siStartInfo.hStdError = g_hChildStd_OUT_Wr;

    siStartInfo.hStdOutput = g_hChildStd_OUT_Wr;

    siStartInfo.hStdInput = g_hChildStd_IN_Rd;

    siStartInfo.dwFlags |= STARTF_USESTDHANDLES;

    CreateProcess(NULL,
        szCmdline,     // command line
        NULL,          // process security attributes
        NULL,          // primary thread security attributes
        TRUE,          // handles are inherited
        CREATE_NEW_PROCESS_GROUP,             // creation flags
        NULL,          // use parent's environment
        NULL,          // use parent's current directory
        &siStartInfo,  // STARTUPINFO pointer
        &piProcInfo);  // receives PROCESS_INFORMATION

    CloseHandle(piProcInfo.hThread);

    WriteFile(g_hChildStd_IN_Wr, decompString.c_str(), decompString.size(), 0, NULL);

    CloseHandle(g_hChildStd_IN_Wr);

    DWORD dwRead;

    CloseHandle(siStartInfo.hStdOutput);

    DWORD test = WaitForSingleObject(piProcInfo.hProcess, timeout_ == 0 ? INFINITE : timeout_);

    test = GetLastError();

    DWORD code;

    GetExitCodeProcess(piProcInfo.hProcess, &code);

    if (code == STILL_ACTIVE)
    {
        TerminateProcess(piProcInfo.hProcess, 0);
    }

    CloseHandle(piProcInfo.hProcess);

    CHAR chBuf[1024];

    for (;;)
    {
        bSuccess = ReadFile(g_hChildStd_OUT_Rd, chBuf, 1024, &dwRead, NULL);
        if (!bSuccess || dwRead == 0) break;
        decompString.append(chBuf,dwRead);
    }

    return decompString;

}

#elif __GNUC__

#include <unistd.h>
#include <sys/wait.h>
#include <poll.h>
#include <wordexp.h>

std::string htd::ExternalPipeTreeDecompositionAlgorithm::Implementation::getDecomp(const htd::IMultiHypergraph & graph) const
{
    std::string graphString = convert(graph);

    std::string decompString;

    int inpipefd[2];

    int outpipefd[2];

    char buf[256];

    if (pipe(inpipefd) < 0)
    {
        return std::string("");
    };

    if (pipe(outpipefd) < 0)
    {
        return std::string("");
    };

    pid_t pid_solver = fork();

    if (pid_solver < 0)
    {
        return std::string("");
    }

    if (pid_solver == 0)
    {
        if (dup2(outpipefd[0], STDIN_FILENO) < 0)
        {
            _Exit(0);
        };

        close(outpipefd[1]);

        if (dup2(inpipefd[1], STDOUT_FILENO) < 0)
        {
            _Exit(0);
        };

        close(inpipefd[0]);

        wordexp_t p;

        wordexp(cmd_.c_str(), &p, 0);

        char ** cmd = p.we_wordv;

        chdir(dir_.c_str());

        execvp(cmd[0], cmd);

        _Exit(0);
    }

    //close unused pipe ends
    close(outpipefd[0]);

    close(inpipefd[1]);

    size_t graphLen = graphString.length();

    if (write(outpipefd[1], graphString.c_str(), graphLen) < (long)graphLen)
    {
        return std::string("");
    }

    close(outpipefd[1]);


    pid_t pid_watcher = 0;

    if (timeout_ > 0)
    {
        pid_watcher = fork();

        if (pid_watcher == 0)
        {
            usleep(timeout_ * 1000);

            kill(pid_solver, SIGTERM);

            usleep(1000000);

            kill(pid_solver, SIGKILL);
        }
    }

    ssize_t outputLen = 256;

    while (outputLen > 0)
    {
        outputLen = read(inpipefd[0], buf, 256);

        decompString.append(buf, outputLen);
    }

    close(inpipefd[0]);

    if (timeout_ > 0 && pid_watcher != 0)
    {
        kill(pid_watcher, SIGKILL);
    }

    return decompString;
}

#endif
#endif /* HTD_EXTERNALPIPETREEDECOMPOSITIONALGORITHM_CPP */

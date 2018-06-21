#ifndef HTD_IEXTERNALTREEDECOMPOSITIONALGORITHM_HPP
#define HTD_IEXTERNALTREEDECOMPOSITIONALGORITHM_HPP

#include <htd/Globals.hpp>
#include <htd/ITreeDecompositionAlgorithm.hpp>
#include <htd/IOrderingAlgorithm.hpp>
#include <htd/IPreprocessedGraph.hpp>

#include <utility>

namespace htd
{
    /**
     * Inteface for a tree decomposition algorithm which calls an external solver.
     */
    class IExternalTreeDecompositionAlgorithm : public virtual htd::ITreeDecompositionAlgorithm
    {
        public:
            /**
             * Setter for the command which is used to call the external solver.
             *
             * @param cmd   The command used to call the external solver.
             */
            virtual void setCommand(const std::string & cmd) = 0;

            /**
             * Getter for the command which is used to call the external solver.
             *
             * @return  The command used to call the external solver.
             */
            virtual std::string getCommand() = 0;

            /**
             * Setter for the timeout of the external solver.
             *
             * @param timeout   The timeout for the external solver in milliseconds.
             */
            virtual void setTimeout(const unsigned int & timeout) = 0;

            /**
             * Getter for the timeout of the external solver.
             *
             * @return  The timeout for the external solver in milliseconds.
             */
            virtual unsigned int getTimeout() = 0;

            /**
             * Setter for the directory in which the solver gets executed.
             *
             * @param dir   The directory in which the solver gets executed.
             */
            virtual void setDirectory(const std::string dir) = 0;

            /**
             * Getter for the directory in which the solver gets executed.
             *
             * @return   The directory in which the solver gets executed.
             */
            virtual std::string getDirectory() = 0;

            virtual IExternalTreeDecompositionAlgorithm * clone() const HTD_OVERRIDE = 0;
    };
}

#endif //HTD_IEXTERNALTREEDECOMPOSITIONALGORITHM_HPP

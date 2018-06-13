#ifndef HTD_IO_IGRAPHEXPORTER_HPP
#define HTD_IO_IGRAPHEXPORTER_HPP

#include <htd/IMultiHypergraph.hpp>
#include <htd/NamedMultiHypergraph.hpp>
#include <htd/ITreeDecomposition.hpp>

#include <iostream>

namespace htd_io
{
    /**
     * Interface for algorithms which can be used to export graphs to streams.
     */
    class IGraphExporter
    {
        public:
            virtual ~IGraphExporter() = 0;

            /**
             *  Write a graph to a given stream.
             *
             *  @param[in] graph            The graph which shall be exported.
             *  @param[out] outputStream    The output stream to which the information shall be written.
             */
            virtual void write(const htd::IMultiHypergraph & graph, std::ostream & outputStream) const = 0;
    };

    inline htd_io::IGraphExporter::~IGraphExporter() { }
}

#endif /* HTD_IO_IGRAPHEXPORTER_HPP */

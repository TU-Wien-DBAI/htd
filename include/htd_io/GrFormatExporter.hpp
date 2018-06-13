#ifndef HTD_IO_GRFORMATEXPORTER_HPP
#define HTD_IO_GRFORMATEXPORTER_HPP

#include <htd_io/PreprocessorDefinitions.hpp>

#include <htd_io/IGraphExporter.hpp>

#include <iostream>

namespace htd_io
{
    /**
     *  Exporter which allows to export graphs in the graph format 'gr'.
     */
    class GrFormatExporter : public htd_io::IGraphExporter
    {
        public:
            HTD_IO_API GrFormatExporter();

            HTD_IO_API ~GrFormatExporter() final;

            HTD_IO_API void write(const htd::IMultiHypergraph & graph, std::ostream & outputStream) const HTD_OVERRIDE;
    };
}

#endif /* HTD_IO_GRFORMATEXPORTER_HPP */

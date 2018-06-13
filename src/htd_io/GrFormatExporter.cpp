#ifndef HTD_IO_GRFORMATEXPORTER_CPP
#define HTD_IO_GRFORMATEXPORTER_CPP

#include <htd_io/GrFormatExporter.hpp>

#include <sstream>
#include <unordered_map>

htd_io::GrFormatExporter::GrFormatExporter() = default;

htd_io::GrFormatExporter::~GrFormatExporter() = default;

void htd_io::GrFormatExporter::write(const htd::IMultiHypergraph & graph, std::ostream & outputStream) const
{
    std::unordered_map<htd::vertex_t, std::size_t> indices;

    outputStream << "p tw " << std::to_string(graph.vertexCount()) << " " << std::to_string(graph.edgeCount()) << "\n";

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

            outputStream << std::to_string(vertex1) << " " << std::to_string(vertex2) << "\n";

            ++it;
        }
    }
}

#endif /* HTD_IO_GRFORMATEXPORTER_CPP */

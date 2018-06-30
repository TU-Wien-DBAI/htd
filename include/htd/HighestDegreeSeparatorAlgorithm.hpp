#ifndef HTD_HTD_HIGHESTDEGREESEPARATORALGORITHM_HPP
#define HTD_HTD_HIGHESTDEGREESEPARATORALGORITHM_HPP

#include <htd/IGraphSeparatorAlgorithm.hpp>
#include <htd/MultiGraph.hpp>

namespace htd
{
	/*
	*	Implementation of the IGraphSeparator algorithm which computes miniaml separating vertex set
	*/
	class HighestDegreeSeparatorAlgorithm : public htd::IGraphSeparatorAlgorithm
	{
	public:
		/*
		*
		*
		*/
		HTD_API HighestDegreeSeparatorAlgorithm(const htd::LibraryInstance * const manager);


		/**
		*  Destructor of a minimal separator algorithm.
		*/
		HTD_API virtual ~HighestDegreeSeparatorAlgorithm();

		HTD_API std::vector<htd::vertex_t> * computeSeparator(const htd::IGraphStructure & graph) const HTD_OVERRIDE;

		HTD_API std::vector<htd::vertex_t> * computeSeparator(htd::MultiGraph * graph) const;

		HTD_API const htd::LibraryInstance * managementInstance(void) const HTD_NOEXCEPT HTD_OVERRIDE;

		HTD_API void setManagementInstance(const htd::LibraryInstance * const manager) HTD_OVERRIDE;

		HTD_API HighestDegreeSeparatorAlgorithm * clone(void) const HTD_OVERRIDE;

	private:
		struct Implementation;

		std::unique_ptr<Implementation> implementation_;
	};

}

#endif /* HTD_HTD_HIGHESTDEGREESEPARATORALGORITHM_HPP */
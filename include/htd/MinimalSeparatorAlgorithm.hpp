#ifndef HTD_HTD_MINIMALSEPARATORALGORITHM_HPP
#define HTD_HTD_MINIMALSEPARATORALGORITHM_HPP

#include <htd/IGraphSeparatorAlgorithm.hpp>

namespace htd
{
	/*
	*	Implementation of the IGraphSeparator algorithm which computes miniaml separating vertex set
	*/
	class MinimalSeparatorAlgorithm : public htd::IGraphSeparatorAlgorithm
	{
	public:
		/*
		*
		*
		*/
		HTD_API MinimalSeparatorAlgorithm(const htd::LibraryInstance * const manager);


		/**
		*  Destructor of a minimal separator algorithm.
		*/
		HTD_API virtual ~MinimalSeparatorAlgorithm();

		HTD_API std::vector<htd::vertex_t> * computeSeparator(const htd::IGraphStructure & graph) const HTD_OVERRIDE;

		HTD_API const htd::LibraryInstance * managementInstance(void) const HTD_NOEXCEPT HTD_OVERRIDE;

		HTD_API void setManagementInstance(const htd::LibraryInstance * const manager) HTD_OVERRIDE;

		HTD_API MinimalSeparatorAlgorithm * clone(void) const HTD_OVERRIDE;

	private:
		struct Implementation;

		std::unique_ptr<Implementation> implementation_;
	};

}

#endif /* HTD_HTD_MINIMALSEPARATORALGORITHM_HPP */
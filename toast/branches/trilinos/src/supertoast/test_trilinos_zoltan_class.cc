#include "test_trilinos_zoltan_class.h"


template<class MT>
void TCompRowMatrixTrilinosTest<MT>::tmp (
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    const Teuchos::RCP<typename matrix_type::node_type>& node)
{
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::arcp;

    const Tpetra::global_size_t numGlobalElements = (const Tpetra::global_size_t)100;
    const size_t numMyElements = 10;
    const int indexBase = 0;

    std::vector<int> myNode(numGlobalElements);
    ArrayView<int> myNodeView(myNode);

    ArrayRCP<size_t> NumNz = arcp<size_t> (numMyElements);

    RCP<matrix_type> A = rcp (new matrix_type (map, NumNz,
					       Tpetra::StaticProfile));
}

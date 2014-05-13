#ifndef __TEST_TRILINOS_ZOLTAN_CLASS_H
#define __TEST_TRILINOS_ZOLTAN_CLASS_H

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>
#include <Ifpack2_Factory.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <zoltan.h>
#include <iostream>
#include "stoastlib.h"
#include "mesh_zoltan.h"

template<class MT> class TCompRowMatrixTrilinosTest
{
public:
    typedef Kokkos::DefaultNode::DefaultNodeType node_type;
    typedef Tpetra::CrsMatrix<double, int, int, node_type> matrix_type;
    typedef Tpetra::Map<int, int, node_type> map_type;
  
    TCompRowMatrixTrilinosTest() {}

    void tmp (const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
	      const Teuchos::RCP<typename matrix_type::node_type>& node);

private:
    Teuchos::RCP<const map_type> map;
};

#endif // !__TEST_TRILINOS_ZOLTAN_CLASS_H

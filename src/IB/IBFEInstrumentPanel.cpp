
/*
 * File:   IBFEInstrumentPanel.cpp
 * Author: cpuelz
 *
 * Created on April 12, 2018, 11:06 AM
 */

#include <algorithm>
#include <fstream>
#include <limits>
#include <map>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string>
#include <utility>
#include <vector>

#include "BasePatchLevel.h"
#include "Box.h"
#include "BoxArray.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "Eigen/Geometry" // IWYU pragma: keep
#include "IBAMR_config.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SideData.h"
#include "SideIndex.h"
#include "boost/array.hpp"
#include "boost/multi_array.hpp"
#include "ibamr/IBFEInstrumentPanel.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LMesh.h"
#include "ibtk/LNode.h"
#include "ibtk/ibtk_utilities.h"
#include "petscvec.h"
#include "tbox/Database.h"

#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

#include "libmesh/boundary_info.h"
#include "libmesh/dense_vector.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/face_tri3.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_function.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/point.h"
#include "libmesh/serial_mesh.h"

#include "ibamr/IBFEMethod.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/IndexUtilities.h"

#include "BoxArray.h"

#include <sstream>

// static functions
namespace
{
double
linear_interp(const Vector& X,
              const Index<NDIM>& i_cell,
              const Vector& X_cell,
              const CellData<NDIM, double>& v,
              const Index<NDIM>& /*patch_lower*/,
              const Index<NDIM>& /*patch_upper*/,
              const double* const /*x_lower*/,
              const double* const /*x_upper*/,
              const double* const dx)
{
    boost::array<bool, NDIM> is_lower;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        is_lower[d] = X[d] < X_cell[d];
    }
    double U = 0.0;
#if (NDIM == 3)
    for (int i_shift2 = (is_lower[2] ? -1 : 0); i_shift2 <= (is_lower[2] ? 0 : 1); ++i_shift2)
    {
#endif
        for (int i_shift1 = (is_lower[1] ? -1 : 0); i_shift1 <= (is_lower[1] ? 0 : 1); ++i_shift1)
        {
            for (int i_shift0 = (is_lower[0] ? -1 : 0); i_shift0 <= (is_lower[0] ? 0 : 1); ++i_shift0)
            {
                const Vector X_center(X_cell[0] + static_cast<double>(i_shift0) * dx[0],
                                      X_cell[1] + static_cast<double>(i_shift1) * dx[1]
#if (NDIM == 3)
                                      ,
                                      X_cell[2] + static_cast<double>(i_shift2) * dx[2]
#endif
                                      );
                const double wgt =
                    (((X[0] < X_center[0] ? X[0] - (X_center[0] - dx[0]) : (X_center[0] + dx[0]) - X[0]) / dx[0]) *
                     ((X[1] < X_center[1] ? X[1] - (X_center[1] - dx[1]) : (X_center[1] + dx[1]) - X[1]) / dx[1])
#if (NDIM == 3)
                     *
                     ((X[2] < X_center[2] ? X[2] - (X_center[2] - dx[2]) : (X_center[2] + dx[2]) - X[2]) / dx[2])
#endif
                         );
                const Index<NDIM> i(i_shift0 + i_cell(0),
                                    i_shift1 + i_cell(1)
#if (NDIM == 3)
                                        ,
                                    i_shift2 + i_cell(2)
#endif
                                        );
                const CellIndex<NDIM> i_c(i);
                U += v(i_c) * wgt;
            }
        }
#if (NDIM == 3)
    }
#endif
    return U;
} // linear_interp

template <int N>
Eigen::Matrix<double, N, 1>
linear_interp(const Vector& X,
              const Index<NDIM>& i_cell,
              const Vector& X_cell,
              const CellData<NDIM, double>& v,
              const Index<NDIM>& /*patch_lower*/,
              const Index<NDIM>& /*patch_upper*/,
              const double* const /*x_lower*/,
              const double* const /*x_upper*/,
              const double* const dx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(v.getDepth() == N);
#endif
    boost::array<bool, NDIM> is_lower;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        is_lower[d] = X[d] < X_cell[d];
    }
    Eigen::Matrix<double, N, 1> U(Eigen::Matrix<double, N, 1>::Zero());
#if (NDIM == 3)
    for (int i_shift2 = (is_lower[2] ? -1 : 0); i_shift2 <= (is_lower[2] ? 0 : 1); ++i_shift2)
    {
#endif
        for (int i_shift1 = (is_lower[1] ? -1 : 0); i_shift1 <= (is_lower[1] ? 0 : 1); ++i_shift1)
        {
            for (int i_shift0 = (is_lower[0] ? -1 : 0); i_shift0 <= (is_lower[0] ? 0 : 1); ++i_shift0)
            {
                const Vector X_center(X_cell[0] + static_cast<double>(i_shift0) * dx[0],
                                      X_cell[1] + static_cast<double>(i_shift1) * dx[1]
#if (NDIM == 3)
                                      ,
                                      X_cell[2] + static_cast<double>(i_shift2) * dx[2]
#endif
                                      );
                const double wgt =
                    (((X[0] < X_center[0] ? X[0] - (X_center[0] - dx[0]) : (X_center[0] + dx[0]) - X[0]) / dx[0]) *
                     ((X[1] < X_center[1] ? X[1] - (X_center[1] - dx[1]) : (X_center[1] + dx[1]) - X[1]) / dx[1])
#if (NDIM == 3)
                     *
                     ((X[2] < X_center[2] ? X[2] - (X_center[2] - dx[2]) : (X_center[2] + dx[2]) - X[2]) / dx[2])
#endif
                         );
                const Index<NDIM> i(i_shift0 + i_cell(0),
                                    i_shift1 + i_cell(1)
#if (NDIM == 3)
                                        ,
                                    i_shift2 + i_cell(2)
#endif
                                        );
                const CellIndex<NDIM> i_c(i);
                for (int k = 0; k < N; ++k)
                {
                    U[k] += v(i_c, k) * wgt;
                }
            }
        }
#if (NDIM == 3)
    }
#endif
    return U;
} // linear_interp

Vector
linear_interp(const Vector& X,
              const Index<NDIM>& i_cell,
              const Vector& X_cell,
              const SideData<NDIM, double>& v,
              const Index<NDIM>& /*patch_lower*/,
              const Index<NDIM>& /*patch_upper*/,
              const double* const /*x_lower*/,
              const double* const /*x_upper*/,
              const double* const dx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(v.getDepth() == 1);
#endif
    Vector U(Vector::Zero());
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        boost::array<bool, NDIM> is_lower;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            if (d == axis)
            {
                is_lower[d] = false;
            }
            else
            {
                is_lower[d] = X[d] < X_cell[d];
            }
        }
#if (NDIM == 3)
        for (int i_shift2 = (is_lower[2] ? -1 : 0); i_shift2 <= (is_lower[2] ? 0 : 1); ++i_shift2)
        {
#endif
            for (int i_shift1 = (is_lower[1] ? -1 : 0); i_shift1 <= (is_lower[1] ? 0 : 1); ++i_shift1)
            {
                for (int i_shift0 = (is_lower[0] ? -1 : 0); i_shift0 <= (is_lower[0] ? 0 : 1); ++i_shift0)
                {
                    const Vector X_side(X_cell[0] + (static_cast<double>(i_shift0) + (axis == 0 ? -0.5 : 0.0)) * dx[0],
                                        X_cell[1] + (static_cast<double>(i_shift1) + (axis == 1 ? -0.5 : 0.0)) * dx[1]
#if (NDIM == 3)
                                        ,
                                        X_cell[2] + (static_cast<double>(i_shift2) + (axis == 2 ? -0.5 : 0.0)) * dx[2]
#endif
                                        );
                    const double wgt =
                        (((X[0] < X_side[0] ? X[0] - (X_side[0] - dx[0]) : (X_side[0] + dx[0]) - X[0]) / dx[0]) *
                         ((X[1] < X_side[1] ? X[1] - (X_side[1] - dx[1]) : (X_side[1] + dx[1]) - X[1]) / dx[1])
#if (NDIM == 3)
                         *
                         ((X[2] < X_side[2] ? X[2] - (X_side[2] - dx[2]) : (X_side[2] + dx[2]) - X[2]) / dx[2])
#endif
                             );
                    const Index<NDIM> i(i_shift0 + i_cell(0),
                                        i_shift1 + i_cell(1)
#if (NDIM == 3)
                                            ,
                                        i_shift2 + i_cell(2)
#endif
                                            );
                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                    U[axis] += v(i_s) * wgt;
                }
            }
#if (NDIM == 3)
        }
#endif
    }
    return U;
} // linear_interp
}

//**********************************************
// constructor
//**********************************************
IBFEInstrumentPanel::IBFEInstrumentPanel(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db, const int part)
    : d_num_meters(0),
      d_part(part),
      d_nodes(),
      d_node_dof_IDs(),
      d_quad_order(),
      d_num_quad_points(),
      d_num_nodes(),
      d_U_dof_idx(),
      d_dX_dof_idx(),
      d_nodeset_IDs_for_meters(),
      d_meter_meshes(),
      d_meter_systems(),
      d_meter_mesh_names(),
      d_quad_point_map(),
      d_flow_values(),
      d_mean_pressure_values(),
      d_flux_stream(),
      d_mean_pressure_stream(),
      d_plot_directory_name(NDIM == 2 ? "viz_inst2d" : "viz_inst3d")
{
    // get input data
    IBFEInstrumentPanel::getFromInput(input_db);

    // make plot directory
    Utilities::recursiveMkdir(d_plot_directory_name);

    // set up file streams
    if (SAMRAI_MPI::getRank() == 0)
    {
        std::ostringstream press_output;
        std::ostringstream flux_output;
        press_output << d_plot_directory_name << "/"
                     << ""
                     << "mean_pressure.dat";
        d_mean_pressure_stream.open(press_output.str().c_str());
        flux_output << d_plot_directory_name << "/"
                    << ""
                    << "flux.dat";
        d_flux_stream.open(flux_output.str().c_str());
        d_mean_pressure_stream.precision(15);
        d_flux_stream.precision(15);
    }
} // constructor

//**********************************************
// destructor
//**********************************************
IBFEInstrumentPanel::~IBFEInstrumentPanel()
{
    // Close the log file stream.
    if (SAMRAI_MPI::getRank() == 0)
    {
        d_mean_pressure_stream.close();
        d_flux_stream.close();
    }
    // delete vectors of pointers
    for (int ii = 0; ii < d_num_meters; ++ii)
    {
        delete d_exodus_io[ii];
        delete d_meter_meshes[ii];
        delete d_meter_systems[ii];
    }
} // destructor

//**********************************************
// initialize data which is independent of the Cartesian hierarchy
//**********************************************
void
IBFEInstrumentPanel::initializeHierarchyIndependentData(IBAMR::IBFEMethod* ib_method_ops)
{
    // get relevant things for corresponding part
    const FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager(d_part);
    const EquationSystems* equation_systems = fe_data_manager->getEquationSystems();
    const MeshBase* mesh = &equation_systems->get_mesh();
    const BoundaryInfo& boundary_info = *mesh->boundary_info;
    const libMesh::Parallel::Communicator& comm_in = mesh->comm();

    // get equation systems from the mesh we will need.
    const System& dX_system = equation_systems->get_system(IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
    const unsigned int dX_sys_num = dX_system.number();
    const System& U_system = equation_systems->get_system(IBFEMethod::VELOCITY_SYSTEM_NAME);
    const unsigned int U_sys_num = U_system.number();

    // some local variables
    std::vector<dof_id_type> nodes;
    std::vector<boundary_id_type> bcs;
    std::vector<std::vector<dof_id_type> > temp_node_dof_IDs;
    std::vector<std::vector<libMesh::Point> > temp_nodes;
    std::vector<libMesh::Point> meter_centroids;
    boundary_info.build_node_list(nodes, bcs);

    // check to make sure there are node sets to work with
    TBOX_ASSERT(nodes.size() > 0 && bcs.size() > 0);
    
    // resize members and local variables
    d_num_meters = d_nodeset_IDs_for_meters.size();
    d_num_quad_points.resize(d_num_meters);
    d_U_dof_idx.resize(d_num_meters);
    d_dX_dof_idx.resize(d_num_meters);
    d_node_dof_IDs.resize(d_num_meters);
    d_nodes.resize(d_num_meters);
    d_num_nodes.resize(d_num_meters);
    d_mean_pressure_values.resize(d_num_meters);
    d_flow_values.resize(d_num_meters);
    temp_node_dof_IDs.resize(d_num_meters);
    temp_nodes.resize(d_num_meters);
    meter_centroids.resize(d_num_meters);

    // populate temp vectors
    for (int ii = 0; ii < nodes.size(); ++ii)
    {
        for (int jj = 0; jj < d_nodeset_IDs_for_meters.size(); ++jj)
        {
            if (d_nodeset_IDs_for_meters[jj] == bcs[ii])
            {
                temp_node_dof_IDs[jj].push_back(nodes[ii]);
                const Node* node = &mesh->node_ref(nodes[ii]);
                temp_nodes[jj].push_back(*node);
                meter_centroids[jj] += *node;
            }
        }
    }

    // loop over meters and sort the nodes
    for (int jj = 0; jj < d_num_meters; ++jj)
    {
        // finish computing centroid
        meter_centroids[jj] /= static_cast<double>(temp_nodes[jj].size());
        // make sure nodes are sorted with a particular orientation
        libMesh::Point temp_point;
        dof_id_type temp_dof_id;
        double max_dist = std::numeric_limits<double>::max();
        d_nodes[jj].push_back(temp_nodes[jj][0]);
        d_node_dof_IDs[jj].push_back(temp_node_dof_IDs[jj][0]);
        max_dist = std::numeric_limits<double>::max();
        for (int kk = 1; kk < temp_nodes[jj].size(); ++kk)
        {
            for (int ll = 0; ll < temp_nodes[jj].size(); ++ll)
            {
                // here we find the closest node to the previous one
                // with some orientation
                libMesh::Point dist = temp_nodes[jj][ll] - d_nodes[jj][kk - 1];

                if (dist.norm() < max_dist)
                {
                    // make sure we haven't already added this node
                    bool added = false;
                    for (int ii = 1; ii < kk + 1; ++ii)
                    {
                        if (temp_nodes[jj][ll] == d_nodes[jj][ii - 1]) added = true;
                    }
                    if (!added)
                    {
                        temp_point = temp_nodes[jj][ll];
                        temp_dof_id = temp_node_dof_IDs[jj][ll];
                        max_dist = dist.norm();
                    }
                }
            }
            d_nodes[jj].push_back(temp_point);
            d_node_dof_IDs[jj].push_back(temp_dof_id);
            max_dist = std::numeric_limits<double>::max();
        }
    } // loop over meters

    // initialize meshes and number of nodes
    for (int jj = 0; jj < d_num_meters; ++jj)
    {
        d_meter_meshes.push_back(new SerialMesh(comm_in, NDIM));
        std::ostringstream id;
        id << d_nodeset_IDs_for_meters[jj];
        d_meter_mesh_names.push_back("meter_mesh_" + id.str());
        d_num_nodes[jj] = d_nodes[jj].size();
    }

    // build the meshes
    for (int ii = 0; ii < d_num_meters; ++ii)
    {
        d_meter_meshes[ii]->set_spatial_dimension(NDIM);
        d_meter_meshes[ii]->set_mesh_dimension(NDIM - 1);
        d_meter_meshes[ii]->reserve_nodes(d_num_nodes[ii]);
        d_meter_meshes[ii]->reserve_elem(d_num_nodes[ii] - 2);

        for (int jj = 0; jj < d_num_nodes[ii]; ++jj)
        {
            d_meter_meshes[ii]->add_point(d_nodes[ii][jj], jj);
        }

        for (int jj = 0; jj < d_num_nodes[ii] - 2; ++jj)
        {
            Elem* elem = new Tri3;
            elem->set_id(jj);
            elem = d_meter_meshes[ii]->add_elem(elem);
            elem->set_node(0) = d_meter_meshes[ii]->node_ptr(0);
            elem->set_node(1) = d_meter_meshes[ii]->node_ptr(jj + 1);
            elem->set_node(2) = d_meter_meshes[ii]->node_ptr(jj + 2);
        }
        d_meter_meshes[ii]->prepare_for_use();
        d_exodus_io.push_back(new ExodusII_IO(*d_meter_meshes[ii]));
    } // loop over meters

    // initialize meter mesh equation systems, for both velocity and displacement
    for (int jj = 0; jj < d_num_meters; ++jj)
    {
        d_meter_systems.push_back(new EquationSystems(*d_meter_meshes[jj]));
        LinearImplicitSystem& velocity_sys =
            d_meter_systems[jj]->add_system<LinearImplicitSystem>(IBFEMethod::VELOCITY_SYSTEM_NAME);
        velocity_sys.add_variable("U_0", static_cast<Order>(1), LAGRANGE);
        velocity_sys.add_variable("U_1", static_cast<Order>(1), LAGRANGE);
        velocity_sys.add_variable("U_2", static_cast<Order>(1), LAGRANGE);
        LinearImplicitSystem& displacement_sys =
            d_meter_systems[jj]->add_system<LinearImplicitSystem>(IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
        displacement_sys.add_variable("dX_0", static_cast<Order>(1), LAGRANGE);
        displacement_sys.add_variable("dX_1", static_cast<Order>(1), LAGRANGE);
        displacement_sys.add_variable("dX_2", static_cast<Order>(1), LAGRANGE);

        // we use these serialized vectors to store the DOFs for the systems on all processes.
        velocity_sys.add_vector("serial solution", false, libMesh::SERIAL);
        displacement_sys.add_vector("serial solution", false, libMesh::SERIAL);

        d_meter_systems[jj]->init();
    }

    // store the number of quadrature points for each meter mesh
    for (int jj = 0; jj < d_num_meters; ++jj)
    {
        const LinearImplicitSystem& displacement_sys =
            d_meter_systems[jj]->get_system<LinearImplicitSystem>(IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
        FEType fe_type = displacement_sys.variable_type(0);
        UniquePtr<FEBase> fe_elem(FEBase::build(NDIM - 1, fe_type));
        QGauss qrule(NDIM - 1, d_quad_order);
        fe_elem->attach_quadrature_rule(&qrule);
        const std::vector<libMesh::Point>& qp_points = fe_elem->get_xyz();
        MeshBase::const_element_iterator el = d_meter_meshes[jj]->active_elements_begin();
        const MeshBase::const_element_iterator end_el = d_meter_meshes[jj]->active_elements_end();
        for (; el != end_el; ++el)
        {
            const Elem* elem = *el;
            fe_elem->reinit(elem);
            d_num_quad_points[jj] += qp_points.size();
        }
    }

    // store dof indices for the velocity and displacement systems that we will use later
    for (int jj = 0; jj < d_num_meters; ++jj)
    {
        for (int ii = 0; ii < d_num_nodes[jj]; ++ii)
        {
            const Node* node = &mesh->node_ref(d_node_dof_IDs[jj][ii]);
            std::vector<dof_id_type> dX_dof_index;
            std::vector<dof_id_type> U_dof_index;
            for (int d = 0; d < NDIM; ++d)
            {
                dX_dof_index.push_back(node->dof_number(dX_sys_num, d, 0));
                U_dof_index.push_back(node->dof_number(U_sys_num, d, 0));
            }
            d_dX_dof_idx[jj].push_back(dX_dof_index);
            d_U_dof_idx[jj].push_back(U_dof_index);
        }
    }
    d_initialized = true;

} // initializeHierarchyIndependentData

//**********************************************
// initialize data which depends on the Cartesian hierarchy
//**********************************************
void
IBFEInstrumentPanel::initializeHierarchyDependentData(IBAMR::IBFEMethod* ib_method_ops,
                                                      const Pointer<PatchHierarchy<NDIM> > hierarchy)
{
    if (!d_initialized)
    {
        initializeHierarchyIndependentData(ib_method_ops);
    }
    if (d_num_meters == 0) return;

    // loop over meters and update system data
    for (int jj = 0; jj < d_num_meters; ++jj)
    {
        // update FE system data for meter_mesh
        updateSystemData(ib_method_ops, jj);
    }

    // get info about levels in AMR mesh
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    // Determine the finest grid spacing in the Cartesian grid hierarchy.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const domainXLower = grid_geom->getXLower();
    const double* const domainXUpper = grid_geom->getXUpper();
    const double* const dx_coarsest = grid_geom->getDx();
    TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
    const Box<NDIM> domain_box = grid_geom->getPhysicalDomain()[0];

    // reset the quad point maps
    d_quad_point_map.clear();
    d_quad_point_map.resize(finest_ln + 1);

    // loop over levels and assign each quadrature point to one level
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        const IntVector<NDIM>& ratio = level->getRatio();
        const Box<NDIM> domain_box_level = Box<NDIM>::refine(domain_box, ratio);
        const Index<NDIM>& domain_box_level_lower = domain_box_level.lower();
        const Index<NDIM>& domain_box_level_upper = domain_box_level.upper();
        boost::array<double, NDIM> dx;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dx[d] = dx_coarsest[d] / static_cast<double>(ratio(d));
        }

        Pointer<PatchLevel<NDIM> > finer_level =
            (ln < finest_ln ? hierarchy->getPatchLevel(ln + 1) : Pointer<BasePatchLevel<NDIM> >(NULL));
        const IntVector<NDIM>& finer_ratio = (ln < finest_ln ? finer_level->getRatio() : IntVector<NDIM>(1));
        const Box<NDIM> finer_domain_box_level = Box<NDIM>::refine(domain_box, finer_ratio);
        const Index<NDIM>& finer_domain_box_level_lower = finer_domain_box_level.lower();
        const Index<NDIM>& finer_domain_box_level_upper = finer_domain_box_level.upper();
        boost::array<double, NDIM> finer_dx;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            finer_dx[d] = dx_coarsest[d] / static_cast<double>(finer_ratio(d));
        }

        for (unsigned int jj = 0; jj < d_num_meters; ++jj)
        {
            const LinearImplicitSystem& displacement_sys =
                d_meter_systems[jj]->get_system<LinearImplicitSystem>(IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
            const NumericVector<double>& displacement_coords = displacement_sys.get_vector("serial solution");
            const DofMap& dof_map = displacement_sys.get_dof_map();
            FEType fe_type = displacement_sys.variable_type(0);

            // set up FE objects
            UniquePtr<FEBase> fe_elem(FEBase::build(NDIM - 1, fe_type));
            QGauss qrule(NDIM - 1, d_quad_order);
            fe_elem->attach_quadrature_rule(&qrule);

            //  for evaluating the displacement system
            const std::vector<Real>& JxW = fe_elem->get_JxW();
            const std::vector<std::vector<Real> >& phi = fe_elem->get_phi();
            const std::vector<libMesh::Point>& qp_points = fe_elem->get_xyz();
            std::vector<dof_id_type> dof_indices;
            DenseMatrix<double> disp_coords;

            // loop over ALL elements in meter mesh, not just the local ones on this process!!
            MeshBase::const_element_iterator el = d_meter_meshes[jj]->active_elements_begin();
            const MeshBase::const_element_iterator end_el = d_meter_meshes[jj]->active_elements_end();
            for (; el != end_el; ++el)
            {
                const Elem* elem = *el;
                fe_elem->reinit(elem);

                // get dofs for displacement system and store in dense matrix
                disp_coords.resize(NDIM, phi.size());
                for (int d = 0; d < NDIM; ++d) // here d is the "variable number"
                {
                    dof_map.dof_indices(elem, dof_indices, d);
                    for (int nn = 0; nn < dof_indices.size(); ++nn)
                    {
                        disp_coords(d, nn) = displacement_coords(dof_indices[nn]);
                    }
                }

                // compute normal vector to element
                const libMesh::Point foo1 = *elem->node_ptr(1) - *elem->node_ptr(0);
                const libMesh::Point foo2 = *elem->node_ptr(2) - *elem->node_ptr(1);
                libMesh::Point foo3 = foo1.cross(foo2).unit();
                Vector normal;
                for (int d = 0; d < NDIM; ++d) normal[d] = foo3(d);

                // loop over quadrature points, compute their physical locations
                // after displacement, and stores their indices.
                for (int qp = 0; qp < qp_points.size(); ++qp)
                {
                    Vector qp_temp;
                    double disp_comp = 0.0;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        disp_comp = 0.0;
                        for (int nn = 0; nn < phi.size(); ++nn)
                        {
                            disp_comp += disp_coords(d, nn) * phi[nn][qp];
                        }
                        // calculating physical location of the quadrature point
                        qp_temp[d] = qp_points[qp](d) + disp_comp;
                    }

                    const Index<NDIM> i = IndexUtilities::getCellIndex(&qp_temp[0],
                                                                       domainXLower,
                                                                       domainXUpper,
                                                                       dx.data(),
                                                                       domain_box_level_lower,
                                                                       domain_box_level_upper);

                    const Index<NDIM> finer_i = IndexUtilities::getCellIndex(&qp_temp[0],
                                                                             domainXLower,
                                                                             domainXUpper,
                                                                             finer_dx.data(),
                                                                             finer_domain_box_level_lower,
                                                                             finer_domain_box_level_upper);

                    if (level->getBoxes().contains(i) &&
                        (ln == finest_ln || !finer_level->getBoxes().contains(finer_i)))
                    {
                        QuadPointStruct q;
                        q.meter_num = jj;
                        q.qp_xyz_current = qp_temp;
                        q.JxW = JxW[qp];
                        q.normal = normal;
                        d_quad_point_map[ln].insert(std::make_pair(i, q));
                    }
                }
            }
        }
    }
} // initializeHierarchyDependentData

//**********************************************
// read instrument data
//**********************************************
void
IBFEInstrumentPanel::readInstrumentData(const int U_data_idx,
                                        const int P_data_idx,
                                        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                        const double data_time)
{
    if (d_num_meters == 0) return;

    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    // Reset the instrument values.
    std::fill(d_flow_values.begin(), d_flow_values.end(), 0.0);
    std::fill(d_mean_pressure_values.begin(), d_mean_pressure_values.end(), 0.0);
    std::vector<double> A(d_num_meters, 0.0);

    int count_qp = 0;

    int count_qp_2 = 0;

    // compute flow and mean pressure on mesh meters
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        count_qp += d_quad_point_map[ln].size();

        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Index<NDIM>& patch_lower = patch_box.lower();
            const Index<NDIM>& patch_upper = patch_box.upper();

            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const x_lower = pgeom->getXLower();
            const double* const x_upper = pgeom->getXUpper();
            const double* const dx = pgeom->getDx();

            Pointer<CellData<NDIM, double> > U_cc_data = patch->getPatchData(U_data_idx);
            Pointer<SideData<NDIM, double> > U_sc_data = patch->getPatchData(U_data_idx);
            Pointer<CellData<NDIM, double> > P_cc_data = patch->getPatchData(P_data_idx);

            for (Box<NDIM>::Iterator b(patch_box); b; b++)
            {
                const Index<NDIM>& i = b();
                std::pair<QuadPointMap::const_iterator, QuadPointMap::const_iterator> qp_range =
                    d_quad_point_map[ln].equal_range(i);
                if (qp_range.first != qp_range.second)
                {
                    const Vector X_cell(x_lower[0] + dx[0] * (static_cast<double>(i(0) - patch_lower(0)) + 0.5),
                                        x_lower[1] + dx[1] * (static_cast<double>(i(1) - patch_lower(1)) + 0.5)
#if (NDIM == 3)
                                            ,
                                        x_lower[2] + dx[2] * (static_cast<double>(i(2) - patch_lower(2)) + 0.5)
#endif
                                            );
                    if (U_cc_data)
                    {
                        for (QuadPointMap::const_iterator it = qp_range.first; it != qp_range.second; ++it)
                        {
                            const int& meter_num = it->second.meter_num;
                            const double& JxW = it->second.JxW;
                            const Vector& X = it->second.qp_xyz_current;
                            const Vector& normal = it->second.normal;
                            const Vector U = linear_interp<NDIM>(
                                X, i, X_cell, *U_cc_data, patch_lower, patch_upper, x_lower, x_upper, dx);
                            d_flow_values[meter_num] += (U.dot(normal)) * JxW;
                        }
                    }
                    if (U_sc_data)
                    {
                        for (QuadPointMap::const_iterator it = qp_range.first; it != qp_range.second; ++it)
                        {
                            const int& meter_num = it->second.meter_num;
                            const double& JxW = it->second.JxW;
                            const Vector& X = it->second.qp_xyz_current;
                            const Vector& normal = it->second.normal;
                            const Vector U =
                                linear_interp(X, i, X_cell, *U_sc_data, patch_lower, patch_upper, x_lower, x_upper, dx);
                            d_flow_values[meter_num] += (U.dot(normal)) * JxW;
                        }
                    }
                    if (P_cc_data)
                    {
                        for (QuadPointMap::const_iterator it = qp_range.first; it != qp_range.second; ++it)
                        {
                            const int& meter_num = it->second.meter_num;
                            const double& JxW = it->second.JxW;
                            const Vector& X = it->second.qp_xyz_current;
                            double P =
                                linear_interp(X, i, X_cell, *P_cc_data, patch_lower, patch_upper, x_lower, x_upper, dx);
                            d_mean_pressure_values[meter_num] += P * JxW;
                            A[meter_num] += JxW;
                            count_qp_2 += 1;
                        }
                    }
                }
            }
        }
    }

    /*perr << " num qp in map = " << count_qp << "\n";
    perr << " num qp after counting = " << count_qp_2 << "\n";
    perr << " actual num qp = " << d_num_quad_points[0] << "\n";
    perr << " flux before sum reduction = " << d_flow_values[0] << "\n";*/

    // Synchronize the values across all processes.
    SAMRAI_MPI::sumReduction(&d_flow_values[0], d_num_meters);
    SAMRAI_MPI::sumReduction(&d_mean_pressure_values[0], d_num_meters);
    SAMRAI_MPI::sumReduction(&A[0], d_num_meters);

    // Normalize the mean pressure.
    for (unsigned int jj = 0; jj < d_num_meters; ++jj)
    {
        d_mean_pressure_values[jj] /= A[jj];
    }

    // we need to compute the flow correction by calculating the contribution
    // from the velocity of each meter mesh.
    for (int jj = 0; jj < d_num_meters; ++jj)
    {
        // perr << " flux before correction = " << d_flow_values[jj] << "\n";

        // get displacement and velocity systems for meter mesh
        const LinearImplicitSystem& velocity_sys =
            d_meter_systems[jj]->get_system<LinearImplicitSystem>(IBFEMethod::VELOCITY_SYSTEM_NAME);
        const NumericVector<double>& velocity_coords = velocity_sys.get_vector("serial solution");
        const DofMap& dof_map = velocity_sys.get_dof_map();
        FEType fe_type = velocity_sys.variable_type(0);

        // set up FE objects
        UniquePtr<FEBase> fe_elem(FEBase::build(NDIM - 1, fe_type));
        QGauss qrule(NDIM - 1, fe_type.default_quadrature_order());
        fe_elem->attach_quadrature_rule(&qrule);

        //  for evaluating the velocity system
        const std::vector<Real>& JxW = fe_elem->get_JxW();
        const std::vector<std::vector<Real> >& phi = fe_elem->get_phi();
        const std::vector<libMesh::Point>& qp_points = fe_elem->get_xyz();
        std::vector<dof_id_type> dof_indices;
        DenseMatrix<double> vel_coords;

        // loop over elements again to compute mass flux and mean pressure
        double flux_correction = 0.0;
        double area = 0.0;
        MeshBase::const_element_iterator el = d_meter_meshes[jj]->active_local_elements_begin();
        const MeshBase::const_element_iterator end_el = d_meter_meshes[jj]->active_local_elements_end();
        for (; el != end_el; ++el)
        {
            const Elem* elem = *el;
            fe_elem->reinit(elem);

            // get dofs for displacement system and store in dense matrix
            vel_coords.resize(NDIM, phi.size());
            for (int d = 0; d < NDIM; ++d) // here d is the "variable number"
            {
                dof_map.dof_indices(elem, dof_indices, d);
                for (int nn = 0; nn < dof_indices.size(); ++nn)
                {
                    vel_coords(d, nn) = velocity_coords(dof_indices[nn]);
                }
            }

            // compute normal vector to element
            const libMesh::Point foo1 = *elem->node_ptr(1) - *elem->node_ptr(0);
            const libMesh::Point foo2 = *elem->node_ptr(2) - *elem->node_ptr(1);
            const libMesh::Point normal = (foo1.cross(foo2)).unit();

            area += 0.5 * (foo1.cross(foo2)).norm();

            // loop over quadrature points
            double vel_comp;
            for (int qp = 0; qp < qp_points.size(); ++qp)
            {
                for (int d = 0; d < NDIM; ++d)
                {
                    vel_comp = 0.0;
                    for (int nn = 0; nn < phi.size(); ++nn)
                    {
                        vel_comp += vel_coords(d, nn) * phi[nn][qp];
                    }
                    flux_correction += vel_comp * normal(d) * JxW[qp];
                }
            }
        }

        const double total_correction = SAMRAI_MPI::sumReduction(flux_correction);
        d_flow_values[jj] -= total_correction;

        // perr << " correction = " << total_correction << "\n";
        // perr << " flux after correction = " << d_flow_values[jj] << "\n";
        // perr << " area = " << A[jj] << "\n";

    } // loop over meters

    // write data
    outputData(data_time);

} // readInstrumentData

//**********************************************
// write out meshes to Exodus file
//**********************************************
void
IBFEInstrumentPanel::outputExodus(const int timestep, const double loop_time)
{
    for (int ii = 0; ii < d_num_meters; ++ii)
    {
        std::ostringstream mesh_output;
        mesh_output << d_plot_directory_name << "/"
                    << "" << d_meter_mesh_names[ii] << ".ex2";
        d_exodus_io[ii]->write_timestep(mesh_output.str(), *d_meter_systems[ii], timestep, loop_time);
    }
} // outputExodus

//**********************************************
// write out mesh nodes to data file
//**********************************************
void
IBFEInstrumentPanel::outputNodes()
{
    for (int ii = 0; ii < d_num_meters; ++ii)
    {
        std::ofstream stuff_stream;
        std::ostringstream node_output;
        node_output << d_plot_directory_name << "/"
                    << "" << d_meter_mesh_names[ii] << "_nodes.dat";
        if (SAMRAI_MPI::getRank() == 0)
        {
            stuff_stream.open(node_output.str().c_str());
            for (int dd = 0; dd < d_nodes[0].size(); ++dd)
            {
                stuff_stream << d_nodes[0][dd](0) << " " << d_nodes[0][dd](1) << " " << d_nodes[0][dd](2) << "\n";
            }
            stuff_stream.close();
        }
    }
} // outputNodes

//**********************************************
// get data from input file
//**********************************************
void
IBFEInstrumentPanel::getFromInput(Pointer<Database> db)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(db);
#endif
    if (db->keyExists("plot_directory_name")) d_plot_directory_name = db->getString("plot_directory_name");
    if (db->keyExists("nodeset_IDs_for_meters"))
        d_nodeset_IDs_for_meters = db->getIntegerArray("nodeset_IDs_for_meters");
    if (db->keyExists("meter_mesh_quad_order"))
        d_quad_order = Utility::string_to_enum<Order>(db->getStringWithDefault("meter_mesh_quad_order", "SECOND"));
    return;
} // getFromInput

//**********************************************
// update system data
//**********************************************
void
IBFEInstrumentPanel::updateSystemData(IBAMR::IBFEMethod* ib_method_ops, const int meter_mesh_number)
{
    // get the coordinate mapping system and velocity systems for the parent mesh
    const FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager(d_part);
    const EquationSystems* equation_systems = fe_data_manager->getEquationSystems();
    const System& dX_system = equation_systems->get_system(IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
    // TO DO: find a better way to do this
    std::vector<double> dX_coords_parent;
    dX_system.update_global_solution(dX_coords_parent);
    const System& U_system = equation_systems->get_system(IBFEMethod::VELOCITY_SYSTEM_NAME);
    std::vector<double> U_coords_parent;
    U_system.update_global_solution(U_coords_parent);

    // get displacement and velocity systems for meter mesh
    LinearImplicitSystem& velocity_sys =
        d_meter_systems[meter_mesh_number]->get_system<LinearImplicitSystem>(IBFEMethod::VELOCITY_SYSTEM_NAME);
    const unsigned int velocity_sys_num = velocity_sys.number();
    NumericVector<double>& velocity_solution = *velocity_sys.solution;
    NumericVector<double>& velocity_coords = velocity_sys.get_vector("serial solution");

    LinearImplicitSystem& displacement_sys =
        d_meter_systems[meter_mesh_number]->get_system<LinearImplicitSystem>(IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
    const unsigned int displacement_sys_num = displacement_sys.number();
    NumericVector<double>& displacement_solution = *displacement_sys.solution;
    NumericVector<double>& displacement_coords = displacement_sys.get_vector("serial solution");

    // loop over all nodes in meter mesh
    for (int ii = 0; ii < d_num_nodes[meter_mesh_number]; ++ii)
    {
        // get node on meter mesh
        const Node* node = &d_meter_meshes[meter_mesh_number]->node_ref(ii);

        // get corresponding dofs on parent mesh
        std::vector<double> U_dofs;
        U_dofs.resize(NDIM);
        std::vector<double> dX_dofs;
        dX_dofs.resize(NDIM);

        for (int d = 0; d < NDIM; ++d)
        {
            U_dofs[d] = U_coords_parent[d_U_dof_idx[meter_mesh_number][ii][d]];
            dX_dofs[d] = dX_coords_parent[d_dX_dof_idx[meter_mesh_number][ii][d]];
        }

        // set dofs in meter mesh to correspond to the same values
        // as in the parent mesh
        for (int d = 0; d < NDIM; ++d)
        {
            const int vel_dof_idx = node->dof_number(velocity_sys_num, d, 0);
            velocity_coords.set(vel_dof_idx, U_dofs[d]);
            const int disp_dof_idx = node->dof_number(displacement_sys_num, d, 0);
            displacement_coords.set(disp_dof_idx, dX_dofs[d]);
        }
    }

    // populate solution vector in system also... why not?
    MeshBase::const_node_iterator node_it = d_meter_meshes[meter_mesh_number]->local_nodes_begin();
    const MeshBase::const_node_iterator end_node_it = d_meter_meshes[meter_mesh_number]->local_nodes_end();
    for (; node_it != end_node_it; ++node_it)
    {
        const Node* node = *node_it;
        for (int d = 0; d < NDIM; ++d)
        {
            const int vel_dof_idx = node->dof_number(velocity_sys_num, d, 0);
            velocity_solution.set(vel_dof_idx, velocity_coords(vel_dof_idx));
            const int disp_dof_idx = node->dof_number(displacement_sys_num, d, 0);
            displacement_solution.set(disp_dof_idx, displacement_coords(disp_dof_idx));
        }
    }
    velocity_solution.close();
    displacement_solution.close();

} // updateSystemData

//**********************************************
// output mesh info
//**********************************************
void
IBFEInstrumentPanel::outputMeterMeshes(const int timestep_num, const double data_time)
{
    // things to do at initial timestep
    if (timestep_num == 1) outputNodes();
    outputExodus(timestep_num, data_time);
} // outputMeterMeshes

//**********************************************
// output data read from the meters
//**********************************************
void
IBFEInstrumentPanel::outputData(const double data_time)
{
    if (SAMRAI_MPI::getRank() == 0)
    {
        d_mean_pressure_stream << data_time;
        d_flux_stream << data_time;
        for (int jj = 0; jj < d_num_meters; ++jj)
        {
            d_mean_pressure_stream << " " << d_mean_pressure_values[jj];
            d_flux_stream << " " << d_flow_values[jj];
        }
        d_mean_pressure_stream << "\n";
        d_flux_stream << "\n";
    }
} // outputData

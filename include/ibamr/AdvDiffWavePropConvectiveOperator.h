#ifndef included_AdvDiffWavePropConvectiveOperator
#define included_AdvDiffWavePropConvectiveOperator

#include <string>
#include <vector>

#include "CellVariable.h"
#include "CoarsenAlgorithm.h"
#include "FaceVariable.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "RefineAlgorithm.h"
#include "RefinePatchStrategy.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/ibamr_enums.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include <ostream>
#include <stddef.h>
#include <string>
#include <vector>

#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellDataFactory.h"
#include "CellVariable.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "CoarsenSchedule.h"
#include "FaceData.h"
#include "FaceVariable.h"
#include "IBAMR_config.h"
#include "Index.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefinePatchStrategy.h"
#include "RefineSchedule.h"
#include "SAMRAIVectorReal.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "ibamr/AdvDiffCenteredConvectiveOperator.h"
#include "ibamr/AdvDiffPhysicalBoundaryUtilities.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CartExtrapPhysBdryOp.h"
#include "ibtk/CartSideRobinPhysBdryOp.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/PhysicalBoundaryUtilities.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <boost/concept_check.hpp>
#include <unsupported/Eigen/MatrixFunctions>

namespace SAMRAI
{
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
namespace xfer
{
template <int DIM>
class CoarsenSchedule;
template <int DIM>
class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI

namespace IBAMR
{
class AdvDiffWavePropConvectiveOperator : public ConvectiveOperator
{
public:
    // Constructor
    AdvDiffWavePropConvectiveOperator(const std::string& object_name,
                                      SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                      const ConvectiveDifferencingType difference_form,
                                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& conc_bc_coefs);
    // Destructor
    ~AdvDiffWavePropConvectiveOperator();

    static SAMRAI::tbox::Pointer<ConvectiveOperator>
    allocate_operator(const std::string& object_name,
                      SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                      const ConvectiveDifferencingType difference_form,
                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& conc_bc_coefs)
    {
        return new AdvDiffWavePropConvectiveOperator(object_name, Q_var, input_db, difference_form, conc_bc_coefs);
    }

    void applyConvectiveOperator(int Q_idx, int Y_idx);

    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out);

    void deallocateOperatorState();

private:
    AdvDiffWavePropConvectiveOperator();

    AdvDiffWavePropConvectiveOperator(const AdvDiffWavePropConvectiveOperator& from);

    AdvDiffWavePropConvectiveOperator& operator=(const AdvDiffWavePropConvectiveOperator& that);

    // Data communication algorithms, operators, and schedules.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > d_coarsen_alg_Q;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > d_coarsen_scheds_Q;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_ghostfill_alg_Q;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefinePatchStrategy<NDIM> > d_ghostfill_strategy_Q;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_ghostfill_scheds_Q;
    std::string d_outflow_bdry_extrap_type;

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln, d_finest_ln;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Q_var;
    unsigned int d_Q_data_depth;
    int d_Q_scratch_idx;

    const std::vector<RobinBcCoefStrategy<NDIM>*> d_conc_bc_coefs;
    // Reconstruction Values
    void calculateWeights();
    int d_k;
    std::vector<double> d_interp_weights_f, d_smooth_weights;
    std::vector<std::vector<double> > d_interp_weights;

    ConvectiveDifferencingType d_difference_form;
};
}

#endif

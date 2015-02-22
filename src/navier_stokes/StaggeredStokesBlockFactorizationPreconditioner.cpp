// Filename: StaggeredStokesBlockFactorizationPreconditioner.cpp
// Created on 22 Sep 2008 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <ostream>
#include <string>

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/MultiblockDataTranslator.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/solv/PoissonSpecifications.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/xfer/VariableFillPattern.h"
#include "ibamr/StaggeredStokesBlockFactorizationPreconditioner.h"
#include "ibamr/StaggeredStokesBlockPreconditioner.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CellNoCornersFillPattern.h"
#include "ibtk/GeneralSolver.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/PoissonSolver.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
static const int SIDEG = 1;

// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "NONE";
static const bool USE_CF_INTERPOLATION = true;
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Timers.
static Timer* t_solve_system;
static Timer* t_initialize_solver_state;
static Timer* t_deallocate_solver_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesBlockFactorizationPreconditioner::StaggeredStokesBlockFactorizationPreconditioner(
    const std::string& object_name,
    Pointer<Database> /*input_db*/,
    const std::string& /*default_options_prefix*/)
    : StaggeredStokesBlockPreconditioner(/*needs_velocity_solver*/ true,
                                         /*needs_pressure_solver*/ true),
      d_P_bdry_fill_op(NULL), d_no_fill_op(NULL), d_U_var(NULL), d_F_U_mod_idx(-1), d_P_var(NULL), d_P_scratch_idx(-1)
{
    GeneralSolver::init(object_name, /*homogeneous_bc*/ true);

    // Present implementation requires zero initial guess and can perform only
    // one iteration.
    d_initial_guess_nonzero = false;
    d_max_iterations = 1;

    // Setup variables.
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    Pointer<VariableContext> context = var_db->getContext(d_object_name + "::CONTEXT");

    const std::string U_var_name = d_object_name + "::U";
    d_U_var = var_db->getVariable(U_var_name);
    if (d_U_var)
    {
        d_F_U_mod_idx = var_db->mapVariableAndContextToIndex(d_U_var, context);
    }
    else
    {
        d_U_var = new SideVariable<double>(DIM, U_var_name);
        d_F_U_mod_idx = var_db->registerVariableAndContext(d_U_var, context, IntVector(DIM, SIDEG));
    }

    const std::string P_var_name = d_object_name + "::P";
    d_P_var = var_db->getVariable(P_var_name);
    if (d_P_var)
    {
        d_P_scratch_idx = var_db->mapVariableAndContextToIndex(d_P_var, context);
    }
    else
    {
        d_P_var = new CellVariable<double>(DIM, P_var_name);
        d_P_scratch_idx = var_db->registerVariableAndContext(d_P_var, context, IntVector(DIM, CELLG));
    }

    // Setup Timers.
    IBAMR_DO_ONCE(t_solve_system = TimerManager::getManager()->getTimer(
                      "IBAMR::StaggeredStokesBlockFactorizationPreconditioner::solveSystem()");
                  t_initialize_solver_state = TimerManager::getManager()->getTimer(
                      "IBAMR::StaggeredStokesBlockFactorizationPreconditioner::initializeSolverState()");
                  t_deallocate_solver_state = TimerManager::getManager()->getTimer(
                      "IBAMR::StaggeredStokesBlockFactorizationPreconditioner::deallocateSolverState("
                      ")"););
    return;
} // StaggeredStokesBlockFactorizationPreconditioner

StaggeredStokesBlockFactorizationPreconditioner::~StaggeredStokesBlockFactorizationPreconditioner()
{
    deallocateSolverState();
    return;
} // ~StaggeredStokesBlockFactorizationPreconditioner

bool StaggeredStokesBlockFactorizationPreconditioner::solveSystem(SAMRAIVectorReal<double>& x,
                                                                  SAMRAIVectorReal<double>& b)
{
    IBAMR_TIMER_START(t_solve_system);

    // Initialize the solver (if necessary).
    const bool deallocate_at_completion = !d_is_initialized;
    if (!d_is_initialized) initializeSolverState(x, b);

    // Determine whether we are solving a steady-state problem.
    const bool steady_state =
        d_U_problem_coefs.cIsZero() ||
        (d_U_problem_coefs.cIsConstant() && MathUtilities<double>::equalEps(d_U_problem_coefs.getCConstant(), 0.0));

    // Get the vector components.
    const int F_U_idx = b.getComponentDescriptorIndex(0);
    const int F_P_idx = b.getComponentDescriptorIndex(1);

    const Pointer<Variable>& F_U_var = b.getComponentVariable(0);
    const Pointer<Variable>& F_P_var = b.getComponentVariable(1);

    Pointer<SideVariable<double> > F_U_sc_var = F_U_var;
    Pointer<CellVariable<double> > F_P_cc_var = F_P_var;

    const int U_idx = x.getComponentDescriptorIndex(0);
    const int P_idx = x.getComponentDescriptorIndex(1);

    const Pointer<Variable>& U_var = x.getComponentVariable(0);
    const Pointer<Variable>& P_var = x.getComponentVariable(1);

    Pointer<SideVariable<double> > U_sc_var = U_var;
    Pointer<CellVariable<double> > P_cc_var = P_var;

    // Setup the component solver vectors.
    Pointer<SAMRAIVectorReal<double> > F_U_mod_vec;
    F_U_mod_vec = new SAMRAIVectorReal<double>(d_object_name + "::F_U_mod", d_hierarchy, d_coarsest_ln, d_finest_ln);
    F_U_mod_vec->addComponent(d_U_var, d_F_U_mod_idx, d_velocity_wgt_idx, d_velocity_data_ops);

    Pointer<SAMRAIVectorReal<double> > U_vec;
    U_vec = new SAMRAIVectorReal<double>(d_object_name + "::U", d_hierarchy, d_coarsest_ln, d_finest_ln);
    U_vec->addComponent(U_sc_var, U_idx, d_velocity_wgt_idx, d_velocity_data_ops);

    Pointer<SAMRAIVectorReal<double> > P_scratch_vec;
    P_scratch_vec =
        new SAMRAIVectorReal<double>(d_object_name + "::P_scratch", d_hierarchy, d_coarsest_ln, d_finest_ln);
    P_scratch_vec->addComponent(d_P_var, d_P_scratch_idx, d_pressure_wgt_idx, d_pressure_data_ops);

    Pointer<SAMRAIVectorReal<double> > F_P_vec;
    F_P_vec = new SAMRAIVectorReal<double>(d_object_name + "::F_P", d_hierarchy, d_coarsest_ln, d_finest_ln);
    F_P_vec->addComponent(F_P_cc_var, F_P_idx, d_pressure_wgt_idx, d_pressure_data_ops);

    Pointer<SAMRAIVectorReal<double> > P_vec;
    P_vec = new SAMRAIVectorReal<double>(d_object_name + "::P", d_hierarchy, d_coarsest_ln, d_finest_ln);
    P_vec->addComponent(P_cc_var, P_idx, d_pressure_wgt_idx, d_pressure_data_ops);

    // Setup the interpolation transaction information.
    Pointer<VariableFillPattern> fill_pattern(new CellNoCornersFillPattern(CELLG, false, false, true));
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent P_transaction_comp(P_idx,
                                                         DATA_REFINE_TYPE,
                                                         USE_CF_INTERPOLATION,
                                                         DATA_COARSEN_TYPE,
                                                         BDRY_EXTRAP_TYPE,
                                                         CONSISTENT_TYPE_2_BDRY,
                                                         d_P_bc_coef,
                                                         fill_pattern);
    InterpolationTransactionComponent P_scratch_transaction_comp(d_P_scratch_idx,
                                                                 DATA_REFINE_TYPE,
                                                                 USE_CF_INTERPOLATION,
                                                                 DATA_COARSEN_TYPE,
                                                                 BDRY_EXTRAP_TYPE,
                                                                 CONSISTENT_TYPE_2_BDRY,
                                                                 d_P_bc_coef,
                                                                 fill_pattern);

    // (1) Solve the pressure sub-problem by applying inv(S^) to F_P, in which
    // S^ is the approximate Schur complement.
    //
    // The Schur complement S is
    //
    //    S = D inv(rho/dt - K*mu*L) G
    //
    // We obtain S^ by assuming that
    //
    //    D inv(rho/dt - K*mu*L) G ~ L_p inv(rho/dt - K*mu*L_p)
    //
    // in which L_p = D*G.
    //
    // We treat two cases:
    //
    // (i) rho/dt = 0.
    //
    // In this case,
    //
    //    inv(S^) = inv(L_p inv(-K*mu*L_p)) = -K*mu
    //
    // so that
    //
    //    P := -K*mu*F_P
    //
    // (ii) rho/dt != 0.
    //
    // In this case, we make the further approximation that
    //
    //    L_p ~ rho L_rho = rho (D (1/rho) G)
    //
    // so that
    //
    //    inv(S^) = (1/dt) inv(L_rho) - K*mu
    //
    // and
    //
    //    P := [(1/dt) inv(L_rho) - K*mu] F_P
    //
    // NOTE: d_U_problem_coefs.getCConstant() == rho/dt
    //       d_U_problem_coefs.getDConstant() == -K*mu
    //
    // in which K depends on the form of the time stepping scheme.
    if (steady_state)
    {
        d_pressure_data_ops->scale(P_idx, d_U_problem_coefs.getDConstant(), F_P_idx);
    }
    else
    {
        d_pressure_solver->setHomogeneousBc(true);
        LinearSolver* p_pressure_solver = dynamic_cast<LinearSolver*>(d_pressure_solver.getPointer());
        if (p_pressure_solver) p_pressure_solver->setInitialGuessNonzero(false);
        d_pressure_solver->solveSystem(*P_scratch_vec, *F_P_vec); // P_scratch_idx := -inv(L_rho)*F_P
        d_pressure_data_ops->linearSum(
            P_idx, -1.0 / getDt(), d_P_scratch_idx, d_U_problem_coefs.getDConstant(), F_P_idx);
    }
    d_P_bdry_fill_op->resetTransactionComponent(P_transaction_comp);
    d_P_bdry_fill_op->fillData(d_pressure_solver->getSolutionTime());
    d_P_bdry_fill_op->resetTransactionComponent(P_scratch_transaction_comp);

    // (2) Solve the velocity sub-problem.
    //
    // U := inv(rho/dt - K*mu*L) * [F_U - G P]
    static const bool cf_bdry_synch = true;
    d_hier_math_ops->grad(d_F_U_mod_idx,
                          d_U_var,
                          cf_bdry_synch,
                          -1.0,
                          P_idx,
                          P_cc_var,
                          d_no_fill_op,
                          d_pressure_solver->getSolutionTime(),
                          1.0,
                          F_U_idx,
                          F_U_sc_var);
    d_velocity_solver->setHomogeneousBc(true);
    LinearSolver* p_velocity_solver = dynamic_cast<LinearSolver*>(d_velocity_solver.getPointer());
    if (p_velocity_solver) p_velocity_solver->setInitialGuessNonzero(false);
    d_velocity_solver->solveSystem(*U_vec, *F_U_mod_vec);

    // Account for nullspace vectors.
    correctNullspace(U_vec, P_vec);

    // Deallocate the solver (if necessary).
    if (deallocate_at_completion) deallocateSolverState();

    IBAMR_TIMER_STOP(t_solve_system);
    return true;
} // solveSystem

void StaggeredStokesBlockFactorizationPreconditioner::initializeSolverState(const SAMRAIVectorReal<double>& x,
                                                                            const SAMRAIVectorReal<double>& b)
{
    IBAMR_TIMER_START(t_initialize_solver_state);

    if (d_is_initialized) deallocateSolverState();

    // Parent class initialization.
    StaggeredStokesBlockPreconditioner::initializeSolverState(x, b);

    // Setup hierarchy operators.
    Pointer<VariableFillPattern> fill_pattern(new CellNoCornersFillPattern(CELLG, false, false, true));
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent P_scratch_component(d_P_scratch_idx,
                                                          DATA_REFINE_TYPE,
                                                          USE_CF_INTERPOLATION,
                                                          DATA_COARSEN_TYPE,
                                                          BDRY_EXTRAP_TYPE,
                                                          CONSISTENT_TYPE_2_BDRY,
                                                          d_P_bc_coef,
                                                          fill_pattern);
    d_P_bdry_fill_op = new HierarchyGhostCellInterpolation();
    d_P_bdry_fill_op->setHomogeneousBc(true);
    d_P_bdry_fill_op->initializeOperatorState(P_scratch_component, d_hierarchy);

    // Allocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel> level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_F_U_mod_idx)) level->allocatePatchData(d_F_U_mod_idx);
        if (!level->checkAllocated(d_P_scratch_idx)) level->allocatePatchData(d_P_scratch_idx);
    }

    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_solver_state);
    return;
} // initializeSolverState

void StaggeredStokesBlockFactorizationPreconditioner::deallocateSolverState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_solver_state);

    // Parent class deallocation.
    StaggeredStokesBlockPreconditioner::deallocateSolverState();

    // Deallocate hierarchy operators.
    d_P_bdry_fill_op.setNull();

    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel> level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_F_U_mod_idx)) level->deallocatePatchData(d_F_U_mod_idx);
        if (level->checkAllocated(d_P_scratch_idx)) level->deallocatePatchData(d_P_scratch_idx);
    }

    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_solver_state);
    return;
} // deallocateSolverState

void StaggeredStokesBlockFactorizationPreconditioner::setInitialGuessNonzero(bool initial_guess_nonzero)
{
    if (initial_guess_nonzero)
    {
        TBOX_ERROR(d_object_name + "::setInitialGuessNonzero()\n"
                   << "  class IBAMR::StaggeredStokesBlockFactorizationPreconditioner requires a "
                      "zero initial guess" << std::endl);
    }
    return;
} // setInitialGuessNonzero

void StaggeredStokesBlockFactorizationPreconditioner::setMaxIterations(int max_iterations)
{
    if (max_iterations != 1)
    {
        TBOX_ERROR(d_object_name + "::setMaxIterations()\n"
                   << "  class IBAMR::StaggeredStokesBlockFactorizationPreconditioner only "
                      "performs a single iteration" << std::endl);
    }
    return;
} // setMaxIterations

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

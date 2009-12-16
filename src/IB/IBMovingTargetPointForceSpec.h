#ifndef included_IBMovingTargetPointForceSpec
#define included_IBMovingTargetPointForceSpec

// Filename: IBMovingTargetPointForceSpec.h
// Last modified: <15.Dec.2009 19:21:06 griffith@boyce-griffiths-mac-pro.local>
// Created on 14 Aug 2008 by Boyce Griffith (boyce@dm-linux.maths.gla.ac.uk)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/DynamicArena.h>
#include <ibtk/Stashable.h>

// SAMRAI INCLUDES
#include <tbox/AbstractStream.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBMovingTargetPointForceSpec encapsulates the data necessary to
 * compute the penalty force generated by a single moving target point.
 */
class IBMovingTargetPointForceSpec
    : public IBTK::Stashable
{
public:
    /*!
     * \brief Register this class and its factory class with the singleton
     * IBTK::StashableManager object.  This method must be called before any
     * IBMovingTargetPointForceSpec objects are created.
     *
     * \note This method is collective on all MPI processes.  This is done to
     * ensure that all processes employ the same stashable ID for the
     * IBMovingTargetPointForceSpec class.
     */
    static void
    registerWithStashableManager();

    /*!
     * \brief Returns a boolean indicating whether the class has been registered
     * with the singleton IBTK::StashableManager object.
     */
    static bool
    getIsRegisteredWithStashableManager();

    /*!
     * \brief Operator new.
     */
    static void*
    operator new(
        std::size_t size);

    /*!
     * \brief Operator delete.
     */
    static void
    operator delete(
        void* ptr,
        std::size_t size);

    /*!
     * \brief Default constructor.
     */
    IBMovingTargetPointForceSpec(
        const int master_idx=-1,
        const double& kappa_target=0.0,
        const double& eta_target=0.0,
        const int spec_fcn_idx=-1,
        const std::vector<double>& periodic_shift=std::vector<double>(NDIM,0.0));

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBMovingTargetPointForceSpec();

    /*!
     * \return A const refernce to the master node index.
     */
    const int&
    getMasterNodeIndex() const;

    /*!
     * \return A non-const reference to the master node index.
     */
    int&
    getMasterNodeIndex();

    /*!
     * \return A const reference to the stiffness of the spring attached to the
     * target point.
     */
    const double&
    getStiffness() const;

    /*!
     * \return A non-const reference to the stiffness of the spring attached to
     * the target point.
     */
    double&
    getStiffness();

    /*!
     * \return A const reference to the damping factor of the spring attached to
     * the target point.
     */
    const double&
    getDamping() const;

    /*!
     * \return A non-const reference to the damping factor of the spring
     * attached to the target point.
     */
    double&
    getDamping();

    /*!
     * \return A const reference to the index of the target point position and
     * velocity specificaiton function associated with the target point.
     */
    const int&
    getPositionAndVelocityFunctionIndex() const;

    /*!
     * \return A non-const reference to the index of the target point position
     * and velocity specification function associated with the target point.
     */
    int&
    getPositionAndVelocityFunctionIndex();

    /*!
     * \return A const reference to the periodic shift associated with the
     * target point.
     */
    const std::vector<double>&
    getPeriodicShift() const;

    /*!
     * \return A non-const reference to the periodic shift associated with the
     * target point.
     */
    std::vector<double>&
    getPeriodicShift();

    /*!
     * \brief Return the unique identifier used to specify the IBTK::StashableFactory
     * object used by the IBTK::StashableManager to extract Stashable objects from
     * data streams.
     */
    virtual int
    getStashableID() const;

    /*!
     * \brief Return an upper bound on the amount of space required to pack the
     * object to a buffer.
     */
    virtual size_t
    getDataStreamSize() const;

    /*!
     * \brief Pack data into the output stream.
     */
    virtual void
    packStream(
        SAMRAI::tbox::AbstractStream& stream);

    /*!
     * \brief Indicate that the object is being shifted across a periodic
     * boundary.
     */
    virtual void
    registerPeriodicShift(
        const SAMRAI::hier::IntVector<NDIM>& offset,
        const std::vector<double>& displacement);

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBMovingTargetPointForceSpec(
        const IBMovingTargetPointForceSpec& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBMovingTargetPointForceSpec&
    operator=(
        const IBMovingTargetPointForceSpec& that);

    /*!
     * Indicates whether the factory has been registered with the
     * IBTK::StashableManager.
     */
    static bool s_registered_factory;

    /*!
     * The stashable ID for this object type.
     */
    static int s_stashable_id;

    /*!
     * Memory arena for allocating objects.
     */
    static IBTK::DynamicArena s_arena;

    /*!
     * Data required to define the target point penalty forces.
     */
    int d_master_idx;
    double d_kappa_target, d_eta_target;
    int d_spec_fcn_idx;
    std::vector<double> d_periodic_shift;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibamr/IBMovingTargetPointForceSpec.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBMovingTargetPointForceSpec

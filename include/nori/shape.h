
#pragma once

#include <nori/object.h>
#include <nori/frame.h>
#include <nori/bbox.h>

NORI_NAMESPACE_BEGIN

class Shape;
class Sampler;

/**
* \brief Intersection data structure
*
* This data structure records local information about a ray-triangle intersection.
* This includes the position, traveled ray distance, uv coordinates, as well
* as well as two local coordinate frames (one that corresponds to the true
* geometry, and one that is used for shading computations).
*/
struct Intersection {
	/// Position of the surface intersection
	Point3f p;
	/// Unoccluded distance along the ray
	float t;
	/// UV coordinates, if any
	Point2f uv;
	/// Shading frame (based on the shading normal)
	Frame shFrame;
	/// Geometric frame (based on the true geometry)
	Frame geoFrame;
	/// Pointer to the associated mesh
	const Shape *shape;

	/// Create an uninitialized intersection record
	Intersection() : shape(nullptr) { }

	/// Transform a direction vector into the local shading frame
	Vector3f toLocal(const Vector3f &d) const {
		return shFrame.toLocal(d);
	}

	/// Transform a direction vector from local to world coordinates
	Vector3f toWorld(const Vector3f &d) const {
		return shFrame.toWorld(d);
	}

	/// Return a human-readable summary of the intersection record
	std::string toString() const;
};

class Shape : public NoriObject
{
public:
	/// Release all memory
	virtual ~Shape();

	/// Initialize internal data structures (called once by the XML parser)
	virtual void activate() override;

	/// Register a child object (e.g. a BSDF) with the shape
	virtual void addChild(NoriObject *obj) override;

	//// Return the centroid of the given triangle
	virtual Point3f getCentroid(uint32_t index) const = 0;

	/// Return the total number of primitives in this shape
	virtual uint32_t getPrimitiveCount() const = 0;

	virtual float getArea() const = 0;

	virtual Point3f sample(Sampler* sampler, Normal3f& normal) const = 0;
	virtual Vector3f sampleSolidAngle(Sampler* sampler, Point3f& x, Normal3f& normal, float& pWi) const = 0;

	/** \brief Ray-triangle intersection test
	*
	* Uses the algorithm by Moeller and Trumbore discussed at
	* <tt>http://www.acm.org/jgt/papers/MollerTrumbore97/code.html</tt>.
	*
	* Note that the test only applies to a single triangle in the mesh.
	* An acceleration data structure like \ref BVH is needed to search
	* for intersections against many triangles.
	*
	* \param index
	*    Index of the triangle that should be intersected
	* \param ray
	*    The ray segment to be used for the intersection query
	* \param t
	*    Upon success, \a t contains the distance from the ray origin to the
	*    intersection point,
	* \param u
	*   Upon success, \c u will contain the 'U' component of the intersection
	*   in barycentric coordinates
	* \param v
	*   Upon success, \c v will contain the 'V' component of the intersection
	*   in barycentric coordinates
	* \return
	*   \c true if an intersection has been detected
	*/
	virtual bool rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const = 0;

	virtual void computeIntersectionInfo(uint32_t index, const Ray3f &ray, Intersection & its) const = 0;

	//// Return an axis-aligned bounding box of the entire shape
	const BoundingBox3f &getBoundingBox() const { return m_bbox; }

	//// Return an axis-aligned bounding box containing the given triangle
	virtual BoundingBox3f getBoundingBox(uint32_t index) const = 0;

	/// Is this shape an area emitter?
	bool isEmitter() const { return m_emitter != nullptr; }

	/// Return a pointer to an attached area emitter instance
	Emitter *getEmitter() { return m_emitter; }

	/// Return a pointer to an attached area emitter instance (const version)
	const Emitter *getEmitter() const { return m_emitter; }

	/// Return a pointer to the BSDF associated with this shape
	const BSDF *getBSDF() const { return m_bsdf; }

	/// Return a human-readable summary of this instance
	virtual std::string toString() const override;

	/**
	* \brief Return the type of object (i.e. Mesh/BSDF/etc.)
	* provided by this instance
	* */
	virtual EClassType getClassType() const override { return EShape; }

protected:
	/// Create an empty shape
	Shape();

protected:
	BSDF          *m_bsdf = nullptr;      ///< BSDF of the surface
	Emitter       *m_emitter = nullptr;   ///< Associated emitter, if any
	BoundingBox3f  m_bbox;                ///< Bounding box of the shape
};

NORI_NAMESPACE_END
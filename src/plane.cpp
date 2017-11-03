
#include <nori/shape.h>
#include <nori/sampler.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

/**
* \brief Plane
*
* Defined by:
*     - center
*     - width and height
      - normalRotation (degX, degY, degZ)
*
* Local space :
*     n = (0,0,1) @ deg 0,0,0
*     c = (0,0,0)
*     w = 2 (w/2 = 1)
*     h = 2
*     xy plane [-1,1]
*
*/
class Plane : public Shape {
public:
    Plane(const PropertyList &propList):
        m_center(propList.getPoint("center", Point3f(0.0f))),
        m_width(propList.getFloat("width", 2.0f)),
        m_height(propList.getFloat("height", 2.0f)),
        m_normalRotation(propList.getVector("normalRotation", Vector3f(0.0f)))
    {
        // create transform from c,w,h,normal
        Eigen::Affine3f objectToWorld;
        objectToWorld.setIdentity();

        // - translate
        Eigen::Affine3f trMat(Eigen::Translation3f(m_center.x(), m_center.y(), m_center.z()));
        objectToWorld *= trMat.matrix();

        // - rotate
        float angleX = degToRad(m_normalRotation.x());
        float angleY = degToRad(m_normalRotation.y());
        float angleZ = degToRad(m_normalRotation.z());
        Eigen::Affine3f rotMat;
        rotMat = Eigen::AngleAxisf(angleX, Vector3f::UnitX())
               * Eigen::AngleAxisf(angleY, Vector3f::UnitY())
               * Eigen::AngleAxisf(angleZ, Vector3f::UnitZ());
        objectToWorld *= rotMat.matrix();

        // - scale
        Eigen::Affine3f scaleMat(Eigen::DiagonalMatrix<float, 3>(m_width/2.0f, m_height/2.0f, 1.0f));
        objectToWorld *= scaleMat.matrix();

        m_objectToWorld = objectToWorld.matrix();

        // update c,w,h
        m_center = m_objectToWorld * Point3f(0.0f,0.0f,0.0f);
        Vector3f dx = m_objectToWorld * Vector3f(2.0f,0.0f,0.0f);
        Vector3f dy = m_objectToWorld * Vector3f(0.0f,2.0f,0.0f);
        m_width = dx.norm();
        m_height = dy.norm();

        // setup normalFrame
        Normal3f normal = m_objectToWorld * Normal3f(0.0f,0.0f,1.0f);
        normal.normalize();
        dx.normalize();
        dy.normalize();
        m_normalFrame = Frame(dx,dy,normal);

        // setup bbox
        m_bbox.expandBy( m_objectToWorld * Point3f(-1.0f,  1.0f, 0.0f) );
        m_bbox.expandBy( m_objectToWorld * Point3f( 1.0f,  1.0f, 0.0f) );
        m_bbox.expandBy( m_objectToWorld * Point3f( 1.0f, -1.0f, 0.0f) );
        m_bbox.expandBy( m_objectToWorld * Point3f(-1.0f, -1.0f, 0.0f) );
    }

    virtual uint32_t getPrimitiveCount() const override {
        return 1;
    }

    virtual BoundingBox3f getBoundingBox(uint32_t index) const override {
        return m_bbox;
    }

    virtual Point3f getCentroid(uint32_t index) const override {
        return m_center;
    }

    virtual float getArea() const override {
        return m_width * m_height;
    }

    virtual Point3f sample(Sampler* sampler, Normal3f& normal) const override {
        Point2f sample = sampler->next2D();
        Point3f p(sample.x(), sample.y(), 0.0f);

        // sample space to local space
        Point3f y = Point3f(p.x()*2.0f, p.y()*2.0f, 0.0f);
        y += Point3f(-1.0f, -1.0f, 0); // samples origin are not at 0,0,0
        // local space to world space
        y = m_objectToWorld * y;

        normal = m_normalFrame.n;
        return y;
    }

    virtual Vector3f sampleSolidAngle(Sampler* sampler, Point3f& x, Normal3f& normal, float& pWi, Point3f& y) const override {
        throw NoriException("Unimplemented Plane::sampleSolidAngle() !!!");
    }

    virtual bool rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const override {
        //see: https://www.cl.cam.ac.uk/teaching/1999/AGraphHCI/SMAG/node2.html#eqn:vectplane

        // ray-plane intersection in local space
        Ray3f ray_(ray);
        ray_.o = m_objectToWorld.inverse() * ray_.o;
        ray_.d = m_objectToWorld.inverse() * ray_.d;

        t = -ray_.o.z() / ray_.d.z();

        if (t > ray_.maxt || t < ray_.mint) {
            return false;
        }

        Point3f P(ray_(t)); // local

        return P.x() >= -1.0f && P.x() <= 1.0f &&
               P.y() >= -1.0f && P.y() <= 1.0f;
    }

    virtual void computeIntersectionInfo(uint32_t index, const Ray3f &ray, Intersection & its) const override {
        its.p = ray(its.t); // world

        its.geoFrame = m_normalFrame;
        its.shFrame = m_normalFrame;
    }

    virtual std::string toString() const override {
        return tfm::format(
        "%s\n"
        "Plane[\n"
        "center = %s\n"
        "width = %d\n"
        "height = %d\n"
        "normalRot = %s\n"
        "]",
        Shape::toString(),
        m_center.toString(),
        m_width,
        m_height,
        m_normalRotation.toString()
        );
    }

    virtual std::string getType() const override { return "plane"; }

protected:
    // params
    Point3f m_center;
    float m_width;
    float m_height;
    const Vector3f m_normalRotation;

    //
    Transform m_objectToWorld;
    Frame m_normalFrame;
};

NORI_REGISTER_CLASS(Plane, "plane");
NORI_NAMESPACE_END

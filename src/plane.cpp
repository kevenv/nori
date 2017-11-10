
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

    // Implementation of "An Area-Preserving Parametrization for Spherical Rectangles"
    struct SphQuad {
        Vector3f o, x, y, z; // local reference system 'R'
        float z0, z0sq; // rectangle coords in 'R'
        float x0, y0, y0sq;
        float x1, y1, y1sq;
        float b0, b1, b0sq, k; // misc precomputed constants
        float S; // solid angle of 'Q'
    };

    virtual Vector3f sampleSolidAngle(Sampler* sampler, Point3f& x, Normal3f& normal, float& pWi, Point3f& y) const override {
        // make sure rect light only emmit from one side
        if( (x - m_center).dot(m_normalFrame.n) < 0.0f ) {
            pWi = 0.0f;
            return 0.0f;
        }

        Point3f s = m_objectToWorld * Point3f(-1.0f,-1.0f,0.0f);
        Vector3f ex = m_objectToWorld * Vector3f(2.0f,0.0f,0.0f);
        Vector3f ey = m_objectToWorld * Vector3f(0.0f,2.0f,0.0f);
        SphQuad squad;
        SphQuadInit(squad, s, ex, ey, x);
        y = SphQuadSample(squad, sampler->next1D(), sampler->next1D());

        Vector3f wo = (y-x).normalized();
        pWi = 1.0f / squad.S;
        normal = m_normalFrame.n;
        return wo;
    }

    Point3f SphQuadSample(const SphQuad& squad, float u, float v) const {
        // 1. compute 'cu'
        float au = u * squad.S + squad.k;
        float fu = (std::cos(au) * squad.b0 - squad.b1) / std::sin(au);
        float cu = 1.0f / std::sqrt(fu*fu + squad.b0sq) * (fu > 0 ? +1.0f : -1.0f);
        cu = clamp(cu, -1.0f, 1.0f); // avoid NaNs

        // 2. compute 'xu'
        float xu = -(cu * squad.z0) / std::sqrt(1.0f - cu*cu);
        xu = clamp(xu, squad.x0, squad.x1); // avoid Infs

        // 3. compute 'yv'
        float d = std::sqrt(xu*xu + squad.z0sq);
        float h0 = squad.y0 / std::sqrt(d*d + squad.y0sq);
        float h1 = squad.y1 / std::sqrt(d*d + squad.y1sq);
        float hv = h0 + v * (h1-h0), hv2 = hv*hv;
        float yv = (hv2 < 1-Epsilon) ? (hv*d)/std::sqrt(1 - hv2) : squad.y1;

        // 4. transform (xu,yv,z0) to world coords
        return (squad.o + xu*squad.x + yv*squad.y + squad.z0*squad.z);
    }

    void SphQuadInit(SphQuad& squad, const Point3f& s, const Vector3f& ex, const Vector3f& ey, const Point3f& o) const {
        squad.o = o;
        float exl = ex.norm(), eyl = ey.norm();

        // compute local reference system ’R’
        squad.x = ex / exl;
        squad.y = ey / eyl;
        squad.z = squad.x.cross(squad.y);

        // compute rectangle coords in local reference system
        Vector3f d = s - o;
        squad.z0 = d.dot(squad.z);

        // flip ’z’ to make it point against ’Q’
        if (squad.z0 > 0) {
            squad.z *= -1;
            squad.z0 *= -1;
        }
        squad.z0sq = squad.z0 * squad.z0;
        squad.x0 = d.dot(squad.x);
        squad.y0 = d.dot(squad.y);
        squad.x1 = squad.x0 + exl;
        squad.y1 = squad.y0 + eyl;
        squad.y0sq = squad.y0 * squad.y0;
        squad.y1sq = squad.y1 * squad.y1;

        // create vectors to four vertices
        Vector3f v00(squad.x0, squad.y0, squad.z0);
        Vector3f v01(squad.x0, squad.y1, squad.z0);
        Vector3f v10(squad.x1, squad.y0, squad.z0);
        Vector3f v11(squad.x1, squad.y1, squad.z0);

        // compute normals to edges
        Vector3f n0 = v00.cross(v10).normalized();
        Vector3f n1 = v10.cross(v11).normalized();
        Vector3f n2 = v11.cross(v01).normalized();
        Vector3f n3 = v01.cross(v00).normalized();

        // compute internal angles (gamma_i)
        float g0 = std::acos(-n0.dot(n1));
        float g1 = std::acos(-n1.dot(n2));
        float g2 = std::acos(-n2.dot(n3));
        float g3 = std::acos(-n3.dot(n0));

        // compute predefined constants
        squad.b0 = n0.z();
        squad.b1 = n2.z();
        squad.b0sq = squad.b0 * squad.b0;
        squad.k = 2*M_PI - g2 - g3;

        // compute solid angle from internal angles
        squad.S = g0 + g1 - squad.k;
        squad.S = std::max(0.0f, squad.S); // make sure it's not negative
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

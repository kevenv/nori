/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToTent(const Point2f &sample) {
    throw NoriException("Warp::squareToTent() is not yet implemented!");
}

float Warp::squareToTentPdf(const Point2f &p) {
    throw NoriException("Warp::squareToTentPdf() is not yet implemented!");
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    throw NoriException("Warp::squareToUniformDisk() is not yet implemented!");
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    throw NoriException("Warp::squareToUniformDiskPdf() is not yet implemented!");
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    throw NoriException("Warp::squareToUniformSphere() is not yet implemented!");
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    throw NoriException("Warp::squareToUniformSpherePdf() is not yet implemented!");
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
	float z = sample.x();
	float r = std::sqrt(std::max(0.0f, 1.0f - z * z));
	float phi = 2 * M_PI * sample.y();
	return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
	return v.z() > 0.0f ? INV_TWOPI : 0.0f;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
	Point2f d = concentricSampleDisk(sample);
	float z = std::sqrt(std::max(0.0f, 1 - d.x()*d.x() - d.y()*d.y()));
	return Vector3f(d.x(), d.y(), z);
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
	return v.z() > 0.0f ? v.z() * INV_PI : 0.0f;
}

Point2f Warp::concentricSampleDisk(const Point2f &sample) {
	// map uniform random numbers to [-1,1]^2
	Point2f uOffset = 2.0f * sample - Vector2f(1, 1);
	// handle degeneracy at the origin
	if (uOffset.x() == 0.0f && uOffset.y() == 0.0f) {
		return Point2f(0.0f, 0.0f);
	}
	// apply concentric mapping to point
	float theta, r;
	if (std::abs(uOffset.x()) > std::abs(uOffset.y())) {
		r = uOffset.x();
		theta = PI_OVER_FOUR * (uOffset.y() / uOffset.x());
	}
	else {
		r = uOffset.y();
		theta = PI_OVER_TWO - PI_OVER_FOUR * (uOffset.x() / uOffset.y());
	}
	return r * Point2f(std::cos(theta), std::sin(theta));
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    throw NoriException("Warp::squareToBeckmann() is not yet implemented!");
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    throw NoriException("Warp::squareToBeckmannPdf() is not yet implemented!");
}

NORI_NAMESPACE_END

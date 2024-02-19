#include "aleph0.h"
#include "primitives.h"

namespace aleph0
{

///////////////////////////////////////////////////////////////////////////////
float GetDistanceBetween(const Point3d& point1, const Point3d& point2)
{
	const float distance = (point2 - point1).GetNorm();

	return distance;
}

///////////////////////////////////////////////////////////////////////////////
float GetDistanceBetween(const Point3d& point, const Line& line)
{
	const Vector3d normalizedDirection = line.direction.GetNormalized();
	const Vector3d originToPoint = point - line.origin;
	const Vector3d closestPointOnLine = line.origin + originToPoint.ScalarProduct(normalizedDirection) * normalizedDirection;
	const Vector3d closestPointOnLineToPoint = point - closestPointOnLine;
	const float distance = closestPointOnLineToPoint.GetNorm();

	return distance;
}

///////////////////////////////////////////////////////////////////////////////
float GetDistanceBetween(const Point3d& point, const Plane& plane)
{
	const float normalNorm = plane.normal.GetNorm();
	const float normalTimesPoint = plane.normal.ScalarProduct(point);
	const float D = -plane.normal.ScalarProduct(plane.origin);
	const float distance = std::abs(normalTimesPoint + D) / normalNorm;

	return distance;
}

///////////////////////////////////////////////////////////////////////////////
Intersection IntersectLineWithLine(const Line& line1, const Line& line2)
{
	Intersection result;

	const Vector3d dir1CrossDir2 = line1.direction.VectorProduct(line2.direction);
	const float dir1CrossDir2NormSq = dir1CrossDir2.GetNormSq();

	if (dir1CrossDir2NormSq < aleph0::epsilon)
	{
		return result;
	}

	const Vector3d origin2ToOrigin1 = line2.origin - line1.origin;
	const Vector3d origin2toOrigin1CrossDir2 = origin2ToOrigin1.VectorProduct(line2.direction);
	const float t = origin2toOrigin1CrossDir2.ScalarProduct(dir1CrossDir2) / dir1CrossDir2NormSq;

	result.point = line1.origin + t * line1.direction;
	result.exists = true;

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Intersection IntersectLineWithPlane(const Line& line, const Plane& plane)
{
	Intersection result;

	const float directionTimesNormal = line.direction.ScalarProduct(plane.normal);

	if (std::abs(directionTimesNormal) < aleph0::epsilon)
	{
		return result;
	}

	const float originTimesNormal = line.origin.ScalarProduct(plane.normal);
	const float D = - plane.origin.ScalarProduct(plane.normal);
	const float t = (- D - originTimesNormal) / directionTimesNormal;

	result.point = line.origin + t * line.direction;

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Intersection IntersectRayWithTriangle_Naive(const Line& line, const Triangle& triangle)
{
	Intersection result;

	const Vector3d v0v1 = triangle.vertices[1] - triangle.vertices[0];
	const Vector3d v0v2 = triangle.vertices[2] - triangle.vertices[0];

	Plane plane;

	plane.normal = v0v1.VectorProduct(v0v2).Normalize();
	plane.origin = triangle.vertices[0];

	const float directionTimesNormal = line.direction.ScalarProduct(plane.normal);

	if (std::abs(directionTimesNormal) < aleph0::epsilon)
	{
		return result;
	}

	const float originTimesNormal = line.origin.ScalarProduct(plane.normal);
	const float D = -plane.origin.ScalarProduct(plane.normal);
	const float t = (-D - originTimesNormal) / directionTimesNormal;

	if (t < 0.f)
	{
		return result;
	}

	result.point = line.origin + t * line.direction;

	const Vector3d v1v2 = triangle.vertices[2] - triangle.vertices[1];
	const Vector3d v2v0 = triangle.vertices[0] - triangle.vertices[2];
	const Vector3d v0Q = result.point - triangle.vertices[0];
	const Vector3d v1Q = result.point - triangle.vertices[1];
	const Vector3d v2Q = result.point - triangle.vertices[2];

	const bool v0v1Check = (v0v1.VectorProduct(v0Q)).ScalarProduct(plane.normal) >= 0.f;
	const bool v1v2Check = (v1v2.VectorProduct(v1Q)).ScalarProduct(plane.normal) >= 0.f;
	const bool v2v0Check = (v2v0.VectorProduct(v2Q)).ScalarProduct(plane.normal) >= 0.f;

	result.exists = v0v1Check && v1v2Check && v2v0Check;

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Intersection IntersectRayWithTriangle_MollerTrumbore(const Line& line, const Triangle& triangle)
{
	Intersection result;

	const Vector3d E1 = triangle.vertices[1] - triangle.vertices[0];
	const Vector3d E2 = triangle.vertices[2] - triangle.vertices[0];

	const Vector3d P = line.direction.VectorProduct(E2);
	const float det = E1.ScalarProduct(P);

	if (std::abs(det) < aleph0::epsilon)
	{
		return result;
	}

	const float invDet = 1.f / det;
	const Vector3d T = line.origin - triangle.vertices[0];
	const float u = invDet * P.ScalarProduct(T);

	if (u < 0.f || u > 1.f)
	{
		return result;
	}

	const Vector3d Q = T.VectorProduct(E1);
	const float v = invDet * Q.ScalarProduct(line.direction);

	if (v < 0.f || u + v > 1.f)
	{
		return result;
	}

	const float t = invDet * Q.ScalarProduct(E2);

	if (t > aleph0::epsilon)
	{
		result.point = line.origin + t * line.direction;
	}

	return result;
}

} // namespace aleph0

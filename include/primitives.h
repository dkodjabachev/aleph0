#pragma once

#include "vector3d.h"

namespace aleph0
{

using Point3d = Vector3d;

///////////////////////////////////////////////////////////////////////////////
struct Line
{
	Point3d origin;
	Vector3d direction;
};

///////////////////////////////////////////////////////////////////////////////
struct Plane
{
	Point3d origin;
	Vector3d normal;
	float constant = 0.f;
};

///////////////////////////////////////////////////////////////////////////////
struct Triangle
{
	Point3d vertices[3];
};

///////////////////////////////////////////////////////////////////////////////
struct Intersection
{
	Point3d point;
	bool exists = false;
};

} // namespace aleph0

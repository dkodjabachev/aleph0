#pragma once

#include "aleph0.h"
#include "vector3d.h"

#include "matrix3x3.h"

namespace aleph0
{
///////////////////////////////////////////////////////////////////////////////
Vector3d::Vector3d()
{
	coordinates_[0] = 0.f;
	coordinates_[1] = 0.f;
	coordinates_[2] = 0.f;
}

///////////////////////////////////////////////////////////////////////////////
Vector3d::Vector3d(float x, float y, float z)
{
	coordinates_[0] = x;
	coordinates_[1] = y;
	coordinates_[2] = z;
}

///////////////////////////////////////////////////////////////////////////////
Vector3d::Vector3d(const Vector2d& other, float z)
{
	coordinates_[0] = other(0);
	coordinates_[1] = other(1);
	coordinates_[2] = z;
}

///////////////////////////////////////////////////////////////////////////////
Vector3d::Vector3d(const Vector3d& other)
{
	coordinates_[0] = other(0);
	coordinates_[1] = other(1);
	coordinates_[2] = other(2);
}

///////////////////////////////////////////////////////////////////////////////
const float& Vector3d::operator()(const int index) const
{
	assert(index >= 0 && index <= 2);

	return coordinates_[index];
}

///////////////////////////////////////////////////////////////////////////////
float& Vector3d::operator()(const int index)
{
	assert(index >= 0 && index <= 2);

	return coordinates_[index];
}

///////////////////////////////////////////////////////////////////////////////
Vector3d Vector3d::operator+(const Vector3d& other) const
{
	Vector3d result;

	result(0) = coordinates_[0] + other(0);
	result(1) = coordinates_[1] + other(1);
	result(2) = coordinates_[2] + other(2);

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Vector3d& Vector3d::operator+=(const Vector3d& other)
{
	coordinates_[0] += other(0);
	coordinates_[1] += other(1);
	coordinates_[2] += other(2);

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Vector3d Vector3d::operator-(const Vector3d& other) const
{
	Vector3d result;

	result(0) = coordinates_[0] - other(0);
	result(1) = coordinates_[1] - other(1);
	result(2) = coordinates_[2] - other(2);

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Vector3d& Vector3d::operator-=(const Vector3d& other)
{
	coordinates_[0] -= other(0);
	coordinates_[1] -= other(1);
	coordinates_[2] -= other(2);

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Vector3d Vector3d::operator*(const float scalar) const
{
	Vector3d result;

	result(0) = coordinates_[0] * scalar;
	result(1) = coordinates_[1] * scalar;
	result(2) = coordinates_[2] * scalar;

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Vector3d& Vector3d::operator*=(const float scalar)
{
	coordinates_[0] *= scalar;
	coordinates_[1] *= scalar;
	coordinates_[2] *= scalar;

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
bool Vector3d::operator==(const Vector3d& other) const
{
	return coordinates_[0] == other(0)
		&& coordinates_[1] == other(1)
		&& coordinates_[2] == other(2);
}

///////////////////////////////////////////////////////////////////////////////
bool Vector3d::operator!=(const Vector3d& other) const
{
	return !(*this == other);
}

///////////////////////////////////////////////////////////////////////////////
float Vector3d::ScalarProduct(const Vector3d& other) const
{
	return coordinates_[0] * other(0) + coordinates_[1] * other(1) + coordinates_[2] * other(2);
}

///////////////////////////////////////////////////////////////////////////////
Vector3d Vector3d::VectorProduct(const Vector3d& other) const
{
	Vector3d result;

	result(0) = coordinates_[1] * other(2) - coordinates_[2] * other(1);
	result(1) = coordinates_[2] * other(0) - coordinates_[0] * other(2);
	result(2) = coordinates_[0] * other(1) - coordinates_[1] * other(0);

	return result;
}

///////////////////////////////////////////////////////////////////////////////
float Vector3d::GetNorm() const
{
	return std::sqrt(ScalarProduct(*this));
}

///////////////////////////////////////////////////////////////////////////////
float Vector3d::GetNormSq() const
{
	return ScalarProduct(*this);
}

///////////////////////////////////////////////////////////////////////////////
Vector3d& Vector3d::Normalize()
{
	const float normSq = GetNormSq();

	assert(normSq > aleph0::epsilon);

	*this *= (1 / std::sqrt(normSq));

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Vector3d Vector3d::GetNormalized() const
{
	Vector3d result(*this);

	result.Normalize();

	return result;
}

///////////////////////////////////////////////////////////////////////////////
bool Vector3d::IsZero() const
{
	const float normSq = GetNormSq();

	return normSq < aleph0::epsilon;
}

///////////////////////////////////////////////////////////////////////////////
bool Vector3d::IsOrthogonalTo(const Vector3d& other) const
{
	return std::abs(ScalarProduct(other)) < aleph0::epsilon;
}

///////////////////////////////////////////////////////////////////////////////
bool Vector3d::IsCollinearTo(const Vector3d& other) const
{
	return VectorProduct(other).GetNorm() < aleph0::epsilon;
}

///////////////////////////////////////////////////////////////////////////////
std::string Vector3d::ToString() const
{
	const std::string result =
		"(" +
		std::to_string(coordinates_[0]) + "," +
		std::to_string(coordinates_[1]) + "," +
		std::to_string(coordinates_[2]) + ")";

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Vector3d operator*(const float scalar, const Vector3d& vector)
{
	return vector * scalar;
}

} // namespace aleph0

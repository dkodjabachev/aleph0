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
Vector3d& Vector3d::RotateAroundOx(const float theta)
{
	Matrix3x3 rotationMatrix;

	rotationMatrix.SetRotationAroundOx(theta);

	*this = rotationMatrix * (*this);

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Vector3d& Vector3d::RotateAroundOy(const float theta)
{
	Matrix3x3 rotationMatrix;

	rotationMatrix.SetRotationAroundOy(theta);

	*this = rotationMatrix * (*this);

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Vector3d& Vector3d::RotateAroundOz(const float theta)
{
	Matrix3x3 rotationMatrix;

	rotationMatrix.SetRotationAroundOz(theta);

	*this = rotationMatrix * (*this);

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Vector3d& Vector3d::RotateAroundAxis(const Vector3d& axis, const float theta)
{
	Matrix3x3 rotationMatrix;

	rotationMatrix.SetRotationAroundAxis(axis, theta);

	*this = rotationMatrix * (*this);

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Vector3d& Vector3d::RotateAroundAxisByRodrigues(const Vector3d& axis, const float theta)
{
	const float cosTheta = std::cosf(theta);
	const float sinTheta = std::sinf(theta);

	*this = *this * cosTheta + axis.VectorProduct(*this) * sinTheta + axis * (axis.ScalarProduct(*this)) * (1 - cosTheta);

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
float Vector3d::GetAngleTo(const Vector3d& vector) const
{
	assert(GetNormSq() > aleph0::epsilon);
	assert(vector.GetNormSq() > aleph0::epsilon);

	const float scalarProduct = ScalarProduct(vector);
	const float theta = std::acosf(scalarProduct / (GetNorm() * vector.GetNorm()));

	return theta;
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

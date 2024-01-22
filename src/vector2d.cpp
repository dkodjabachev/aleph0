#pragma once

#include "aleph0.h"
#include "vector2d.h"

namespace aleph0
{

///////////////////////////////////////////////////////////////////////////////
Vector2d::Vector2d()
{
	coordinates_[0] = 0.f;
	coordinates_[1] = 0.f;
}

///////////////////////////////////////////////////////////////////////////////
Vector2d::Vector2d(float x, float y)
{
	coordinates_[0] = x;
	coordinates_[1] = y;
}

///////////////////////////////////////////////////////////////////////////////
Vector2d::Vector2d(const Vector2d& other)
{
	coordinates_[0] = other(0);
	coordinates_[1] = other(1);
}

///////////////////////////////////////////////////////////////////////////////
const float& Vector2d::operator()(const int index) const
{
	assert(index >= 0 && index <= 1);

	return coordinates_[index];
}

///////////////////////////////////////////////////////////////////////////////
float& Vector2d::operator()(const int index)
{
	assert(index >= 0 && index <= 1);

	return coordinates_[index];
}

///////////////////////////////////////////////////////////////////////////////
Vector2d Vector2d::operator+(const Vector2d& other) const
{
	Vector2d result;

	result(0) = coordinates_[0] + other(0);
	result(1) = coordinates_[1] + other(1);

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Vector2d& Vector2d::operator+=(const Vector2d& other)
{
	coordinates_[0] += other(0);
	coordinates_[1] += other(1);

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Vector2d Vector2d::operator-(const Vector2d& other) const
{
	Vector2d result;

	result(0) = coordinates_[0] - other(0);
	result(1) = coordinates_[1] - other(1);

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Vector2d& Vector2d::operator-=(const Vector2d& other)
{
	coordinates_[0] -= other(0);
	coordinates_[1] -= other(1);

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Vector2d Vector2d::operator*(const float scalar) const
{
	Vector2d result;

	result(0) = coordinates_[0] * scalar;
	result(1) = coordinates_[1] * scalar;

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Vector2d& Vector2d::operator*=(const float scalar)
{
	coordinates_[0] *= scalar;
	coordinates_[1] *= scalar;

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
bool Vector2d::operator==(const Vector2d& other) const
{
	return coordinates_[0] == other(0) && coordinates_[1] == other(1);
}

///////////////////////////////////////////////////////////////////////////////
bool Vector2d::operator!=(const Vector2d& other) const
{
	return !(*this == other);
}

///////////////////////////////////////////////////////////////////////////////
float Vector2d::ScalarProduct(const Vector2d& other) const
{
	return coordinates_[0] * other(0) + coordinates_[1] * other(1);
}

///////////////////////////////////////////////////////////////////////////////
float Vector2d::GetNorm() const
{
	return std::sqrt(ScalarProduct(*this));
}

///////////////////////////////////////////////////////////////////////////////
float Vector2d::GetNormSq() const
{
	return ScalarProduct(*this);
}

///////////////////////////////////////////////////////////////////////////////
Vector2d& Vector2d::Normalize()
{
	const float normSq = GetNormSq();

	assert(normSq > aleph0::epsilon);

	*this *= (1 / std::sqrt(normSq));

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Vector2d Vector2d::GetNormalized() const
{
	Vector2d result(*this);

	result.Normalize();

	return result;
}

///////////////////////////////////////////////////////////////////////////////
bool Vector2d::IsZero() const
{
	const float normSq = GetNormSq();

	return normSq < aleph0::epsilon;
}

///////////////////////////////////////////////////////////////////////////////
bool Vector2d::IsOrthogonalTo(const Vector2d& other) const
{
	return std::abs(ScalarProduct(other)) < aleph0::epsilon;
}

///////////////////////////////////////////////////////////////////////////////
bool Vector2d::IsCollinearTo(const Vector2d& other) const
{
	const float det = coordinates_[0] * other(1) - coordinates_[1] * other(0);

	return std::abs(det) < aleph0::epsilon;
}

///////////////////////////////////////////////////////////////////////////////
std::string Vector2d::ToString() const
{
	const std::string result =
		"(" +
		std::to_string(coordinates_[0]) + "," +
		std::to_string(coordinates_[1]) + ")";

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Vector2d operator*(const float scalar, const Vector2d& vector)
{
	return vector * scalar;
}

} // namespace aleph0

#pragma once

// Aleph0 library headers
#include "aleph0.h"

namespace aleph0
{

class Vector2d
{
public:
	///////////////////////////////////////////////////////////////////////////////
	Vector2d()
	{
		coordinates_[0] = 0.f;
		coordinates_[1] = 0.f;
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector2d(float x, float y)
	{
		coordinates_[0] = x;
		coordinates_[1] = y;
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector2d(const Vector2d& other)
	{
		coordinates_[0] = other(0);
		coordinates_[1] = other(1);
	}

public:
	///////////////////////////////////////////////////////////////////////////////
	const float& operator()(const int index) const
	{
		assert(index >= 0 && index <= 1);

		return coordinates_[index];
	}

	///////////////////////////////////////////////////////////////////////////////
	float& operator()(const int index)
	{
		assert(index >= 0 && index <= 1);

		return coordinates_[index];
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector2d operator+(const Vector2d& other) const
	{
		Vector2d result;

		result(0) = coordinates_[0] + other(0);
		result(1) = coordinates_[1] + other(1);

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector2d& operator+=(const Vector2d& other)
	{
		coordinates_[0] += other(0);
		coordinates_[1] += other(1);

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector2d operator-(const Vector2d& other) const
	{
		Vector2d result;

		result(0) = coordinates_[0] - other(0);
		result(1) = coordinates_[1] - other(1);

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector2d& operator-=(const Vector2d& other)
	{
		coordinates_[0] -= other(0);
		coordinates_[1] -= other(1);

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector2d operator*(const float scalar) const
	{
		Vector2d result;

		result(0) = coordinates_[0] * scalar;
		result(1) = coordinates_[1] * scalar;

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector2d& operator*=(const float scalar)
	{
		coordinates_[0] *= scalar;
		coordinates_[1] *= scalar;

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	bool operator==(const Vector2d& other) const
	{
		return coordinates_[0] == other(0) && coordinates_[1] == other(1);
	}

	///////////////////////////////////////////////////////////////////////////////
	bool operator!=(const Vector2d& other) const
	{
		return !(*this == other);
	}

	///////////////////////////////////////////////////////////////////////////////
	float ScalarProduct(const Vector2d& other) const
	{
		return coordinates_[0] * other(0) + coordinates_[1] * other(1);
	}

	///////////////////////////////////////////////////////////////////////////////
	float GetNorm() const
	{
		return std::sqrt(ScalarProduct(*this));
	}

	///////////////////////////////////////////////////////////////////////////////
	float GetNormSq() const
	{
		return ScalarProduct(*this);
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector2d& Normalize()
	{
		const float normSq = GetNormSq();

		assert(normSq > aleph0::epsilon);

		*this *= (1 / std::sqrt(normSq));

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector2d GetNormalized() const
	{
		Vector2d result(*this);

		result.Normalize();

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	bool IsZero() const
	{
		const float normSq = GetNormSq();

		return normSq < aleph0::epsilon;
	}

	///////////////////////////////////////////////////////////////////////////////
	bool IsOrthogonalTo(const Vector2d& other) const
	{
		return std::abs(ScalarProduct(other)) < aleph0::epsilon;
	}

	///////////////////////////////////////////////////////////////////////////////
	bool IsCollinearTo(const Vector2d& other) const
	{
		const float det = coordinates_[0] * other(1) - coordinates_[1] * other(0);

		return std::abs(det) < aleph0::epsilon;
	}

	///////////////////////////////////////////////////////////////////////////////
	std::string ToString() const
	{
		const std::string result =
			"(" +
			std::to_string(coordinates_[0]) + "," +
			std::to_string(coordinates_[1]) + ")";

		return result;
	}

private:
	float coordinates_[2];
}; // class Vector2d

///////////////////////////////////////////////////////////////////////////////
Vector2d operator*(const float scalar, const Vector2d& vector)
{
	return vector * scalar;
}

} // namespace aleph0

#pragma once

// Aleph0 library headers
#include "aleph0.h"
#include "vector2d.h"

namespace aleph0
{

class Vector3d
{
public:
	///////////////////////////////////////////////////////////////////////////////
	Vector3d()
	{
		coordinates_[0] = 0.f;
		coordinates_[1] = 0.f;
		coordinates_[2] = 0.f;
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector3d(float x, float y, float z = 0.f)
	{
		coordinates_[0] = x;
		coordinates_[1] = y;
		coordinates_[2] = z;
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector3d(const Vector2d& other, float z = 0.f)
	{
		coordinates_[0] = other(0);
		coordinates_[1] = other(1);
		coordinates_[2] = z;
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector3d(const Vector3d& other)
	{
		coordinates_[0] = other(0);
		coordinates_[1] = other(1);
		coordinates_[2] = other(2);
	}

public:
	///////////////////////////////////////////////////////////////////////////////
	const float& operator()(const int index) const
	{
		assert(index >= 0 && index <= 2);

		return coordinates_[index];
	}

	///////////////////////////////////////////////////////////////////////////////
	float& operator()(const int index)
	{
		assert(index >= 0 && index <= 2);

		return coordinates_[index];
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector3d operator+(const Vector3d& other) const
	{
		Vector3d result;

		result(0) = coordinates_[0] + other(0);
		result(1) = coordinates_[1] + other(1);
		result(2) = coordinates_[2] + other(2);

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector3d& operator+=(const Vector3d& other)
	{
		coordinates_[0] += other(0);
		coordinates_[1] += other(1);
		coordinates_[2] += other(2);

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector3d operator-(const Vector3d& other) const
	{
		Vector3d result;

		result(0) = coordinates_[0] - other(0);
		result(1) = coordinates_[1] - other(1);
		result(2) = coordinates_[2] - other(2);

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector3d& operator-=(const Vector3d& other)
	{
		coordinates_[0] -= other(0);
		coordinates_[1] -= other(1);
		coordinates_[2] -= other(2);

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector3d operator*(const float scalar) const
	{
		Vector3d result;

		result(0) = coordinates_[0] * scalar;
		result(1) = coordinates_[1] * scalar;
		result(2) = coordinates_[2] * scalar;

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector3d& operator*=(const float scalar)
	{
		coordinates_[0] *= scalar;
		coordinates_[1] *= scalar;
		coordinates_[2] *= scalar;

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	bool operator==(const Vector3d& other) const
	{
		return coordinates_[0] == other(0)
			&& coordinates_[1] == other(1)
			&& coordinates_[2] == other(2);
	}

	///////////////////////////////////////////////////////////////////////////////
	bool operator!=(const Vector3d& other) const
	{
		return !(*this == other);
	}

	///////////////////////////////////////////////////////////////////////////////
	float ScalarProduct(const Vector3d& other) const
	{
		return coordinates_[0] * other(0) + coordinates_[1] * other(1) + coordinates_[2] * other(2);
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector3d VectorProduct(const Vector3d& other) const
	{
		Vector3d result;

		result(0) = coordinates_[1] * other(2) - coordinates_[2] * other(1);
		result(1) = coordinates_[2] * other(0) - coordinates_[0] * other(2);
		result(2) = coordinates_[0] * other(1) - coordinates_[1] * other(0);

		return result;
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
	Vector3d& Normalize()
	{
		const float normSq = GetNormSq();

		assert(normSq > aleph0::epsilon);

		*this *= (1 / std::sqrt(normSq));

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Vector3d GetNormalized() const
	{
		Vector3d result(*this);

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
	bool IsOrthogonalTo(const Vector3d& other) const
	{
		return std::abs(ScalarProduct(other)) < aleph0::epsilon;
	}

	///////////////////////////////////////////////////////////////////////////////
	bool IsCollinearTo(const Vector3d& other) const
	{
		return VectorProduct(other).GetNorm() < aleph0::epsilon;
	}

	///////////////////////////////////////////////////////////////////////////////
	std::string ToString() const
	{
		const std::string result =
			"(" +
			std::to_string(coordinates_[0]) + "," +
			std::to_string(coordinates_[1]) + "," +
			std::to_string(coordinates_[2]) + ")";

		return result;
	}

private:
	float coordinates_[3];
}; // class Vector3d

///////////////////////////////////////////////////////////////////////////////
Vector3d operator*(const float scalar, const Vector3d& vector)
{
	return vector * scalar;
}

} // namespace aleph0

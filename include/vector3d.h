#pragma once

#include "aleph0.h"
#include "vector2d.h"

namespace aleph0
{

class Vector3d
{
public:
	///////////////////////////////////////////////////////////////////////////////
	Vector3d();

	///////////////////////////////////////////////////////////////////////////////
	Vector3d(float x, float y, float z = 0.f);

	///////////////////////////////////////////////////////////////////////////////
	Vector3d(const Vector2d& other, float z = 0.f);

	///////////////////////////////////////////////////////////////////////////////
	Vector3d(const Vector3d& other);

public:
	///////////////////////////////////////////////////////////////////////////////
	const float& operator()(const int index) const;

	///////////////////////////////////////////////////////////////////////////////
	float& operator()(const int index);

	///////////////////////////////////////////////////////////////////////////////
	Vector3d operator+(const Vector3d& other) const;

	///////////////////////////////////////////////////////////////////////////////
	Vector3d& operator+=(const Vector3d& other);

	///////////////////////////////////////////////////////////////////////////////
	Vector3d operator-(const Vector3d& other) const;

	///////////////////////////////////////////////////////////////////////////////
	Vector3d& operator-=(const Vector3d& other);

	///////////////////////////////////////////////////////////////////////////////
	Vector3d operator*(const float scalar) const;

	///////////////////////////////////////////////////////////////////////////////
	Vector3d& operator*=(const float scalar);

	///////////////////////////////////////////////////////////////////////////////
	bool operator==(const Vector3d& other) const;

	///////////////////////////////////////////////////////////////////////////////
	bool operator!=(const Vector3d& other) const;

	///////////////////////////////////////////////////////////////////////////////
	float ScalarProduct(const Vector3d& other) const;

	///////////////////////////////////////////////////////////////////////////////
	Vector3d VectorProduct(const Vector3d& other) const;

	///////////////////////////////////////////////////////////////////////////////
	float GetNorm() const;

	///////////////////////////////////////////////////////////////////////////////
	float GetNormSq() const;

	///////////////////////////////////////////////////////////////////////////////
	Vector3d& Normalize();

	///////////////////////////////////////////////////////////////////////////////
	Vector3d GetNormalized() const;

	///////////////////////////////////////////////////////////////////////////////
	bool IsZero() const;

	///////////////////////////////////////////////////////////////////////////////
	bool IsOrthogonalTo(const Vector3d& other) const;

	///////////////////////////////////////////////////////////////////////////////
	bool IsCollinearTo(const Vector3d& other) const;

	///////////////////////////////////////////////////////////////////////////////
	Vector3d& RotateAroundOx(const float theta);

	///////////////////////////////////////////////////////////////////////////////
	Vector3d& RotateAroundOy(const float theta);

	///////////////////////////////////////////////////////////////////////////////
	Vector3d& RotateAroundOz(const float theta);

	///////////////////////////////////////////////////////////////////////////////
	Vector3d& RotateAroundAxis(const Vector3d& axis, const float theta);

	///////////////////////////////////////////////////////////////////////////////
	Vector3d& RotateAroundAxisByRodrigues(const Vector3d& axis, const float theta);

	///////////////////////////////////////////////////////////////////////////////
	float GetAngleTo(const Vector3d& vector) const;

	///////////////////////////////////////////////////////////////////////////////
	std::string ToString() const;

private:
	float coordinates_[3];
}; // class Vector3d

///////////////////////////////////////////////////////////////////////////////
Vector3d operator*(const float scalar, const Vector3d& vector);

} // namespace aleph0

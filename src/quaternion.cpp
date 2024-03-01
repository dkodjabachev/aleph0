#include "aleph0.h"
#include "quaternion.h"

#include "vector3d.h"

namespace aleph0
{

///////////////////////////////////////////////////////////////////////////////
Quaternion::Quaternion()
{
	w = 0.f;
	x = 0.f;
	y = 0.f;
	z = 0.f;
}

///////////////////////////////////////////////////////////////////////////////
Quaternion::Quaternion(float w_, float x_, float y_, float z_)
{
	w = w_;
	x = x_;
	y = y_;
	z = z_;
}

///////////////////////////////////////////////////////////////////////////////
Quaternion::Quaternion(const Quaternion& other)
{
	w = other.w;
	x = other.x;
	y = other.y;
	z = other.z;
}

///////////////////////////////////////////////////////////////////////////////
Quaternion& Quaternion::SetRotationAroundAxis(const Vector3d& axis, const float theta)
{
	const float halfTheta = 0.5f * theta;
	const float cosHalfTheta = std::cos(halfTheta);
	const float sinHalfTheta = std::sin(halfTheta);

	w = cosHalfTheta;
	x = sinHalfTheta * axis(0);
	y = sinHalfTheta * axis(1);
	z = sinHalfTheta * axis(2);

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Quaternion Quaternion::operator+(const Quaternion& other) const
{
	Quaternion result;

	result.w = w + other.w;
	result.x = x + other.x;
	result.y = y + other.y;
	result.z = z + other.z;

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Quaternion& Quaternion::operator+=(const Quaternion& other)
{
	*this = *this + other;

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Quaternion Quaternion::operator-(const Quaternion& other) const
{
	Quaternion result;

	result.w = w - other.w;
	result.x = x - other.x;
	result.y = y - other.y;
	result.z = z - other.z;

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Quaternion& Quaternion::operator-=(const Quaternion& other)
{
	*this = *this - other;

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Quaternion Quaternion::operator*(const Quaternion& other) const
{
	Quaternion result;

	result.w = w * other.w - x * other.x - y * other.y - z * other.z;
	result.x = (other.x * w) + (other.w * x) + (other.z * y) - (other.y * z);
	result.y = (other.y * w) + (other.w * y) + (other.x * z) - (other.z * x);
	result.z = (other.z * w) + (other.w * z) + (other.y * x) - (other.x * y);

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Quaternion& Quaternion::operator*=(const Quaternion& other)
{
	*this = *this * other;

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Quaternion Quaternion::operator*(const float scalar) const
{
	Quaternion result;

	result.w *= scalar;
	result.x *= scalar;
	result.y *= scalar;
	result.z *= scalar;

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Quaternion& Quaternion::operator*=(const float scalar)
{
	*this = *this * scalar;

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Quaternion Quaternion::GetConjugate() const
{
	Quaternion result(*this);

	result.x *= -1;
	result.y *= -1;
	result.z *= -1;

	return result;
}

///////////////////////////////////////////////////////////////////////////////
float Quaternion::GetNorm() const
{
	const float norm = std::sqrt(w*w + x*x + y*y + z*z);

	return norm;
}

///////////////////////////////////////////////////////////////////////////////
float Quaternion::GetNormSq() const
{
	const float norm = w * w + x * x + y * y + z * z;

	return norm;
}

///////////////////////////////////////////////////////////////////////////////
Quaternion& Quaternion::Normalize()
{
	const float normSq = GetNormSq();

	assert(normSq > aleph0::epsilon);

	*this *= (1.f / std::sqrt(normSq));

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Quaternion Quaternion::GetNormalized() const
{
	Quaternion result(*this);

	result.Normalize();

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Quaternion Quaternion::GetInverse() const
{
	const float normSq = GetNormSq();

	assert(normSq > aleph0::epsilon);

	Quaternion result = GetConjugate() * (1.f / normSq);

	return result;
}

} // namespace aleph0

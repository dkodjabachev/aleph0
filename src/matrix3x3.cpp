#include "aleph0.h"
#include "matrix3x3.h"

#include "matrix2x2.h"
#include "vector3d.h"

namespace aleph0
{

///////////////////////////////////////////////////////////////////////////////
Matrix3x3::Matrix3x3()
{
	SetZeroMatrix();
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3::Matrix3x3(const Matrix3x3& other)
{
	memcpy(data_.elements_, other.data_.elements_, 9 * sizeof(float));
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3::Matrix3x3(const float elements[9])
{
	memcpy(data_.elements_, elements, 9 * sizeof(float));
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3& Matrix3x3::SetZeroMatrix()
{
	memset(data_.elements_, 0, 9 * sizeof(float));

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3& Matrix3x3::SetIdentityMatrix()
{
	SetZeroMatrix();

	m(0, 0) = 1.f;
	m(1, 1) = 1.f;
	m(2, 2) = 1.f;

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3& Matrix3x3::SetRotationAroundOx(const float theta)
{
	SetIdentityMatrix();

	const float cosTheta = std::cosf(theta);
	const float sinTheta = std::sinf(theta);

	m(1, 1) = cosTheta;
	m(1, 2) = -sinTheta;
	m(2, 1) = sinTheta;
	m(2, 2) = cosTheta;

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3& Matrix3x3::SetRotationAroundOy(const float theta)
{
	SetIdentityMatrix();

	const float cosTheta = std::cosf(theta);
	const float sinTheta = std::sinf(theta);

	m(0, 0) = cosTheta;
	m(0, 2) = -sinTheta;
	m(2, 0) = sinTheta;
	m(2, 2) = cosTheta;

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3& Matrix3x3::SetRotationAroundOz(const float theta)
{
	SetIdentityMatrix();

	const float cosTheta = std::cosf(theta);
	const float sinTheta = std::sinf(theta);

	m(0, 0) = cosTheta;
	m(0, 1) = -sinTheta;
	m(1, 0) = sinTheta;
	m(1, 1) = cosTheta;

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3& Matrix3x3::SetRotationAroundAxis(const Vector3d& axis, const float theta)
{
	assert(std::abs(axis.GetNormSq() - 1.f) > aleph0::epsilon);

	SetIdentityMatrix();

	const float cosTheta = std::cosf(theta);
	const float sinTheta = std::sinf(theta);

	m(0, 0) = cosTheta + axis(0) * axis(0) * (1 - cosTheta);
	m(0, 1) = axis(0) * axis(1) * (1 - cosTheta) - axis(2) * sinTheta;
	m(0, 2) = axis(0) * axis(2) * (1 - cosTheta) + axis(1) * sinTheta;

	m(1, 0) = axis(0) * axis(1) * (1 - cosTheta) + axis(2) * sinTheta;
	m(1, 1) = cosTheta + axis(1) * axis(1) * (1 - cosTheta);
	m(1, 2) = axis(1) * axis(2) * (1 - cosTheta) - axis(0) * sinTheta;

	m(2, 0) = axis(0) * axis(2) * (1 - cosTheta) - axis(1) * sinTheta;
	m(2, 1) = axis(1) * axis(2) * (1 - cosTheta) + axis(0) * sinTheta;
	m(2, 2) = cosTheta + axis(2) * axis(2) * (1 - cosTheta);

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
const float& Matrix3x3::operator()(const int row_index, const int col_index) const
{
	return m(row_index, col_index);
}

///////////////////////////////////////////////////////////////////////////////
float& Matrix3x3::operator()(const int row_index, const int col_index)
{
	return m(row_index, col_index);
}

///////////////////////////////////////////////////////////////////////////////
const float& Matrix3x3::operator()(const int index) const
{
	assert(index >= 0 && index < 9);

	return data_.elements_[index];
}

///////////////////////////////////////////////////////////////////////////////
float& Matrix3x3::operator()(const int index)
{
	assert(index >= 0 && index < 9);

	return data_.elements_[index];
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3& Matrix3x3::operator=(const Matrix3x3& other)
{
	memcpy(data_.elements_, other.data_.elements_, 9 * sizeof(float));

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3 Matrix3x3::operator+(const Matrix3x3& other) const
{
	Matrix3x3 result;

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			result(i, j) = m(i, j) + other(i, j);
		}
	}

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3& Matrix3x3::operator+=(const Matrix3x3& other)
{
	*this = *this + other;

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3 Matrix3x3::operator-(const Matrix3x3& other) const
{
	Matrix3x3 result;

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			result(i, j) = m(i, j) - other(i, j);
		}
	}

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3& Matrix3x3::operator-=(const Matrix3x3& other)
{
	*this = *this - other;

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3 Matrix3x3::operator*(const float scalar) const
{
	Matrix3x3 result;

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			result(i, j) = m(i, j) * scalar;
		}
	}

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3& Matrix3x3::operator*=(const float scalar)
{
	*this = *this * scalar;

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3 Matrix3x3::operator*(const Matrix3x3& other) const
{
	return NaiveMatrixProduct(*this, other);
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3& Matrix3x3::operator*=(const Matrix3x3& other)
{
	*this = NaiveMatrixProduct(*this, other);

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Vector3d Matrix3x3::operator*(const Vector3d& vector)
{
	Vector3d result;

	result(0) = m(0, 0) * vector(0) + m(0, 1) * vector(1) + m(0, 2) * vector(2);
	result(1) = m(1, 0) * vector(0) + m(1, 1) * vector(1) + m(1, 2) * vector(2);
	result(2) = m(2, 0) * vector(0) + m(2, 1) * vector(1) + m(2, 2) * vector(2);

	return result;
}

///////////////////////////////////////////////////////////////////////////////
bool Matrix3x3::operator==(const Matrix3x3& other) const
{
	return IsApproximatelyEqual(other, 0.0f);
}

///////////////////////////////////////////////////////////////////////////////
bool Matrix3x3::operator!=(const Matrix3x3& other) const
{
	return !(*this == other);
}

///////////////////////////////////////////////////////////////////////////////
bool Matrix3x3::IsApproximatelyEqual(const Matrix3x3& other, const float epsilon) const
{
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			if (std::abs(m(i, j) - other(i, j)) > epsilon)
			{
				return false;
			}
		}
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3 Matrix3x3::NaiveMatrixProduct(const Matrix3x3& A, const Matrix3x3& B)
{
	Matrix3x3 result;

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			for (int k = 0; k < 3; ++k)
			{
				result(i, j) += A(i, k) * B(k, j);
			}
		}
	}

	return result;
}

//////////////////////////////////////////////////////////////////////////////
float Matrix3x3::Trace() const
{
	const float trace = m(0, 0) + m(1, 1) + m(2, 2);

	return trace;
}

///////////////////////////////////////////////////////////////////////////////
float Matrix3x3::Determinant() const
{
	const float det =
		m(0, 0) * m(1, 1) * m(2, 2) +
		m(0, 1) * m(1, 2) * m(2, 0) +
		m(0, 2) * m(1, 0) * m(2, 1) -
		m(0, 2) * m(1, 1) * m(2, 0) -
		m(0, 0) * m(1, 2) * m(2, 1) -
		m(0, 1) * m(1, 0) * m(2, 2);

	return det;
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3& Matrix3x3::Transpose()
{
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			if (i < j)
			{
				std::swap(m(i, j), m(j, i));
			}
		}
	}

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3 Matrix3x3::Transposed() const
{
	Matrix3x3 result(*this);

	result.Transpose();

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3& Matrix3x3::Invert()
{
	const float det = Determinant();

	assert(det != 0);

	Matrix3x3 cofactor_matrix;

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			cofactor_matrix(i, j) = std::powf(-1.f, (float)(i + j)) * Submatrix(i, j).Determinant();
		}
	}

	*this = cofactor_matrix.Transpose() * (1 / det);

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3 Matrix3x3::Inverse() const
{
	Matrix3x3 result(*this);

	result.Invert();

	return result;
}

///////////////////////////////////////////////////////////////////////////////
void Matrix3x3::GetRotationAxisAngle(Vector3d& axis, float& theta) const
{
	const float x = m(2, 1) - m(1, 2);
	const float y = m(0, 2) - m(2, 0);
	const float z = m(1, 0) - m(0, 1);
	const float r = std::sqrt(x * x + y * y + z * z);
	const float trace = Trace();
	
	theta = std::atan2(r, trace - 1);
	axis(0) = x;
	axis(1) = y;
	axis(2) = z;
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2 Matrix3x3::Submatrix(int skip_row_index, int skip_col_index) const
{
	Matrix2x2 result;
	int index = 0;

	for (int i = 0; i < 3; ++i)
	{
		if (i == skip_row_index)
		{
			continue;
		}

		for (int j = 0; j < 3; ++j)
		{
			if (j == skip_col_index)
			{
				continue;
			}

			result(index++) = m(i, j);
		}
	}

	return result;
}

///////////////////////////////////////////////////////////////////////////////
std::string Matrix3x3::ToString() const
{
	std::string result;

	for (int i = 0; i < 3; ++i)
	{
		result += "\n(";

		for (int j = 0; j < 3; ++j)
		{
			result += std::to_string(m(i, j));

			if (j + 1 < 3)
			{
				result += ",";
			}
			else
			{
				result += ")";
			}
		}
	}

	return result;
}

///////////////////////////////////////////////////////////////////////////////
const float& Matrix3x3::m(const int row_index, const int col_index) const
{
	assert(row_index >= 0 && row_index < 3);
	assert(col_index >= 0 && col_index < 3);

	return data_.array_[row_index][col_index];
}

///////////////////////////////////////////////////////////////////////////////
float& Matrix3x3::m(const int row_index, const int col_index)
{
	assert(row_index >= 0 && row_index < 3);
	assert(col_index >= 0 && col_index < 3);

	return data_.array_[row_index][col_index];
}

///////////////////////////////////////////////////////////////////////////////
Matrix3x3 operator*(const float scalar, const Matrix3x3& matrix)
{
	return matrix * scalar;
}

} // namespace aleph0

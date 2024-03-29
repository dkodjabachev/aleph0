#include "aleph0.h"
#include "matrix2x2.h"

#include "vector2d.h"

namespace aleph0
{

///////////////////////////////////////////////////////////////////////////////
Matrix2x2::Matrix2x2()
{
	SetZeroMatrix();
}

Matrix2x2::Matrix2x2(const Matrix2x2& other)
{
	memcpy(data_.elements_, other.data_.elements_, 4 * sizeof(float));
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2::Matrix2x2(float a, float b, float c, float d)
{
	m(0, 0) = a;
	m(0, 1) = b;
	m(1, 0) = c;
	m(1, 1) = d;
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2& Matrix2x2::SetZeroMatrix()
{
	memset(data_.elements_, 0, 4 * sizeof(float));

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2& Matrix2x2::SetIdentityMatrix()
{
	m(0, 0) = 1.f;
	m(0, 1) = 0.f;
	m(1, 0) = 0.f;
	m(1, 1) = 1.f;

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2& Matrix2x2::SetRotationAroundOrigin(const float theta)
{
	const float cosTheta = std::cosf(theta);
	const float sinTheta = std::sinf(theta);

	m(0, 0) = cosTheta;
	m(0, 1) = -sinTheta;
	m(1, 0) = sinTheta;
	m(1, 1) = cosTheta;

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
const float& Matrix2x2::operator()(const int row_index, const int col_index) const
{
	return m(row_index, col_index);
}

///////////////////////////////////////////////////////////////////////////////
float& Matrix2x2::operator()(const int row_index, const int col_index)
{
	return m(row_index, col_index);
}

///////////////////////////////////////////////////////////////////////////////
const float& Matrix2x2::operator()(const int index) const
{
	assert(index >= 0 && index < 4);

	return data_.elements_[index];
}

///////////////////////////////////////////////////////////////////////////////
float& Matrix2x2::operator()(const int index)
{
	assert(index >= 0 && index < 4);

	return data_.elements_[index];
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2& Matrix2x2::operator=(const Matrix2x2& other)
{
	m(0, 0) = other(0, 0);
	m(0, 1) = other(0, 1);
	m(1, 0) = other(1, 0);
	m(1, 1) = other(1, 1);

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2 Matrix2x2::operator+(const Matrix2x2& other) const
{
	Matrix2x2 result;

	result(0, 0) = m(0, 0) + other(0, 0);
	result(0, 1) = m(0, 1) + other(0, 1);
	result(1, 0) = m(1, 0) + other(1, 0);
	result(1, 1) = m(1, 1) + other(1, 1);

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2& Matrix2x2::operator+=(const Matrix2x2& other)
{
	*this = *this + other;

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2 Matrix2x2::operator-(const Matrix2x2& other) const
{
	Matrix2x2 result;

	result(0, 0) = m(0, 0) - other(0, 0);
	result(0, 1) = m(0, 1) - other(0, 1);
	result(1, 0) = m(1, 0) - other(1, 0);
	result(1, 1) = m(1, 1) - other(1, 1);

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2& Matrix2x2::operator-=(const Matrix2x2& other)
{
	*this = *this - other;

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2 Matrix2x2::operator*(const float scalar) const
{
	Matrix2x2 result;

	result(0, 0) = m(0, 0) * scalar;
	result(0, 1) = m(0, 1) * scalar;
	result(1, 0) = m(1, 0) * scalar;
	result(1, 1) = m(1, 1) * scalar;

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2& Matrix2x2::operator*=(const float scalar)
{
	*this = *this * scalar;

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2 Matrix2x2::operator*(const Matrix2x2& other) const
{
	return NaiveMatrixProduct(*this, other);
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2& Matrix2x2::operator*=(const Matrix2x2& other)
{
	*this = NaiveMatrixProduct(*this, other);

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Vector2d Matrix2x2::operator*(const Vector2d& vector)
{
	Vector2d result;

	result(0) = m(0, 0) * vector(0) + m(0, 1) * vector(1);
	result(1) = m(1, 0) * vector(0) + m(1, 1) * vector(1);

	return result;
}

///////////////////////////////////////////////////////////////////////////////
bool Matrix2x2::operator==(const Matrix2x2& other) const
{
	return IsApproximatelyEqual(other, 0.0f);
}

///////////////////////////////////////////////////////////////////////////////
bool Matrix2x2::operator!=(const Matrix2x2& other) const
{
	return !(*this == other);
}

///////////////////////////////////////////////////////////////////////////////
bool Matrix2x2::IsApproximatelyEqual(const Matrix2x2& other, const float epsilon) const
{
	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < 2; ++j)
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
Matrix2x2 Matrix2x2::NaiveMatrixProduct(const Matrix2x2& A, const Matrix2x2& B)
{
	Matrix2x2 result;

	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			for (int k = 0; k < 2; ++k)
			{
				result(i, j) += A(i, k) * B(k, j);
			}
		}
	}

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2 Matrix2x2::UnrolledMatrixProduct(const Matrix2x2& A, const Matrix2x2& B)
{
	Matrix2x2 result;

	// Same as above, but without the three for loops, i.e. with loops unrolled.
	// Notice that there are 8 multiplications in total.
	result(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0);
	result(0, 1) = A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1);
	result(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0);
	result(1, 1) = A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1);

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2 Matrix2x2::StrassenMatrixProduct(const Matrix2x2& A, const Matrix2x2& B)
{
	// Strassen's algorithm for 2x2 matrices. Note that we have decreased the number
	// of multiplications to 7, at the cost of more addition and subtraction operations.
	const float a_11 = A(0, 0);
	const float a_12 = A(0, 1);
	const float a_21 = A(1, 0);
	const float a_22 = A(1, 1);
	const float b_11 = B(0, 0);
	const float b_12 = B(0, 1);
	const float b_21 = B(1, 0);
	const float b_22 = B(1, 1);

	const float m_1 = (a_11 + a_22) * (b_11 + b_22);
	const float m_2 = (a_21 + a_22) * b_11;
	const float m_3 = a_11 * (b_12 - b_22);
	const float m_4 = a_22 * (b_21 - b_11);
	const float m_5 = (a_11 + a_12) * b_22;
	const float m_6 = (a_21 - a_11) * (b_11 + b_12);
	const float m_7 = (a_12 - a_22) * (b_21 + b_22);

	const float c_11 = (m_1 + m_4 - m_5 + m_7);
	const float c_12 = (m_3 + m_5);
	const float c_21 = (m_2 + m_4);;
	const float c_22 = (m_1 - m_2 + m_3 + m_6);

	const Matrix2x2 result(c_11, c_12, c_21, c_22);

	return result;
}

///////////////////////////////////////////////////////////////////////////////
float Matrix2x2::Trace() const
{
	const float trace = m(0, 0) + m(1, 1);

	return trace;
}

///////////////////////////////////////////////////////////////////////////////
float Matrix2x2::Determinant() const
{
	const float det = m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);

	return det;
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2& Matrix2x2::Transpose()
{
	std::swap(m(0, 1), m(1, 0));

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2 Matrix2x2::Transposed() const
{
	Matrix2x2 result(*this);

	result.Transpose();

	return result;
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2& Matrix2x2::Invert()
{
	const float det = Determinant();

	assert(det != 0);

	Matrix2x2 cofactor_matrix;

	cofactor_matrix(0, 0) = m(1, 1);
	cofactor_matrix(0, 1) = -m(1, 0);
	cofactor_matrix(1, 0) = -m(0, 1);
	cofactor_matrix(1, 1) = m(0, 0);

	*this = cofactor_matrix.Transpose() * (1 / det);

	return *this;
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2 Matrix2x2::Inverse() const
{
	Matrix2x2 result(*this);

	result.Invert();

	return result;
}

///////////////////////////////////////////////////////////////////////////////
std::string Matrix2x2::ToString() const
{
	const std::string result =
		"\n" + data_.rows_[0].ToString() +
		"\n" + data_.rows_[1].ToString();

	return result;
}

///////////////////////////////////////////////////////////////////////////////
const float& Matrix2x2::m(const int row_index, const int col_index) const
{
	assert(row_index >= 0 && row_index < 2);
	assert(col_index >= 0 && col_index < 2);

	return data_.array_[row_index][col_index];
}

///////////////////////////////////////////////////////////////////////////////
float& Matrix2x2::m(const int row_index, const int col_index)
{
	assert(row_index >= 0 && row_index < 2);
	assert(col_index >= 0 && col_index < 2);

	return data_.array_[row_index][col_index];
}

///////////////////////////////////////////////////////////////////////////////
Matrix2x2 operator*(const float scalar, const Matrix2x2& matrix)
{
return matrix * scalar;
}

} // namespace aleph0

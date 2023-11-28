#pragma once

// Standard C++ library headers
#include <cassert>
#include <cmath>
#include <string>
#include <utility>

// Aleph0 library headers
#include "matrix2x2.h"

namespace aleph0
{

class Matrix3x3
{
public:
	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3()
	{
		SetZeroMatrix();
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3(const Matrix3x3& other)
	{
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				m(i,j) = other(i,j);
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3(const float elements[9])
	{
		for (int i = 0; i < 9; ++i)
		{
			data_.elements_[i] = elements[i];
		}
	}

public:
	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3& SetZeroMatrix()
	{
		memset(data_.elements_, 0, 9 * sizeof(float));

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3& SetIdentityMatrix()
	{
		SetZeroMatrix();

		m(0,0) = 1.f;
		m(1,1) = 1.f;
		m(2,2) = 1.f;

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	const float& operator()(const int row_index, const int col_index) const
	{
		return m(row_index, col_index);
	}

	///////////////////////////////////////////////////////////////////////////////
	float& operator()(const int row_index, const int col_index)
	{
		return m(row_index, col_index);
	}

	///////////////////////////////////////////////////////////////////////////////
	const float& operator()(const int index) const
	{
		assert(index >= 0 && index < 9);

		return data_.elements_[index];
	}

	///////////////////////////////////////////////////////////////////////////////
	float& operator()(const int index)
	{
		assert(index >= 0 && index < 9);

		return data_.elements_[index];
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3 operator+(const Matrix3x3& other) const
	{
		Matrix3x3 result;

		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				result(i,j) = m(i,j) + other(i,j);
			}
		}

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3& operator+=(const Matrix3x3& other)
	{
		*this = *this + other;

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3 operator-(const Matrix3x3& other) const
	{
		Matrix3x3 result;

		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				result(i,j) = m(i,j) - other(i,j);
			}
		}

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3& operator-=(const Matrix3x3& other)
	{
		*this = *this - other;

		return *this;
	}


	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3 operator*(const float scalar) const
	{
		Matrix3x3 result;

		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				result(i,j) = m(i,j) *  scalar;
			}
		}

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3& operator*=(const float scalar)
	{
		*this = *this * scalar;

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3 operator*(const Matrix3x3& other) const
	{
		return NaiveMatrixProduct(*this, other);
	}

	///////////////////////////////////////////////////////////////////////////////
	bool operator==(const Matrix3x3& other) const
	{
		return IsApproximatelyEqual(other);
	}

	///////////////////////////////////////////////////////////////////////////////
	bool operator!=(const Matrix3x3& other) const
	{
		return !(*this == other);
	}

	///////////////////////////////////////////////////////////////////////////////
	bool IsApproximatelyEqual(const Matrix3x3& other, const float epsilon = 0.f) const
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
	static Matrix3x3 NaiveMatrixProduct(const Matrix3x3& A, const Matrix3x3& B)
	{
		Matrix3x3 result;

		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				for (int k = 0; k < 3; ++k)
				{
					result(i,j) += A(i,k) * B(k,j);
				}
			}
		}

		return result;
	}

	//////////////////////////////////////////////////////////////////////////////
	float Trace() const
	{
		const float trace = m(0,0) + m(1,1) + m(2,2) ;

		return trace;
	}

	///////////////////////////////////////////////////////////////////////////////
	float Determinant() const
	{
		const float det =
			m(0,0) * m(1,1) * m(2,2) +
			m(0,1) * m(1,2) * m(2,0) +
			m(0,2) * m(1,0) * m(2,1) -
			m(0,2) * m(1,1) * m(2,0) -
			m(0,0) * m(1,2) * m(2,1) -
			m(0,1) * m(1,0) * m(2,2);

		return det;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3& Transpose()
	{
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				if (i < j)
				{
					std::swap(m(i,j), m(j,i));
				}
			}
		}

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3 Transposed() const
	{
		Matrix3x3 result(*this);

		result.Transpose();

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3& Invert()
	{
		const float det = Determinant();

		assert(det != 0);

		Matrix3x3 cofactor_matrix;
			
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				cofactor_matrix(i,j) = std::powf(-1.f, (float)(i + j)) * Submatrix(i,j).Determinant();
			}
		}

		*this = cofactor_matrix.Transpose() * (1 / det);

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3 Inverse() const
	{
		Matrix3x3 result(*this);

		result.Invert();

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix2x2 Submatrix(int skip_row_index, int skip_col_index) const
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

				result(index++) = m(i,j);
			}
		}

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	std::string ToString() const
	{
		std::string result;

		for (int i = 0; i < 3; ++i)
		{
			result += "\n(";

			for (int j = 0; j < 3; ++j)
			{
				result += std::to_string(m(i,j));

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

private:
	///////////////////////////////////////////////////////////////////////////////
	const float& m(const int row_index, const int col_index) const
	{
		assert(row_index >= 0 && row_index < 3);
		assert(col_index >= 0 && col_index < 3);

		return data_.array_[row_index][col_index];
	}

	///////////////////////////////////////////////////////////////////////////////
	float& m(const int row_index, const int col_index)
	{
		assert(row_index >= 0 && row_index < 3);
		assert(col_index >= 0 && col_index < 3);

		return data_.array_[row_index][col_index];
	}

private:
	///////////////////////////////////////////////////////////////////////////////
	union Data
	{
		Data()
		{
			memset(elements_, 0, 9 * sizeof(float));
		}

		// We can represent our 3x3 matrix in a number of ways:

		// 1) An array of nine floats.
		float elements_[9];

		// 2) As a 3x3 array of floats.
		float array_[3][3];
	};

	Data data_;
}; // class Matrix3x3

///////////////////////////////////////////////////////////////////////////////
Matrix3x3 operator*(const float scalar, const Matrix3x3& matrix)
{
	return matrix * scalar;
}

} // namespace aleph0
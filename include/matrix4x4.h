#pragma once

// Aleph0 library headers
#include "aleph0.h"
#include "matrix2x2.h"
#include "matrix3x3.h"

namespace aleph0
{

class Matrix4x4
{
public:
	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4()
	{
		SetZeroMatrix();
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4(const Matrix4x4& other)
	{
		memcpy(data_.elements_, other.data_.elements_, 16 * sizeof(float));
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4(const float elements[16])
	{
		memcpy(data_.elements_, elements, 16 * sizeof(float));
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4(const Matrix2x2& M_11, const Matrix2x2& M_12, const Matrix2x2& M_21, const Matrix2x2& M_22)
	{
		m(0,0) = M_11(0,0);
		m(0,1) = M_11(0,1);
		m(1,0) = M_11(1,0);
		m(1,1) = M_11(1,1);

		m(0,2) = M_12(0,0);
		m(0,3) = M_12(0,1);
		m(1,2) = M_12(1,0);
		m(1,3) = M_12(1,1);

		m(2,0) = M_21(0,0);
		m(2,1) = M_21(0,1);
		m(3,0) = M_21(1,0);
		m(3,1) = M_21(1,1);

		m(2,2) = M_22(0,0);
		m(2,3) = M_22(0,1);
		m(3,2) = M_22(1,0);
		m(3,3) = M_22(1,1);
	}

public:
	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4& SetZeroMatrix()
	{
		memset(data_.elements_, 0, 16 * sizeof(float));

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4& SetIdentityMatrix()
	{
		SetZeroMatrix();

		m(0,0) = 1.f;
		m(1,1) = 1.f;
		m(2,2) = 1.f;
		m(3,3) = 1.f;

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
		assert(index >= 0 && index < 16);

		return data_.elements_[index];
	}

	///////////////////////////////////////////////////////////////////////////////
	float& operator()(const int index)
	{
		assert(index >= 0 && index < 16);

		return data_.elements_[index];
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4& operator=(const Matrix4x4& other)
	{
		memcpy(data_.elements_, other.data_.elements_, 16 * sizeof(float));

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4 operator+(const Matrix4x4& other) const
	{
		Matrix4x4 result;

		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				result(i,j) = m(i,j) + other(i,j);
			}
		}

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4& operator+=(const Matrix4x4& other)
	{
		*this = *this + other;

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4 operator-(const Matrix4x4& other) const
	{
		Matrix4x4 result;

		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				result(i,j) = m(i,j) - other(i,j);
			}
		}

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4& operator-=(const Matrix4x4& other)
	{
		*this = *this - other;

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4 operator*(const float scalar) const
	{
		Matrix4x4 result;

		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				result(i,j) = m(i,j) * scalar;
			}
		}

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4& operator*=(const float scalar)
	{
		*this = *this * scalar;

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4 operator*(const Matrix4x4& other) const
	{
		return NaiveMatrixProduct(*this, other);
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4& operator*=(const Matrix4x4& other)
	{
		*this = NaiveMatrixProduct(*this, other);

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	bool operator==(const Matrix4x4& other) const
	{
		return IsApproximatelyEqual(other, 0.0f);
	}

	///////////////////////////////////////////////////////////////////////////////
	bool operator!=(const Matrix4x4& other) const
	{
		return !(*this == other);
	}

	///////////////////////////////////////////////////////////////////////////////
	bool IsApproximatelyEqual(const Matrix4x4& other, const float epsilon = aleph0::epsilon) const
	{
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
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
	static Matrix4x4 NaiveMatrixProduct(const Matrix4x4& A, const Matrix4x4& B)
	{
		Matrix4x4 result;

		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				for (int k = 0; k < 4; ++k)
				{
					result(i,j) += A(i,k) * B(k,j);
				}
			}
		}

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	static Matrix4x4 StrassenMatrixProduct(const Matrix4x4& A, const Matrix4x4& B)
	{
		// Strassen's algorithm for 4x4 matrices. Compare it with the corresponding
		// implementation for 2x2 matrices.
		const Matrix2x2 A_11 = Matrix2x2(A(0,0), A(0,1), A(1,0), A(1,1));
		const Matrix2x2 A_12 = Matrix2x2(A(0,2), A(0,3), A(1,2), A(1,3));
		const Matrix2x2 A_21 = Matrix2x2(A(2,0), A(2,1), A(3,0), A(3,1));
		const Matrix2x2 A_22 = Matrix2x2(A(2,2), A(2,3), A(3,2), A(3,3));
		const Matrix2x2 B_11 = Matrix2x2(B(0,0), B(0,1), B(1,0), B(1,1));
		const Matrix2x2 B_12 = Matrix2x2(B(0,2), B(0,3), B(1,2), B(1,3));
		const Matrix2x2 B_21 = Matrix2x2(B(2,0), B(2,1), B(3,0), B(3,1));
		const Matrix2x2 B_22 = Matrix2x2(B(2,2), B(2,3), B(3,2), B(3,3));

		const Matrix2x2 M_1 = (A_11 + A_22) * (B_11 + B_22);
		const Matrix2x2 M_2 = (A_21 + A_22) * B_11;
		const Matrix2x2 M_3 = A_11 * (B_12 - B_22);
		const Matrix2x2 M_4 = A_22 * (B_21 - B_11);
		const Matrix2x2 M_5 = (A_11 + A_12) * B_22;
		const Matrix2x2 M_6 = (A_21 - A_11) * (B_11 + B_12);
		const Matrix2x2 M_7 = (A_12 - A_22) * (B_21 + B_22);

		const Matrix2x2 C_11 = (M_1 + M_4 - M_5 + M_7);
		const Matrix2x2 C_12 = (M_3 + M_5);
		const Matrix2x2 C_21 = (M_2 + M_4);;
		const Matrix2x2 C_22 = (M_1 - M_2 + M_3 + M_6);

		const Matrix4x4 result(C_11, C_12, C_21, C_22);

		return result;
	}

	//////////////////////////////////////////////////////////////////////////////
	float Trace() const
	{
		const float trace = m(0,0) + m(1,1) + m(2,2) + m(3,3);

		return trace;
	}

	///////////////////////////////////////////////////////////////////////////////
	float Determinant() const
	{
		const float a_1 = m(0,0) * m(1,1) - m(0,1) * m(1,0);
		const float a_2 = m(0,0) * m(1,2) - m(0,2) * m(1,0);
		const float a_3 = m(0,0) * m(1,3) - m(0,3) * m(1,0);
		const float a_4 = m(0,1) * m(1,2) - m(0,2) * m(1,1);
		const float a_5 = m(0,1) * m(1,3) - m(0,3) * m(1,1);
		const float a_6 = m(0,2) * m(1,3) - m(0,3) * m(1,2);
		const float b_1 = m(2,0) * m(3,1) - m(2,1) * m(3,0);
		const float b_2 = m(2,0) * m(3,2) - m(2,2) * m(3,0);
		const float b_3 = m(2,0) * m(3,3) - m(2,3) * m(3,0);
		const float b_4 = m(2,1) * m(3,2) - m(2,2) * m(3,1);
		const float b_5 = m(2,1) * m(3,3) - m(2,3) * m(3,1);
		const float b_6 = m(2,2) * m(3,3) - m(2,3) * m(3,2);
		const float det = a_1 * b_6 - a_2 * b_5 + a_3 * b_4 + a_4 * b_3 - a_5 * b_2 + a_6 * b_1;

		return det;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4& Transpose()
	{
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
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
	Matrix4x4 Transposed() const
	{
		Matrix4x4 result(*this);

		result.Transpose();

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4& Invert()
	{
		const float det = Determinant();

		assert(det != 0);

		Matrix4x4 cofactor_matrix;

		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				cofactor_matrix(i,j) = std::powf(-1.f, (float)(i + j)) * Submatrix(i,j).Determinant();
			}
		}

		*this = cofactor_matrix.Transpose() * (1 / det);

		return *this;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4 Inverse() const
	{
		Matrix4x4 result(*this);

		result.Invert();

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3 Submatrix(int skip_row_index, int skip_col_index) const
	{
		assert(skip_row_index >= 0 && skip_row_index < 4);
		assert(skip_col_index >= 0 && skip_col_index < 4);

		Matrix3x3 result;
		int index = 0;

		for (int i = 0; i < 4; ++i)
		{
			if (i == skip_row_index)
			{
				continue;
			}

			for (int j = 0; j < 4; ++j)
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
	std::string ToString() const
	{
		std::string result;

		for (int i = 0; i < 4; ++i)
		{
			result += "\n(";
				
			for (int j = 0; j < 4; ++j)
			{
				result += std::to_string(m(i,j));

				if (j + 1 < 4)
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
		assert(row_index >= 0 && row_index < 4);
		assert(col_index >= 0 && col_index < 4);

		return data_.array_[row_index][col_index];
	}

	///////////////////////////////////////////////////////////////////////////////
	float& m(const int row_index, const int col_index)
	{
		assert(row_index >= 0 && row_index < 4);
		assert(col_index >= 0 && col_index < 4);

		return data_.array_[row_index][col_index];
	}

private:
	///////////////////////////////////////////////////////////////////////////////
	union Data
	{
		Data()
		{
			memset(elements_, 0, 16 * sizeof(float));
		}

		// We can represent our 4x4 matrix in a number of ways:

		// 1) An array of sixteen floats.
		float elements_[16];

		// 2) As a 4x4 array of floats.
		float array_[4][4];
	};

	Data data_;
}; // class Matrix4x4

///////////////////////////////////////////////////////////////////////////////
Matrix4x4 operator*(const float scalar, const Matrix4x4& matrix)
{
	return matrix * scalar;
}

} // namespace aleph0

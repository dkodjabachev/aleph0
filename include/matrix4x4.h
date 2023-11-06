#pragma once

// Standard C++ library headers
#include <cassert>
#include <string>

namespace aleph0
{

class Matrix4x4
{
public:
	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4()
	{
		// The default constructor of the union Data will do the job
		// of initializing the matrix entries to zero.
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4(const Matrix4x4& other)
	{
		for (size_t i = 0; i < 4; ++i)
		{
			for (size_t j = 0; j < 4; ++j)
			{
				(*this)[i][j] = other[i][j];
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4(const Matrix2x2& M_11, const Matrix2x2& M_12, const Matrix2x2& M_21, const Matrix2x2& M_22)
	{
		data_.array_[0][0] = M_11[0][0];
		data_.array_[0][1] = M_11[0][1];
		data_.array_[1][0] = M_11[1][0];
		data_.array_[1][1] = M_11[1][1];

		data_.array_[0][2] = M_12[0][0];
		data_.array_[0][3] = M_12[0][1];
		data_.array_[1][2] = M_12[1][0];
		data_.array_[1][3] = M_12[1][1];

		data_.array_[2][0] = M_21[0][0];
		data_.array_[2][1] = M_21[0][1];
		data_.array_[3][0] = M_21[1][0];
		data_.array_[3][1] = M_21[1][1];

		data_.array_[2][2] = M_22[0][0];
		data_.array_[2][3] = M_22[0][1];
		data_.array_[3][2] = M_22[1][0];
		data_.array_[3][3] = M_22[1][1];
	}

public:
	///////////////////////////////////////////////////////////////////////////////
	void SetZeroMatrix()
	{
		memset(data_.elements_, 0, 16 * sizeof(float));
	}

	///////////////////////////////////////////////////////////////////////////////
	void SetIdentityMatrix()
	{
		SetZeroMatrix();

		(*this)[0][0] = 1.f;
		(*this)[1][1] = 1.f;
		(*this)[2][2] = 1.f;
		(*this)[3][3] = 1.f;
	}

	///////////////////////////////////////////////////////////////////////////////
	const float* operator[](const size_t row_index) const
	{
		assert(row_index < 4);

		return data_.array_[row_index];
	}

	///////////////////////////////////////////////////////////////////////////////
	float* operator[](const size_t row_index)
	{
		assert(row_index < 4);

		return data_.array_[row_index];
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4 operator+(const Matrix4x4& other) const
	{
		Matrix4x4 result;

		for (size_t i = 0; i < 4; ++i)
		{
			for (size_t j = 0; j < 4; ++j)
			{
				result[i][j] = (*this)[i][j] + other[i][j];
			}
		}

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4& operator+=(const Matrix4x4& other)
	{
		(*this) = (*this) + other;

		return (*this);
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4 operator-(const Matrix4x4& other) const
	{
		Matrix4x4 result;

		for (size_t i = 0; i < 4; ++i)
		{
			for (size_t j = 0; j < 4; ++j)
			{
				result[i][j] = (*this)[i][j] - other[i][j];
			}
		}

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4& operator-=(const Matrix4x4& other)
	{
		(*this) = (*this) - other;

		return (*this);
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4 operator*(const float scalar) const
	{
		Matrix4x4 result;

		for (size_t i = 0; i < 16; ++i)
		{
			result.data_.elements_[i] = data_.elements_[i] * scalar;
		}

		return result;
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4& operator*(const float scalar)
	{
		(*this) = (*this) * scalar;

		return (*this);
	}

	///////////////////////////////////////////////////////////////////////////////
	Matrix4x4 operator*(const Matrix4x4& other) const
	{
		return NaiveMatrixProduct(*this, other);
	}

	///////////////////////////////////////////////////////////////////////////////
	static Matrix4x4 NaiveMatrixProduct(const Matrix4x4& A, const Matrix4x4& B)
	{
		Matrix4x4 result;

		for (size_t i = 0; i < 4; ++i)
		{
			for (size_t j = 0; j < 4; ++j)
			{
				for (size_t k = 0; k < 4; ++k)
				{
					result[i][j] += A[i][k] * B[k][j];
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
		const Matrix2x2 A_11 = Matrix2x2(A[0][0], A[0][1], A[1][0], A[1][1]);
		const Matrix2x2 A_12 = Matrix2x2(A[0][2], A[0][3], A[1][2], A[1][3]);
		const Matrix2x2 A_21 = Matrix2x2(A[2][0], A[2][1], A[3][0], A[3][1]);
		const Matrix2x2 A_22 = Matrix2x2(A[2][2], A[2][3], A[3][2], A[3][3]);
		const Matrix2x2 B_11 = Matrix2x2(B[0][0], B[0][1], B[1][0], B[1][1]);
		const Matrix2x2 B_12 = Matrix2x2(B[0][2], B[0][3], B[1][2], B[1][3]);
		const Matrix2x2 B_21 = Matrix2x2(B[2][0], B[2][1], B[3][0], B[3][1]);
		const Matrix2x2 B_22 = Matrix2x2(B[2][2], B[2][3], B[3][2], B[3][3]);

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

	///////////////////////////////////////////////////////////////////////////////
	std::string ToString() const
	{
		std::string result;

		for (size_t i = 0; i < 4; ++i)
		{
			result += "\n(";
				
			for (size_t j = 0; j < 4; ++j)
			{
				result += std::to_string((*this)[i][j]);

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

} // namespace aleph0
#pragma once

#include "aleph0.h"
#include "matrix2x2.h"

namespace aleph0
{

class Vector3d;

class Matrix3x3
{
public:
	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3();

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3(const Matrix3x3& other);

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3(const float elements[9]);

public:
	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3& SetZeroMatrix();

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3& SetIdentityMatrix();

	///////////////////////////////////////////////////////////////////////////////
	const float& operator()(const int row_index, const int col_index) const;

	///////////////////////////////////////////////////////////////////////////////
	float& operator()(const int row_index, const int col_index);

	///////////////////////////////////////////////////////////////////////////////
	const float& operator()(const int index) const;

	///////////////////////////////////////////////////////////////////////////////
	float& operator()(const int index);

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3& operator=(const Matrix3x3& other);

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3 operator+(const Matrix3x3& other) const;

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3& operator+=(const Matrix3x3& other);

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3 operator-(const Matrix3x3& other) const;

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3& operator-=(const Matrix3x3& other);

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3 operator*(const float scalar) const;

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3& operator*=(const float scalar);

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3 operator*(const Matrix3x3& other) const;

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3& operator*=(const Matrix3x3& other);

	///////////////////////////////////////////////////////////////////////////////
	Vector3d operator*(const Vector3d& vector);

	///////////////////////////////////////////////////////////////////////////////
	bool operator==(const Matrix3x3& other) const;

	///////////////////////////////////////////////////////////////////////////////
	bool operator!=(const Matrix3x3& other) const;

	///////////////////////////////////////////////////////////////////////////////
	bool IsApproximatelyEqual(const Matrix3x3& other, const float epsilon = aleph0::epsilon) const;

	///////////////////////////////////////////////////////////////////////////////
	static Matrix3x3 NaiveMatrixProduct(const Matrix3x3& A, const Matrix3x3& B);

	//////////////////////////////////////////////////////////////////////////////
	float Trace() const;

	///////////////////////////////////////////////////////////////////////////////
	float Determinant() const;

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3& Transpose();

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3 Transposed() const;

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3& Invert();

	///////////////////////////////////////////////////////////////////////////////
	Matrix3x3 Inverse() const;

	///////////////////////////////////////////////////////////////////////////////
	void GetAxisAngle(Vector3d& axis, float& angle) const;

	///////////////////////////////////////////////////////////////////////////////
	Matrix2x2 Submatrix(int skip_row_index, int skip_col_index) const;

	///////////////////////////////////////////////////////////////////////////////
	std::string ToString() const;

private:
	///////////////////////////////////////////////////////////////////////////////
	const float& m(const int row_index, const int col_index) const;

	///////////////////////////////////////////////////////////////////////////////
	float& m(const int row_index, const int col_index);

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
Matrix3x3 operator*(const float scalar, const Matrix3x3& matrix);

} // namespace aleph0

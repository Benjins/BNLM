#ifndef BNLM_CORE_H_
#define BNLM_CORE_H_

#pragma once

#include "CppUtils/assert.h"
#include "CppUtils/macros.h"
#include "CppUtils/vector.h"
#include "CppUtils/bitset.h"

#include <math.h>

namespace BNLM {

template<typename _T>
struct MatrixDynamicBlock {
	// NOTE: Could have dangling pointers, this is just meant for temporaries
	_T* data = nullptr;
	int rows = 0;
	int cols = 0;

	_T operator()(int r, int c) const {
		ASSERT(r >= 0);
		ASSERT(c >= 0);
		ASSERT(r < rows);
		ASSERT(c < cols);
		return data[cols * r + c];
	}

	_T& operator()(int r, int c) {
		ASSERT(r >= 0);
		ASSERT(c >= 0);
		ASSERT(r < rows);
		ASSERT(c < cols);
		return data[cols * r + c];
	}
};

template<typename _T>
struct VectorDynamicBlock {
	// NOTE: Could have dangling pointers, this is just meant for temporaries
	_T* data = nullptr;
	int dims = 0;


	_T operator()(int d) const {
		ASSERT(d >= 0);
		ASSERT(d < dims);
		return data[d];
	}

	_T& operator()(int d) {
		ASSERT(d >= 0);
		ASSERT(d < dims);
		return data[d];
	}
};

// Bleeeegghh
template<typename _T, int _Rows, int _Cols>
struct Matrix;

template<typename _T, int _Dims>
struct Vector {
	static_assert(_Dims > 0, "Check vector dimensions");
	_T data[_Dims];

	_T operator()(int d) const {
		ASSERT(d >= 0);
		ASSERT(d < _Dims);
		return data[d];
	}

	_T& operator()(int d) {
		ASSERT(d >= 0);
		ASSERT(d < _Dims);
		return data[d];
	}

	operator VectorDynamicBlock<_T>() {
		VectorDynamicBlock<_T> blk;
		blk.data = data;
		blk.dims = _Dims;
		return blk;
	}

	operator VectorDynamicBlock<const _T>() const {
		VectorDynamicBlock<_T> blk;
		blk.data = data;
		blk.dims = _Dims;
		return blk;
	}

	Vector<_T, _Dims + 1> homo() const {
		Vector<_T, _Dims + 1> retVal;
		BNS_FOR_I(_Dims) {
			retVal.data[i] = data[i];
		}

		retVal.data[_Dims] = 1;
		return retVal;
	}

	static Vector<_T, _Dims> Zero() {
		Vector<_T, _Dims> zero;

		BNS_FOR_I(_Dims) {
			zero.data[i] = 0;
		}

		return zero;
	}

	Vector<_T, _Dims - 1> hnorm() const{
		Vector<_T, _Dims - 1> retVal;
		float s = 1.0f / data[_Dims - 1];
		BNS_FOR_I(_Dims - 1) {
			retVal.data[i] = data[i] * s;
		}
		return retVal;
	}

	template<int _NewDims>
	Vector<_T, _NewDims> subvec(int startIdx) const {
		static_assert(_NewDims <= _Dims, "Check dimensions on subvec");
		static_assert(_NewDims > 0, "Check dimensions on subvec");
		ASSERT(startIdx >= 0 && startIdx + _NewDims <= _Dims);
		Vector<_T, _NewDims> retVal;
		BNS_FOR_I(_NewDims) {
			retVal.data[i] = data[i + startIdx];
		}
		return retVal;
	}

	Vector<_T, _Dims> operator+(const Vector<_T, _Dims>& v) const {
		Vector<_T, _Dims> retVal;
		BNS_FOR_I(_Dims) {
			retVal(i) = data[i] + v(i);
		}
		return retVal;
	}

	void operator+=(const Vector<_T, _Dims>& v) {
		*this = *this + v;
	}

	void operator*=(const float s) {
		*this = *this * s;
	}

	void operator-=(const Vector<_T, _Dims>& v) {
		*this = *this - v;
	}

	void operator/=(const float s) {
		*this = *this / s;
	}

	Vector<_T, _Dims> operator*(const float s) const {
		Vector<_T, _Dims> retVal;
		BNS_FOR_I(_Dims) {
			retVal(i) = data[i] * s;
		}
		return retVal;
	}

	// TODO: Multiplay reciprocal? :p
	Vector<_T, _Dims> operator/(const float s) const {
		Vector<_T, _Dims> retVal;
		BNS_FOR_I(_Dims) {
			retVal(i) = data[i] / s;
		}
		return retVal;
	}

	Vector<_T, _Dims> operator-(const Vector<_T, _Dims>& v) const {
		Vector<_T, _Dims> retVal;
		BNS_FOR_I(_Dims) {
			retVal(i) = data[i] - v(i);
		}
		return retVal;
	}

	_T SquareMag() const {
		_T retVal = 0;
		BNS_FOR_I(_Dims) {
			retVal += BNS_SQR(data[i]);
		}
		return retVal;
	}

	_T Mag() const {
		_T sqrMag = SquareMag();
		return sqrt(sqrMag);
	}

	_T Normalize() {
		_T norm = Mag();

		BNS_FOR_I(_Dims) {
			data[i] /= norm;
		}

		return norm;
	}

	Vector<_T, _Dims> Normalized() const {
		return (*this) / Mag();
	}

	bool operator==(const Vector<_T, _Dims>& other) const {
		BNS_FOR_I(_Dims) {
			if (other.data[i] != data[i]) {
				return false;
			}
		}

		return true;
	}

	bool operator!=(const Vector<_T, _Dims>& other) const {
		BNS_FOR_I(_Dims) {
			if (other.data[i] == data[i]) {
				return false;
			}
		}

		return true;
	}
};

template<typename _T>
struct RemoveConst {
	typedef _T type;
};

template<typename _T>
struct RemoveConst<const _T> {
	typedef _T type;
};

template<typename _T, int _Rows, int _Cols>
struct MatrixBlockBase {
	_T* dataStart;
	// stride is in elements, not bytes
	int stride;

	typedef typename RemoveConst<_T>::type _MutableT;

	MatrixBlockBase(_T* _dataStart, int _stride) {
		dataStart = _dataStart;
		stride = _stride;
	}

	// access operators
	_T operator()(int r, int c) const {
		ASSERT(r >= 0);
		ASSERT(c >= 0);
		ASSERT(r < _Rows);
		ASSERT(c < _Cols);
		return dataStart[stride * r + c];
	}

	_T& operator()(int r, int c) {
		ASSERT(r >= 0);
		ASSERT(c >= 0);
		ASSERT(r < _Rows);
		ASSERT(c < _Cols);
		return dataStart[stride * r + c];
	}

	// multiply vec+matrix operator
	Vector<_MutableT, _Rows> operator*(const Vector<_MutableT, _Cols>& vec) const {
		Vector<_MutableT, _Rows> retVal;
		BNS_FOR_J(_Rows) {
			_MutableT val = 0;
			BNS_FOR_I(_Cols) {
				val += (vec.data[i] * dataStart[stride * j + i]);
			}
			retVal.data[j] = val;
		}

		return retVal;
	}

	// TODO: transpose?
	//MatrixBlock<_T, _Cols, _Rows> transpose() const {
	//	MatrixBlock<_T, _Cols, _Rows> t(dataStart, stride);
	//}

	/*
	void operator=(const Vector<_MutableT, _Rows>& vec) {
		static_assert(_Cols == 1, "JLDNG");
		BNS_FOR_I(_Rows) {
			dataStart[stride * i] = vec.data[i];
		}
	}

	void operator=(const Vector<_MutableT, _Cols>& vec) {
		static_assert(_Rows == 1, "JLDNG");
		BNS_FOR_I(_Cols) {
			dataStart[i] = vec.data[i];
		}
	}
	*/
};

template<typename _T, int _Rows, int _Cols>
struct MatrixBlock : MatrixBlockBase<_T, _Rows, _Cols> {
	MatrixBlock(_T* _dataStart, int _stride) : MatrixBlockBase<_T, _Rows, _Cols>(_dataStart, _stride) { }

	typedef typename RemoveConst<_T>::type _MutableT;

	void operator=(const Matrix<_MutableT, _Rows, _Cols>& mat) {
		BNS_FOR_J(_Rows) {
			BNS_FOR_I(_Cols) {
				this->dataStart[this->stride * j + i] = mat.data[j * _Cols + i];
			}
		}
	}
};


template<typename _T, int _Rows>
struct MatrixBlock<_T, _Rows, 1> : MatrixBlockBase<_T, _Rows, 1> {
	MatrixBlock(_T* _dataStart, int _stride) : MatrixBlockBase<_T, _Rows, 1>(_dataStart, _stride) { }

	typedef typename RemoveConst<_T>::type _MutableT;

	void operator=(const Vector<_MutableT, _Rows>& vec) {
		BNS_FOR_I(_Rows) {
			this->dataStart[this->stride * i] = vec.data[i];
		}
	}
};

template<typename _T, int _Cols>
struct MatrixBlock<_T, 1, _Cols> : MatrixBlockBase<_T, 1, _Cols> {
	MatrixBlock(_T* _dataStart, int _stride) : MatrixBlockBase<_T, 1, _Cols>(_dataStart, _stride) { }

	typedef typename RemoveConst<_T>::type _MutableT;

	void operator=(const Vector<_MutableT, _Cols>& vec) {
		BNS_FOR_I(_Cols) {
			this->dataStart[i] = vec.data[i];
		}
	}
};


// Stored row major (i.e. row0 | row1 | row2 ...)
template<typename _T, int _Rows, int _Cols>
struct Matrix {
	_T data[_Rows * _Cols];

	Matrix() { }

	Matrix(const MatrixBlock<_T, _Rows, _Cols>& blk) {
		BNS_FOR_J(_Cols) {
			BNS_FOR_I(_Rows) {
				(*this)(j, i) = blk(j, i);
			}
		}
	}

	Matrix(const MatrixBlock<const _T, _Rows, _Cols>& blk) {
		BNS_FOR_J(_Cols) {
			BNS_FOR_I(_Rows) {
				(*this)(j, i) = blk(j, i);
			}
		}
	}

	operator MatrixDynamicBlock<_T>() {
		MatrixDynamicBlock<_T> blk;
		blk.data = data;
		blk.rows = _Rows;
		blk.cols = _Cols;
		return blk;
	}

	operator MatrixDynamicBlock<const _T>() const {
		MatrixDynamicBlock<_T> blk;
		blk.data = data;
		blk.rows = _Rows;
		blk.cols = _Cols;
		return blk;
	}
	
	void ZeroOut() {
		BNS_FOR_I(_Rows * _Cols) {
			data[i] = 0;
		}
	}
	
	void LoadIdentity() {
		ZeroOut();
		BNS_FOR_I(BNS_MIN(_Rows, _Cols)) {
			data[i * _Cols + i] = 1.0f;
		}
	}
	
	_T operator()(int r, int c) const {
		ASSERT(r >= 0);
		ASSERT(c >= 0);
		ASSERT(r < _Rows);
		ASSERT(c < _Cols);
		return data[_Cols * r + c];
	}
	
	_T& operator()(int r, int c) {
		ASSERT(r >= 0);
		ASSERT(c >= 0);
		ASSERT(r < _Rows);
		ASSERT(c < _Cols);
		return data[_Cols * r + c];
	}
	
	// TODO: Template trickery?? :p
	Matrix<_T, _Cols, _Rows> transpose() const {
		Matrix<_T, _Cols, _Rows> t;
		BNS_FOR_J(_Rows) {
			BNS_FOR_I(_Cols) {
				t(i, j) = (*this)(j, i);
			}
		}
		return t;
	}

	Matrix<_T, _Rows, _Cols> operator*(const float s) const {
		Matrix<_T, _Rows, _Cols> retVal;
		BNS_FOR_I(_Rows * _Cols) {
			retVal.data[i] = data[i] * s;
		}
		return retVal;
	}

	Vector<_T, _Rows> operator*(const Vector<_T, _Cols>& vec) const {
		Vector<_T, _Rows> retVal;
		BNS_FOR_J(_Rows) {
			_T val = 0;
			BNS_FOR_I(_Cols) {
				val += (data[j * _Cols + i] * vec.data[i]);
			}
			retVal(j) = val;
		}

		return retVal;
	}

	template<int _OtherCols>
	Matrix<_T, _Rows, _OtherCols> operator*(const Matrix<_T, _Cols, _OtherCols>& other) const {
		Matrix<_T, _Rows, _OtherCols> retVal;
		
		BNS_FOR_I(_OtherCols) {
			BNS_FOR_J(_Rows) {
				_T val = 0;
				BNS_FOR_NAME(k, _Cols) {
					val += ((*this)(j, k) * other(k, i));
				}

				retVal(j, i) = val;
			}
		}
		
		return retVal;
	}

	template<int _NewR, int _NewC>
	MatrixBlock<_T, _NewR, _NewC> block(int rStart, int cStart) {
		MatrixBlock<_T, _NewR, _NewC> blk(&data[rStart * _Cols + cStart], _Cols);
		return blk;
	}

	template<int _NewR, int _NewC>
	MatrixBlock<const _T, _NewR, _NewC> block(int rStart, int cStart) const {
		MatrixBlock<const _T, _NewR, _NewC> blk(&data[rStart * _Cols + cStart], _Cols);
		return blk;
	}

	static Matrix<_T, _Rows, _Cols> Identity() {
		Matrix<_T, _Rows, _Cols> mat;
		mat.LoadIdentity();
		return mat;
	}
	
	static Matrix<_T, _Rows, _Cols> Zero() {
		Matrix<_T, _Rows, _Cols> mat;
		mat.ZeroOut();
		return mat;
	}
};

struct Vector2f : Vector<float, 2> {
	Vector2f() { }
	Vector2f(const Vector<float, 2>& orig) {
		data[0] = orig.data[0];
		data[1] = orig.data[1];
	}
	Vector2f(float _x, float _y) {
		data[0] = _x;
		data[1] = _y;
	}

	static Vector2f XAxis() { return Vector2f(1, 0); }
	static Vector2f YAxis() { return Vector2f(0, 1); }
	
	float x() const { return data[0]; }
	float& x() { return data[0]; }
	float y() const { return data[1]; }
	float& y() { return data[1]; }
};

struct Vector3f : Vector<float, 3> {
	Vector3f() { }
	Vector3f(const Vector<float, 3>& orig) {
		data[0] = orig.data[0];
		data[1] = orig.data[1];
		data[2] = orig.data[2];
	}
	Vector3f(float _x, float _y, float _z) {
		data[0] = _x;
		data[1] = _y;
		data[2] = _z;
	}

	static Vector3f XAxis() { return Vector3f(1, 0, 0); }
	static Vector3f YAxis() { return Vector3f(0, 1, 0); }
	static Vector3f ZAxis() { return Vector3f(0, 0, 1); }

	float x() const { return data[0]; }
	float& x() { return data[0]; }
	float y() const { return data[1]; }
	float& y() { return data[1]; }
	float z() const { return data[2]; }
	float& z() { return data[2]; }
};

typedef Matrix<float, 2, 2> Matrix2f;
typedef Matrix<float, 3, 3> Matrix3f;
typedef Matrix<float, 4, 4> Matrix4f;

// Quaternion
template<typename _T>
struct Quaternion {
	union {
		struct {
			_T w;
			_T x;
			_T y;
			_T z;
		};

		_T wxyz[4];
	};

	Quaternion() { }

	Quaternion(_T _w, _T _x, _T _y, _T _z) {
		w = _w;
		x = _x;
		y = _y;
		z = _z;
	}

	Quaternion(Vector3f vec) {
		w = 0.0f;
		x = vec.x();
		y = vec.y();
		z = vec.z();
	}

	Quaternion(Vector3f axis, _T angleDegrees) {
		_T angleRadians = angleDegrees * BNS_DEG2RAD;
		_T halfAngle = angleRadians / 2;
		w = cos(halfAngle);

		Vector3f normalizedAxis = axis.Normalized();
		_T sinHalfAngle = sin(halfAngle);

		y = normalizedAxis.y() * sinHalfAngle;
		z = normalizedAxis.z() * sinHalfAngle;
		x = normalizedAxis.x() * sinHalfAngle;
	}

	Quaternion<_T> Normalized() const {
		_T sqrMag = w * w + x * x + y * y + z * z;
		_T s = 1 / sqrt(sqrMag);

		Quaternion<_T> ret;
		ret.w = w * s;
		ret.x = x * s;
		ret.y = y * s;
		ret.z = z * s;
		return ret;
	}

	Quaternion<_T> operator*(const Quaternion<_T>& multQuat) const {
		return Quaternion(w*multQuat.w - x * multQuat.x - y * multQuat.y - z * multQuat.z,
			              x*multQuat.w + w * multQuat.x - z * multQuat.y + y * multQuat.z,
			              y*multQuat.w + z * multQuat.x + w * multQuat.y - x * multQuat.z,
			              z*multQuat.w - y * multQuat.x + x * multQuat.y + w * multQuat.z);
	}

	Quaternion<_T> Conjugate() const {
		return Quaternion<_T>(w, -x, -y, -z);
	}

	static Quaternion<_T> Identity() {
		return Quaternion(1, 0, 0, 0);
	}
};

typedef Quaternion<float> Quaternionf;

Vector3f RotateVectorByQuaternion(Vector3f vec, Quaternionf quat) {
	Quaternionf rotate = quat.Normalized();
	Quaternionf vectorQuat = Quaternionf(0.0f, vec.x(), vec.y(), vec.z());
	Quaternionf conj = rotate.Conjugate();

	Quaternionf finalResult = rotate * (vectorQuat * conj);

	return Vector3f(finalResult.x, finalResult.y, finalResult.z);
}

// TODO: Rename to Quaternion to RotationMatrix
Matrix3f QuaternionToAngleAxis(Quaternionf quat) {
	Matrix3f mat;
	mat.block<3, 1>(0, 0) = RotateVectorByQuaternion(Vector3f::XAxis(), quat);
	mat.block<3, 1>(0, 1) = RotateVectorByQuaternion(Vector3f::YAxis(), quat);
	mat.block<3, 1>(0, 2) = RotateVectorByQuaternion(Vector3f::ZAxis(), quat);

	return mat;
}

Matrix3f AxisAngle (Vector3f axis, float degrees) {

	Quaternionf quat = Quaternionf(axis, degrees);

	return QuaternionToAngleAxis(quat);
}

Matrix2f RotationMat(float degrees) {
	float theta = degrees * BNS_DEG2RAD;
	float cosTheta = cos(theta);
	float sinTheta = sin(theta);

	Matrix2f mat;
	mat(0, 0) = cosTheta;
	mat(0, 1) = -sinTheta;
	mat(1, 0) = sinTheta;
	mat(1, 1) = cosTheta;

	return mat;
}

template<typename _T, int _Dims>
_T DotProduct(const Vector<_T, _Dims>& a, const Vector<_T, _Dims>& b) {
	_T acc = 0;
	BNS_FOR_I(_Dims) {
		acc += a.data[i] * b.data[i];
	}
	
	return acc;
}

Vector3f CrossProduct(const Vector3f& a, const Vector3f& b) {
	return Vector3f((a.y()*b.z() - a.z()*b.y()),
				    (a.z()*b.x() - a.x()*b.z()),
				    (a.x()*b.y() - a.y()*b.x()));
}

// TODO
template<typename _T>
struct VectorX {
	_T* data = nullptr;
	int dims = 0;

	// TODO: Capacity
	// int capacity = 0;

	VectorX() { }

	VectorX(int d) {
		SetSize(d);
	}

	operator VectorDynamicBlock<_T>() {
		VectorDynamicBlock<_T> blk;
		blk.data = data;
		blk.dims = dims;
		return blk;
	}

	operator VectorDynamicBlock<const _T>() const {
		VectorDynamicBlock<_T> blk;
		blk.data = data;
		blk.dims = dims;
		return blk;
	}

	VectorX(const VectorX<_T>& orig) {
		SetSize(orig.dims);
		if (orig.data != nullptr) {
			BNS_MEMCPY(data, orig.data, sizeof(_T) * dims);
		}
	}



	~VectorX() {
		Destroy();
	}

	void operator=(const VectorX<_T>& orig) {
		// TODO:
		// Uhhhh...how is self-assignment normally handled?
		if (data != orig.data) {
			SetSize(orig.dims);
			if (orig.data != nullptr) {
				BNS_MEMCPY(data, orig.data, sizeof(_T) * dims);
			}
		}
	}

	static inline VectorX<_T> Zero(int d) {
		VectorX<_T> vec(d);
		vec.ZeroOut();
		return vec;
	}

	void inline ZeroOut() {
		BNS_FOR_I(dims) {
			data[i] = 0;
		}
	}

	void Destroy() {
		if (data != nullptr) {
			delete[] data;
			data = nullptr;
		}
	}

	void SetSize(int d) {
		ASSERT(d >= 0);
		if (dims != d) {
			Destroy();
			data = new _T[d];
			dims = d;
		}
	}

	_T operator()(int d) const {
		ASSERT(d >= 0);
		ASSERT(d < dims);
		return data[d];
	}

	_T& operator()(int d) {
		ASSERT(d >= 0);
		ASSERT(d < dims);
		return data[d];
	}
};

// TODO
template<typename _T>
struct MatrixX {
	_T* data = nullptr;
	int rows = 0;
	int cols = 0;

	// TODO: Capacity, etc.
	//int capacity = 0;

	MatrixX() { }

	MatrixX(int r, int c) {
		SetSize(r, c);
	}

	MatrixX(const MatrixX<_T>& orig) {
		SetSize(orig.rows, orig.cols);
		if (orig.data != nullptr) {
			BNS_MEMCPY(data, orig.data, sizeof(_T) * rows * cols);
		}
	}

	operator MatrixDynamicBlock<_T>() {
		MatrixDynamicBlock<_T> blk;
		blk.data = data;
		blk.rows = rows;
		blk.cols = cols;
		return blk;
	}

	operator MatrixDynamicBlock<const _T>() const {
		MatrixDynamicBlock<_T> blk;
		blk.data = data;
		blk.rows = rows;
		blk.cols = cols;
		return blk;
	}

	void operator=(const MatrixX<_T>& orig) {
		SetSize(orig.rows, orig.cols);
		if (orig.data != nullptr) {
			BNS_MEMCPY(data, orig.data, sizeof(_T) * rows * cols);
		}
	}

	static inline MatrixX Zero(int r, int c) {
		MatrixX mat(r, c);
		mat.ZeroOut();
		return mat;
	}

	static inline MatrixX Identity(int r, int c) {
		MatrixX mat(r, c);
		mat.LoadIdentity();
		return mat;
	}

	inline void ZeroOut() {
		BNS_FOR_I(rows*cols) {
			data[i] = 0;
		}
	}

	inline void LoadIdentity() {
		BNS_FOR_I(rows) {
			BNS_FOR_J(cols) {
				data[i * cols + j] = (i == j) ? 1 : 0;
			}
		}
	}

	MatrixX<_T> transpose() const {
		MatrixX T(cols, rows);
		BNS_FOR_I(rows) {
			BNS_FOR_J(cols) {
				T(j, i) = data[i * cols + j];
			}
		}
		return T;
	}

	void Destroy() {
		if (data != nullptr) {
			delete[] data;
			data = nullptr;
		}
	}

	void SetSize(int r, int c) {
		if (rows != r || cols != c) {
			Destroy();

			data = new _T[r*c];
			rows = r;
			cols = c;
		}
	}

	_T operator()(int r, int c) const {
		ASSERT(r >= 0);
		ASSERT(c >= 0);
		ASSERT(r < rows);
		ASSERT(c < cols);
		return data[cols * r + c];
	}

	_T& operator()(int r, int c) {
		ASSERT(r >= 0);
		ASSERT(c >= 0);
		ASSERT(r < rows);
		ASSERT(c < cols);
		return data[cols * r + c];
	}

	// multiply vec+matrix operator
	VectorX<_T> operator*(const VectorX<_T>& vec) {
		ASSERT(vec.dims == cols);
		VectorX<_T> retVal(rows);
		BNS_FOR_J(rows) {
			_T val = 0;
			BNS_FOR_I(cols) {
				val += (vec.data[i] * data[cols * j + i]);
			}
			retVal.data[j] = val;
		}

		return retVal;
	}

	MatrixX<_T> operator*(const MatrixX<_T>& other) const {
		MatrixX<_T> retVal(rows, other.cols);

		BNS_FOR_I(other.cols) {
			BNS_FOR_J(rows) {
				_T val = 0;
				BNS_FOR_NAME(k, cols) {
					val += ((*this)(j, k) * other(k, i));
				}

				retVal(j, i) = val;
			}
		}

		return retVal;
	}

	~MatrixX() {
		Destroy();
	}
};

typedef MatrixX<float> MatrixXf;
typedef VectorX<float> VectorXf;

// TODO: Fix this cholesky decomp, it's broken
// TODO: Cholesky decomp for dynamic matrices
/*
template<int _Dim>
void CholeskyDecomposition(const Matrix<float, _Dim, _Dim>* mat, Matrix<float, _Dim, _Dim>* outL) {
	outL->ZeroOut();
	BNS_FOR_NAME(k, _Dim) {
		for (int i = 0; i < k; i++) {
			float acc = 0.0f;
			BNS_FOR_J(i) {
				acc += (*outL)(i, j) * (*outL)(k, j);
			}

			(*outL)(k, i) = (*mat)(k, i) - acc;
		}

		float acc = 0.0f;
		BNS_FOR_J(k) {
			acc += BNS_SQR((*outL)(k, j));
		}

		(*outL)(k, k) = sqrt((*mat)(k, k) - acc);
	}
}
*/

// TODO: Fix SVD Decomp U output
// TODO: Dynamic SVD decomp
int EigenDecomp_maxind(int k, const MatrixDynamicBlock<float>& mat) {
	ASSERT(mat.cols == mat.rows);
	int m = k + 1;
	for (int i = k + 2; i < mat.cols; i++) {
		if (BNS_ABS(mat(k, i)) > BNS_ABS(mat(k, m))) {
			m = i;
		}
	}

	return m;
}

void EigenDecomp_update(int k, float t, VectorDynamicBlock<float> outVals, BitSet* changed, int* state) {
	float y = outVals(k);
	(outVals)(k) = y + t;
	if (changed->GetBit(k) && y == (outVals)(k)) {
		changed->SetBit(k, false);
		(*state)--;
	}
	else if (!changed->GetBit(k) && y != (outVals)(k)) {
		changed->SetBit(k, true);
		(*state)++;
	}
}

void EigenDecomp_rotate(MatrixDynamicBlock<float> mat, int k, int l, int i, int j, float c, float s) {
	float newKL = c * mat(k, l) - s * mat(i, j);
	float newIJ = s * mat(k, l) + c * mat(i, j);
	mat(k, l) = newKL;
	mat(i, j) = newIJ;
}

void EigenDecomposition_internal(VectorDynamicBlock<float> outVals, MatrixDynamicBlock<float> outVecs, MatrixDynamicBlock<float> matScratch, VectorDynamicBlock<int> ind);

template<int _Dim>
void EigenDecomposition(const Matrix<float, _Dim, _Dim>* matOrig, Vector<float, _Dim>* outVals, Matrix<float, _Dim, _Dim>* outVecs) {
	Matrix<float, _Dim, _Dim> matScratch = *matOrig;
	Vector<int, _Dim> ind;
	outVecs->LoadIdentity();
	EigenDecomposition_internal(*outVals, *outVecs, matScratch, ind);
}

void EigenDecomposition(const MatrixX<float>* matOrig, VectorX<float>* outVals, MatrixX<float>* outVecs) {
	ASSERT(matOrig->cols == matOrig->rows);
	MatrixX<float> matScratch = *matOrig;
	VectorX<int> ind(matScratch.cols);

	outVals->SetSize(matOrig->cols);
	outVecs->SetSize(matOrig->rows, matOrig->cols);
	outVecs->LoadIdentity();
	EigenDecomposition_internal(*outVals, *outVecs, matScratch, ind);
}

void EigenDecomposition_internal(VectorDynamicBlock<float> outVals, MatrixDynamicBlock<float> outVecs, MatrixDynamicBlock<float> matScratch, VectorDynamicBlock<int> ind) {
	ASSERT(matScratch.cols == matScratch.rows);
	
	const int _Dim = matScratch.cols;

	int k, l, m, state;
	float s, c, t, p, y, d, r;
	
	BitSet changed;
	changed.EnsureCapacity(_Dim);

	state = _Dim;
	BNS_FOR_NAME(k, _Dim) {
		ind(k) = EigenDecomp_maxind(k, matScratch);
		outVals(k) = matScratch(k, k);
		changed.SetBit(k, true);
	}

	while (state != 0) {
		m = 0;
		for (int k = 1; k < _Dim - 1; k++) {
			if (BNS_ABS(matScratch(k, ind(k))) > BNS_ABS(matScratch(m, ind(m)))) {
				m = k;
			}
		}
		k = m;
		l = ind(m);
		p = matScratch(k, l);
		if (BNS_ABS(p) <= 0.000000001f) {
			break;
		}

		y = (outVals(l) - outVals(k)) * 0.5f;
		d = BNS_ABS(y) + sqrt(p * p + y * y);

		r = sqrt(p * p + d * d);
		if (r != 0.0f) {
			c = d / r;
			s = p / r;
		}
		else {
			c = 1.0f;
			s = 0.0f;
		}

		if (d != 0.0f) {
			t = p * p / d;
		}
		else {
			t = 0.0f;
		}

		if (y < 0) {
			s = -s;
			t = -t;
		}

		matScratch(k, l) = 0.0f;
		EigenDecomp_update(k, -t, outVals, &changed, &state);
		EigenDecomp_update(l, t, outVals, &changed, &state);

		for (int i = 0; i < k; i++) {
			EigenDecomp_rotate(matScratch, i, k, i, l, c, s);
		}

		for (int i = k + 1; i < l; i++) {
			EigenDecomp_rotate(matScratch, k, i, i, l, c, s);
		}

		for (int i = l + 1; i < _Dim; i++) {
			EigenDecomp_rotate(matScratch, k, i, l, i, c, s);
		}

		BNS_FOR_I(_Dim) {
			float newKI = c * outVecs(i, k) - s * outVecs(i, l);
			float newLI = s * outVecs(i, k) + c * outVecs(i, l);

			outVecs(i, k) = newKI;
			outVecs(i, l) = newLI;
		}

		ind(k) = EigenDecomp_maxind(k, matScratch);
		ind(l) = EigenDecomp_maxind(l, matScratch);
	}

	// Put eigen values in descending order, from greatest to least
	// TODO: Absoulte value?
	for (k = 0; k < _Dim - 1; k++) {
		m = k;
		for (int l = k + 1; l < _Dim; l++) {
			if (outVals(l) > outVals(m)) {
				m = l;
			}
		}
		if (m != k) {
			float tmp = outVals(m);
			outVals(m) = outVals(k);
			outVals(k) = tmp;

			BNS_FOR_I(_Dim) {
				tmp = outVecs(i, m);
				outVecs(i, m) = outVecs(i, k);
				outVecs(i, k) = tmp;
			}
		}
	}
}

// TODO: Broken, don't know why plz don't use
template<int _R, int _C>
void SingularValueDecomposition(const Matrix<float, _R, _C>& mat, Matrix<float, _R, _R>* outU, Vector<float, _C>* outS, Matrix<float, _C, _C>* outV) {
	Matrix<float, _C, _R> trans = mat.transpose();

	{
		Matrix<float, _C, _C> matTransMat = trans * mat;
		EigenDecomposition(&matTransMat, outS, outV);
	}

	BNS_FOR_I(_C) {
		(*outS)(i) = sqrtf((*outS)(i));
	}

	{
		Matrix<float, _R, _R> matMatTrans = mat * trans;
		Vector<float, _R> leftSingularLol;
		EigenDecomposition(&matMatTrans, &leftSingularLol, outU);


	}
}

void SingularValueDecomposition(const MatrixX<float>& mat, MatrixX<float>* outU, VectorX<float>* outS, MatrixX<float>* outV) {
	MatrixX<float> trans = mat.transpose();

	{
		MatrixX<float> matTransMat = trans * mat;
		EigenDecomposition(&matTransMat, outS, outV);
	}

	BNS_FOR_I(outS->dims) {
		(*outS)(i) = sqrtf((*outS)(i));
	}

	{
		MatrixX<float> matMatTrans = mat * trans;
		VectorX<float> leftSingularLol;
		EigenDecomposition(&matMatTrans, &leftSingularLol, outU);


	}
}

}

#endif

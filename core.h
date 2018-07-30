#ifndef BNLM_CORE_H_
#define BNLM_CORE_H_

#pragma once

#include "../CppUtils/assert.h"
#include "../CppUtils/macros.h"
#include "../CppUtils/vector.h"

namespace BNLM {

// Bleeeegghh
template<typename _T, int _Rows, int _Cols>
struct Matrix;

template<typename _T, int _Dims>
struct Vector {
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

	Vector<_T, _Dims + 1> homo() const {
		Vector<_T, _Dims + 1> retVal;
		BNS_FOR_I(_Dims) {
			retVal.data[i] = data[i];
		}

		retVal.data[_Dims] = 1.0f;
		return retVal;
	}

	// TODO: Issue with _Dims = 1?
	Vector<_T, _Dims - 1> hnorm() const{
		Vector<_T, _Dims - 1> retVal;
		float s = 1.0f / data[_Dims - 1];
		BNS_FOR_I(_Dims - 1) {
			retVal.data[i] = data[i] * s;
		}
		return retVal;
	}

	template<int _NewDims>
	Vector<_T, _NewDims> subvec(int startIdx) {
		static_assert(_NewDims <= _Dims, "Check dimensions on subvec");
		static_assert(_NewDims > 0, "Check dimensions on subvec");
		ASSERT(startIdx >= 0 && startIdx + _NewDims <= _Dims);
		Vector<_T, _NewDims> retVal;
		BNS_FOR_I(_NewDims) {
			retVal.data[i] = data[i + startIdx];
		}
		return retVal;
	}

	Vector<_T, _Dims> operator+(const Vector<_T, _Dims>& v) {
		Vector<_T, _Dims> retVal;
		BNS_FOR_I(_Dims) {
			retVal(i) = data[i] + v(i);
		}
		return retVal;
	}

	Vector<_T, _Dims> operator*(const float s) {
		Vector<_T, _Dims> retVal;
		BNS_FOR_I(_Dims) {
			retVal(i) = data[i] * s;
		}
		return retVal;
	}

	// TODO: Multiplay reciprocal? :p
	Vector<_T, _Dims> operator/(const float s) {
		Vector<_T, _Dims> retVal;
		BNS_FOR_I(_Dims) {
			retVal(i) = data[i] / s;
		}
		return retVal;
	}

	Vector<_T, _Dims> operator-(const Vector<_T, _Dims>& v) {
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
struct MatrixBlock {
	_T* dataStart;
	// stride is in elements, not bytes
	int stride;

	typedef typename RemoveConst<_T>::type _MutableT;

	MatrixBlock(_T* _dataStart, int _stride) {
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
	Vector<_MutableT, _Rows> operator*(const Vector<_MutableT, _Cols>& vec) {
		Vector<_MutableT, _Rows> retVal;
		BNS_FOR_J(_Rows) {
			_MutableT val = 0;
			BNS_FOR_I(_Cols) {
				val += (vec.data[i] * dataStart[stride * j + i]);
			}
		}

		return retVal;
	}

	// TODO: transpose?
	//MatrixBlock<_T, _Cols, _Rows> transpose() const {
	//	MatrixBlock<_T, _Cols, _Rows> t(dataStart, stride);
	//}

	void operator=(const Vector<_MutableT, _Rows>& vec) {
		static_assert(_Cols == 1, "JLDNG");
		BNS_FOR_I(_Rows) {
			dataStart[stride * i] = vec.data[i];
		}
	}

	// assign?
	void operator=(const Matrix<_MutableT, _Rows, _Cols>& mat) {
		BNS_FOR_J(_Rows) {
			BNS_FOR_I(_Cols) {
				dataStart[stride * j + i] = mat.data[j * _Cols + i];
			}
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
	Matrix<_T, _Rows, _OtherCols> operator*(const Matrix<_T, _Cols, _OtherCols>& other) {
		Matrix<_T, _Rows, _OtherCols> retVal;
		
		BNS_FOR_I(_OtherCols) {
			BNS_FOR_J(_Rows) {
				_T val;
				BNS_FOR_NAME(k, _B) {
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
	Vector2f(const Vector<float, 2>& orig) {
		data[0] = orig.data[0];
		data[1] = orig.data[1];
	}
	Vector2f(float _x, float _y) {
		data[0] = _x;
		data[1] = _y;
	}
	
	float x() const { return data[0]; }
	float& x() { return data[0]; }
	float y() const { return data[1]; }
	float& y() { return data[1]; }
};

struct Vector3f : Vector<float, 3> {
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

//Matrix3f AngleAxis(float, Vector3f)

// Matrix2f RotationMat(float 

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

}

#endif

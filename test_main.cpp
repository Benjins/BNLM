
#include "core.h"


#include "CppUtils/assert.h"
#include "CppUtils/bitset.h"
#include "CppUtils/strings.h"
#include "CppUtils/macros.h"
#include "CppUtils/testing.h"

#include <time.h>
#include <stdlib.h>

#define TEST_EPS 0.001f

#define ASSERT_APPROX(x0, x1) do { ASSERT(BNS_ABS(x0 - x1) < TEST_EPS) } while(0)
#define ASSERT_APPROX_WITH_EPS(x0, x1, eps) do { ASSERT(BNS_ABS(x0 - x1) < eps) } while(0)
#define ASSERT_APPROX_V3(v0, v1) do { ASSERT_APPROX(v0.x(), v1.x()); ASSERT_APPROX(v0.y(), v1.y()); ASSERT_APPROX(v0.z(), v1.z()); } while(0)

CREATE_TEST_CASE("Basic BNLM") {
	{
		int currTime = time(NULL);
		printf("Seeding RNG with seed %d\n", currTime);
		srand(currTime);
	}

	{
		BNLM::Vector3f v1(1.0f, 2.0f, -3.0f);
		ASSERT(v1.x() ==  1.0f);
		ASSERT(v1.y() ==  2.0f);
		ASSERT(v1.z() == -3.0f);

		ASSERT(v1(0) == 1.0f);
		ASSERT(v1(1) == 2.0f);
		ASSERT(v1(2) == -3.0f);
	}

	{
		BNLM::Vector3f v1(1.0f, 1.0f, 1.0f);
		BNLM::Matrix<float, 3, 3> mat = BNLM::Matrix<float, 3, 3>::Identity();
		BNLM::Vector3f v2 = mat * v1;
		ASSERT(v1 == v2);
	}

	{
		BNLM::Vector3f v1(1.0f, 1.0f, 1.0f);
		BNLM::Matrix<float, 3, 3> mat = BNLM::Matrix<float, 3, 3>::Identity() * 4.0f;
		BNLM::Vector3f v2 = mat * v1;
		ASSERT(v1 * 4.0f == v2);
		ASSERT(v2.x() == 4.0f);
		ASSERT(v2.y() == 4.0f);
		ASSERT(v2.z() == 4.0f);
	}

	{
		BNLM::Vector3f v1(1.0f, 2.0f, 5.0f);
		BNLM::Matrix<float, 3, 3> mat = BNLM::Matrix<float, 3, 3>::Identity() * 3.0f;
		BNLM::Vector3f v2 = mat * v1;
		ASSERT(v2.x() ==  3.0f);
		ASSERT(v2.y() ==  6.0f);
		ASSERT(v2.z() == 15.0f);
	}

	{
		BNLM::Matrix<float, 3, 3> mat = BNLM::Matrix<float, 3, 3>::Identity();
		mat(0, 1) = -1.0f;
		mat(0, 2) =  2.0f;
		mat(1, 2) =  3.0f;
		BNLM::Matrix<float, 2, 2> mat2 = mat.block<2, 2>(0, 1);
		ASSERT(mat2(0, 0) == -1.0f);
		ASSERT(mat2(0, 1) ==  2.0f);
		ASSERT(mat2(1, 0) ==  1.0f);
		ASSERT(mat2(1, 1) ==  3.0f);
	}

	{
		BNLM::Matrix<float, 3, 3> mat = BNLM::Matrix<float, 3, 3>::Identity();
		mat(0, 1) = -1.0f;
		mat(0, 2) = 2.0f;
		mat(1, 2) = 3.0f;
		BNLM::Vector2f v1;
		v1(0) = 2.0f;
		v1(1) = 1.0f;

		BNLM::Vector2f v2 = mat.block<2, 2>(0, 1) * v1;

		ASSERT(v2(0) == 0.0f);
		ASSERT(v2(1) == 5.0f);
	}

	{
		BNLM::Vector3f v1(1.0f, 2.0f, 5.0f);
		BNLM::Vector2f v2 = v1.subvec<2>(0);
		ASSERT(v2.x() == 1.0f);
		ASSERT(v2.y() == 2.0f);

		BNLM::Vector2f v3 = v1.subvec<2>(1);
		ASSERT(v3.x() == 2.0f);
		ASSERT(v3.y() == 5.0f);
	}

	{
		BNLM::Matrix<float, 3, 3> mat = BNLM::Matrix<float, 3, 3>::Identity();
		mat(0, 1) = -1.0f;
		mat(0, 2) = 2.0f;
		mat(1, 2) = 3.0f;
		BNLM::Matrix<float, 2, 2> mat2 = mat.block<2, 2>(0, 1);
		ASSERT(mat2(0, 0) == -1.0f);
		ASSERT(mat2(0, 1) == 2.0f);
		ASSERT(mat2(1, 0) == 1.0f);
		ASSERT(mat2(1, 1) == 3.0f);

		mat.block<2, 2>(0, 0) = BNLM::Matrix2f::Identity();
		ASSERT(mat(0, 0) == 1.0f);
		ASSERT(mat(0, 1) == 0.0f);
		ASSERT(mat(1, 0) == 0.0f);
		ASSERT(mat(1, 1) == 1.0f);
	}

	{
		BNLM::Matrix3f mat = BNLM::Matrix3f::Identity();
		mat.block<2, 1>(0, 2) = BNLM::Vector2f(3.0f, 2.0f);
		ASSERT(mat(0, 2) == 3.0f);
		ASSERT(mat(1, 2) == 2.0f);
	}

	{
		BNLM::Matrix<float, 3, 3> mat1 = BNLM::Matrix<float, 3, 3>::Identity();
		BNLM::Matrix<float, 3, 3> mat2;



		BNS_FOR_NAME(iter, 1000) {
			BNS_FOR_I(3) {
				BNS_FOR_J(3) {
					float val = rand();
					// No divide by zero
					val = val / (float)(rand() + 1);
					mat2(i, j) = val;
				}
			}

			BNLM::Matrix<float, 3, 3> mat3 = mat1 * mat2;
			BNLM::Matrix<float, 3, 3> mat4 = mat2 * mat1;

			BNS_FOR_I(3) {
				BNS_FOR_J(3) {
					ASSERT(mat3(i, j) == mat2(i, j));
					ASSERT(mat4(i, j) == mat2(i, j));
				}
			}
		}

		{
			BNLM::Quaternionf quat = BNLM::Quaternionf::Identity();

			BNLM::Vector3f v0 = BNLM::Vector3f(1.5f, 2.5f, -0.5f);

			BNLM::Vector3f rotV0 = BNLM::RotateVectorByQuaternion(v0, quat);

			ASSERT_APPROX_V3(rotV0, v0);
		}

		{
			BNLM::Quaternionf quat = BNLM::Quaternionf(BNLM::Vector3f(0.0f, 0.0f, 1.0f), 90.0f);

			BNLM::Vector3f v0 = BNLM::Vector3f(1.5f, 2.5f, -0.5f);

			BNLM::Vector3f rotV0 = BNLM::RotateVectorByQuaternion(v0, quat);
			BNLM::Vector3f expRotV0 = BNLM::Vector3f(-2.5f, 1.5f, -0.5f);

			ASSERT_APPROX_V3(rotV0, expRotV0);
		}

		{
			BNLM::Matrix3f rot = BNLM::AxisAngle(BNLM::Vector3f(0.0f, 0.0f, 1.0f), 90.0f);

			BNLM::Vector3f v0 = BNLM::Vector3f(1.5f, 2.5f, -0.5f);

			BNLM::Vector3f rotV0 = rot * v0;
			BNLM::Vector3f expRotV0 = BNLM::Vector3f(-2.5f, 1.5f, -0.5f);

			ASSERT_APPROX_V3(rotV0, expRotV0);
		}

		{
			BNLM::Matrix2f rot = BNLM::Matrix2f::Identity();
			rot.block<2, 1>(0, 0) = BNLM::Vector2f(2.5f, 3.5f);
			ASSERT(rot(0, 0) == 2.5f);
			ASSERT(rot(1, 0) == 3.5f);
			rot.block<1, 2>(0, 0) = BNLM::Vector2f(5.5f, 6.5f);
			ASSERT(rot(0, 0) == 5.5f);
			ASSERT(rot(0, 1) == 6.5f);
		}
	}

	{
		BNLM::Vector<float, 5> vec = BNLM::Vector<float, 5>::Zero();
		BNS_FOR_I(5) {
			ASSERT(vec(i) == 0.0f);
		}
	}

	{
		BNLM::Vector3f maybeZ = BNLM::CrossProduct(BNLM::Vector3f::XAxis(), BNLM::Vector3f::YAxis());
		ASSERT_APPROX_V3(maybeZ, BNLM::Vector3f::ZAxis());
	}

	{
		BNLM::Vector3f maybeNegZ = BNLM::CrossProduct(BNLM::Vector3f::YAxis(), BNLM::Vector3f::XAxis());
		ASSERT_APPROX_V3(maybeNegZ, -BNLM::Vector3f::ZAxis());
	}

	{
		BNLM::Vector3f maybeZero = BNLM::CrossProduct(BNLM::Vector3f::XAxis(), BNLM::Vector3f::XAxis());
		BNLM::Vector3f zeroVec = BNLM::Vector3f::Zero();
		ASSERT_APPROX_V3(maybeZero, zeroVec);
	}

	{
		BNLM::VectorXf vec(5);
		BNS_FOR_I(5) {
			vec(i) = i + 0.5f;
		}

		BNS_FOR_I(5) {
			ASSERT(vec(i) == i + 0.5f);
		}
	}

	{
		BNLM::VectorXf vec(5);
		BNS_FOR_I(5) {
			vec(i) = i + 0.5f;
		}

		BNS_FOR_I(5) {
			ASSERT(vec(i) == i + 0.5f);
		}
	}

	{
		BNLM::MatrixXf mat(5, 2);
		BNS_FOR_I(5) {
			BNS_FOR_J(2) {
				mat(i, j) = i * 0.5f + j * 2.25f;
			}
		}

		BNS_FOR_I(5) {
			BNS_FOR_J(2) {
				ASSERT(mat(i, j) == i * 0.5f + j * 2.25f);
			}
		}
	}

	{
		BNLM::MatrixXf mat(5, 2);
		BNS_FOR_I(5) {
			BNS_FOR_J(2) {
				mat(i, j) = i * 0.5f + j * 2.25f;
			}
		}
	}

	{
		BNLM::MatrixXf mat(5, 2);
		mat.LoadIdentity();
		mat(4, 0) = 2.0f;
		mat(4, 1) = 1.0f;

		BNLM::VectorXf vec1(2);
		vec1(0) = 2.5f;
		vec1(1) = -4.5f;

		BNLM::VectorXf vec2 = mat * vec1;
		ASSERT(vec2.dims == 5);
		ASSERT(vec2(0) == 2.5f);
		ASSERT(vec2(1) == -4.5f);
		ASSERT(vec2(2) == 0.0f);
		ASSERT(vec2(3) == 0.0f);
		ASSERT(vec2(4) == 0.5f);
	}

	{
		BNLM::VectorXf vec1(5);
		vec1.ZeroOut();

		BNLM::VectorXf vec2(7);
		vec2 = vec1;
		ASSERT(vec2.dims == vec1.dims);
		ASSERT(vec2.dims == 5);
		vec1 = BNLM::VectorXf::Zero(17);
		ASSERT(vec1.dims == 17);
		BNS_FOR_I(17) {
			ASSERT(vec1(i) == 0);
		}
		vec2 = vec2;
		ASSERT(vec2.dims == 5);
	}

	{
		BNLM::MatrixXf mat1(5, 7);
		BNLM::MatrixXf mat2(43, 9);
		mat2 = mat1;
		ASSERT(mat2.rows == 5);
		ASSERT(mat2.cols == 7);
		mat1 = BNLM::MatrixXf::Identity(2, 7);
		ASSERT(mat1.rows == 2);
		ASSERT(mat1.cols == 7);
		ASSERT(mat1(0, 0) == 1.0f);
		ASSERT(mat1(0, 1) == 0.0f);
		ASSERT(mat1(1, 1) == 1.0f);
		mat2 = BNLM::MatrixXf::Zero(8, 19);
		ASSERT(mat2.rows == 8);
		ASSERT(mat2.cols == 19);
		BNS_FOR_I(8) {
			BNS_FOR_J(19) {
				ASSERT(mat2(i, j) == 0);
			}
		}
	}
#if 0
	{
		BNLM::Matrix3f Chol;
		Chol(0, 0) = 4;
		Chol(0, 1) = 12;
		Chol(0, 2) = -16;

		Chol(1, 0) = 12;
		Chol(1, 1) = 37;
		Chol(1, 2) = -43;

		Chol(0, 0) = -16;
		Chol(0, 1) = -43;
		Chol(0, 2) = 98;

		BNLM::Matrix3f Best;
		BNLM::CholeskyDecomposition(&Chol, &Best);
		int xc = 0;
		(void)xc;
	}

	{
		BNLM::Matrix3f B = BNLM::Matrix3f::Identity();
		B(0, 0) = 2.0f;
		B(1, 0) = 4.0f;
		B(2, 0) = 1.5f;
		B(2, 1) = -3.0f;

		BNLM::Matrix3f BT = B.transpose();
		BNLM::Matrix3f BBT = B * BT;

		BNLM::Matrix3f Best;
		BNLM::CholeskyDecomposition(&BBT, &Best);
		BNS_FOR_I(3) {
			BNS_FOR_J(3) {
				ASSERT_APPROX(B(i, j), Best(i, j));
			}
		}
	}
#endif

	{
		BNLM::Matrix3f mat = BNLM::Matrix3f::Identity();
		mat(0, 0) = 2.0f;
		mat(0, 1) = 3.0f;
		mat(0, 2) = -1.0f;

		mat(1, 0) = 5.0f;
		mat(1, 1) = -3.0f;
		mat(1, 2) = -2.0f;

		mat(2, 0) = 4.0f;
		mat(2, 1) = 1.0f;
		mat(2, 2) = 2.0f;

		BNLM::Matrix3f U, V;
		BNLM::Vector3f sigma;
		BNLM::SingularValueDecomposition(mat, &U, &sigma, &V);

		//printf("U matrix:\n");
		//printf("  %f %f %f\n", U(0, 0), U(0, 1), U(0, 2));
		//printf("  %f %f %f\n", U(1, 0), U(1, 1), U(1, 2));
		//printf("  %f %f %f\n\n", U(2, 0), U(2, 1), U(2, 2));
		//
		//printf("singular vals: (%f %f %f)\n\n", sigma(0), sigma(1), sigma(2));
		//
		//printf("V matrix:\n");
		//printf("  %f %f %f\n", V(0, 0), V(0, 1), V(0, 2));
		//printf("  %f %f %f\n", V(1, 0), V(1, 1), V(1, 2));
		//printf("  %f %f %f\n\n", V(2, 0), V(2, 1), V(2, 2));

		BNLM::Matrix3f diag = BNLM::Matrix3f::Identity();
		BNS_FOR_I(3) {
			diag(i, i) = sigma(i);
		}
		BNLM::Matrix3f reconstMat = U * diag * V.transpose();

		BNS_FOR_I(3) {
			BNS_FOR_J(3) {
				//ASSERT_APPROX(reconstMat(i, j), mat(i, j));
			}
		}
	}

	{
		BNLM::MatrixXf mat(3, 3);
		mat.LoadIdentity();
		mat(0, 0) = 2.0f;
		mat(0, 1) = 3.0f;
		mat(0, 2) = -1.0f;

		mat(1, 0) = 5.0f;
		mat(1, 1) = -3.0f;
		mat(1, 2) = -2.0f;

		mat(2, 0) = 4.0f;
		mat(2, 1) = 1.0f;
		mat(2, 2) = 2.0f;

		BNLM::MatrixXf U, V;
		BNLM::VectorXf sigma;
		BNLM::SingularValueDecomposition(mat, &U, &sigma, &V);

		ASSERT(sigma.dims == 3);
		ASSERT(U.rows == 3);
		ASSERT(U.cols == 3);
		ASSERT(V.rows == 3);
		ASSERT(V.cols == 3);

		BNLM::MatrixXf diag = BNLM::MatrixXf::Identity(3, 3);
		BNS_FOR_I(3) {
			diag(i, i) = sigma(i);
		}
		BNLM::MatrixXf reconstMat = U * diag * V.transpose();

		BNS_FOR_I(3) {
			BNS_FOR_J(3) {
				//ASSERT_APPROX(reconstMat(i, j), mat(i, j));
			}
		}
	}

	return 0;
}

CREATE_TEST_CASE("SVD minim") {
	BNLM::Matrix2f mat = BNLM::Matrix2f::Identity();
	mat(0, 0) = 2.0f;
	mat(0, 1) = 3.0f;

	mat(1, 0) = 0.0f;
	mat(1, 1) = -3.0f;

	BNLM::Matrix2f U, V;
	BNLM::Vector2f sigma;
	BNLM::SingularValueDecomposition(mat, &U, &sigma, &V);

	//BNS_FOR_I(2) {
	//	BNS_FOR_J(2) {
	//		U(i, j) *= -1.0f;
	//	}
	//}

	BNLM::Matrix2f diag = BNLM::Matrix2f::Identity();
	BNS_FOR_I(2) {
		diag(i, i) = sigma(i);
	}
	BNLM::Matrix2f reconstMat = U * diag * V.transpose();

	printf("mat: %f %f | %f %f\n", mat(0, 0), mat(0, 1), mat(1, 0), mat(1, 1));
	printf("rec: %f %f | %f %f\n", reconstMat(0, 0), reconstMat(0, 1), reconstMat(1, 0), reconstMat(1, 1));

	BNS_FOR_I(2) {
		BNS_FOR_J(2) {
			//ASSERT_APPROX(reconstMat(i, j), mat(i, j));
		}
	}

	return 0;
}


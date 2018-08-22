
#include "core.h"


#include "../CppUtils/assert.cpp"
#include "../CppUtils/macros.h"

#include <time.h>
#include <stdlib.h>

#define TEST_EPS 0.001f

#define ASSERT_APPROX(x0, x1) do { ASSERT(BNS_ABS(x0 - x1) < TEST_EPS) } while(0)
#define ASSERT_APPROX_V3(v0, v1) do { ASSERT_APPROX(v0.x(), v1.x()); ASSERT_APPROX(v0.y(), v1.y()); ASSERT_APPROX(v0.z(), v1.z()); } while(0)

int main() {
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

	return 0;
}


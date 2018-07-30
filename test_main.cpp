
#include "core.h"


#include "../CppUtils/assert.cpp"


int main() {

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

	return 0;
}


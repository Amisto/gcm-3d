#ifndef _GCM_ANISOTROPIC_MATRIX_3D_H
#define _GCM_ANISOTROPIC_MATRIX_3D_H  1

#include <assert.h>

#include "util/matrixes.h"
#include "Exception.h"

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

using namespace gcm;

namespace gcm {
	
	typedef struct
	{
		union 
		{
			float c[21];
			struct
			{
				float c11, c12, c13, c14, c15, c16;
				float      c22, c23, c24, c25, c26;
				float 		    c33, c34, c35, c36;
				float 				 c44, c45, c46;
				float 					  c55, c56;
				float 						   c66;
			};
		};
	} AnisotropicNumbers;
	
	class AnisotropicMatrix3D
	{
	public:
		AnisotropicMatrix3D();
		~AnisotropicMatrix3D();
		void prepare_matrix(const AnisotropicNumbers &numbers, float rho, int stage);
		float max_lambda();
		void DecompositeIt(gsl_matrix* a, gsl_matrix* u, gsl_matrix* l, gsl_matrix* u1);
		void RealChecker(gsl_vector_complex* a, gsl_matrix* l);
		void RealChecker(gsl_matrix_complex* a, gsl_matrix* u);
		
		gcm_matrix A;
		gcm_matrix L;
		gcm_matrix U;
		gcm_matrix U1;

	private:
		void CreateAx(const AnisotropicNumbers &numbers, float rho);
		void CreateAy(const AnisotropicNumbers &numbers, float rho);
		void CreateAz(const AnisotropicNumbers &numbers, float rho);
		void zero_all();
		void GcmToGsl(const gcm_matrix &a, gsl_matrix* b);
		void GslToGcm(gsl_matrix* a, gcm_matrix& b);
		void setZero(gcm_matrix& a);
		void setZero(gsl_matrix* a);
	};
}

#endif

#pragma once
#include <vector>
#include "math.h"

using namespace std;

inline double* cross( double* u,
	 double* v)
{

	double *vec = new double[3];
	vec[0] = *(u+1) *  *(v+2) - *(u+2) *  *(v+1);
	vec[1] = *(u+2) * *(v) - *(u) * *(v+2);
	vec[2] = *(u) * *(v+1) - *(u+1) * *(v);
	
	return &vec[0];
}


   template<typename T>
   inline T dot(T *uvec,  T *vvec)
   {
	   return *(uvec) * *(vvec) + *(uvec+1) * *(vvec+1) + *(uvec+2) * *(vvec+2);
   }


   template<typename T>
   inline T* dividecon(T* uvec, double &CONST)
   {
	
	   T *POM = new T[3];
	   POM[0] = *uvec / CONST;
	   POM[1] = *(uvec+1) / CONST;
	   POM[2] = *(uvec+2) / CONST;

	   return &POM[0];
   }

	template<typename T>
	inline void multconst(T* res, T* uvec, double &CONST)
	{
		res[0] = *uvec * CONST;
		res[1] = *(uvec + 1) * CONST;
		res[2] = *(uvec + 2) * CONST;
	}

   template<typename T>
   inline T norm2(T* vec)
   {
	   return dot<T>(vec, vec);
   }

   template<typename T>
   inline T norm(T* vec)
   {
	   return sqrt(norm2<T>(vec));
   }

	template<typename T>
	inline T norm2( T* uvec, T* vvec)
	{
		const T tmp0 = *uvec - *vvec;
		const T tmp1 = *(uvec+1) - *(vvec+1);
		const T tmp2 = *(uvec+2) - *(vvec+2);
		return tmp0 * tmp0 + tmp1 * tmp1 + tmp2 * tmp2;
	}

	 template<typename T>
	 inline  T norm(T *uvec, T *vvec)
	 {
		 return sqrt(norm2(uvec, vvec));
	 }


	  template<typename T>
	  inline void subtract(T* res, T* uvec, T* vvec)
	  {
		  res[0] = (*uvec) - (*vvec);
		  res[1] = *(uvec + 1) - *(vvec + 1);
		  res[2] = *(uvec + 2) - *(vvec + 2);
	  }

	  template<typename T>
	  inline  T* add(T* uvec, T* vvec)
	  {

		  T *zbroj = new T[3];

		  zbroj[0] = *uvec + *vvec;
		  zbroj[1] = *(uvec+1) + *(vvec+1);
		  zbroj[2] = *(uvec+2) + *(vvec+2);

		  return &zbroj[0];
	  }

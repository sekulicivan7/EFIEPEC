#pragma once
#include <vector>
#include "math.h"

using namespace std;

inline void cross(double* res, double* u,
	 double* v)
{

	res[0] = *(u+1) *  *(v+2) - *(u+2) *  *(v+1);
	res[1] = *(u+2) * *(v) - *(u) * *(v+2);
	res[2] = *(u) * *(v+1) - *(u+1) * *(v);
	
}


   template<typename T>
   inline T dot(T *uvec,  T *vvec)
   {
	   return *(uvec) * *(vvec) + *(uvec+1) * *(vvec+1) + *(uvec+2) * *(vvec+2);
   }


   template<typename T>
   inline void dividecon(T* res, T* uvec, double &CONST)
   {
	

	   res[0] = *uvec / CONST;
	   res[1] = *(uvec+1) / CONST;
	   res[2] = *(uvec+2) / CONST;


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
	  inline  void add(T* res, T* uvec, T* vvec)
	  {

		  res[0] = *uvec + *vvec;
		  res[1] = *(uvec+1) + *(vvec+1);
		  res[2] = *(uvec+2) + *(vvec+2);

	  }

#ifndef GCM_TETR_MESH_H_
#define GCM_TETR_MESH_H_

#include "../Mesh.h"
#include "../../elem/TetrFirstOrder.h"
#include "../../elem/TriangleFirstOrder.h"

using namespace gcm;

namespace gcm
{
	class TetrMesh: public Mesh {
		
	public:
		TetrMesh();
		virtual ~TetrMesh();
		
		virtual int getTetrsNumber() = 0;
		virtual int getTriangleNumber() = 0;
		
		virtual void addTetr(TetrFirstOrder* tetr) = 0;
		
		virtual TetrFirstOrder* getTetr(int index) = 0;
		
		virtual TetrFirstOrder* getTetrByLocalIndex(int index) = 0;
		
		virtual TriangleFirstOrder* getTriangle(int index) = 0;
		
		virtual void createTetrs(int number) = 0;
		
		virtual void createTriangles(int number) = 0;
		
		/*void preProcessGeometry();
		
		float getRecommendedTimeStep();

		float getMinH();
		
		void doNextPartStep(float tau, int stage);
		
		void checkTopology(float tau);*/
	};
}
#endif
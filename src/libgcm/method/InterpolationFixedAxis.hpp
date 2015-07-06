#ifndef GCM_INTERPOLATION_FIXED_AXIS_H
#define    GCM_INTERPOLATION_FIXED_AXIS_H

#include "libgcm/mesh/tetr/TetrMeshFirstOrder.hpp"
#include "libgcm/mesh/tetr/TetrMeshSecondOrder.hpp"
#include "libgcm/method/NumericalMethod.hpp"
#include "libgcm/mesh/Mesh.hpp"
//#include "libgcm/util/AnisotropicMatrix3D.hpp"
////#include "libgcm/util/ElasticMatrix3D.hpp"
#include "libgcm/util/Types.hpp"
#include "libgcm/node/CalcNode.hpp"
#include "libgcm/Logging.hpp"
#include "libgcm/Exception.hpp"


namespace gcm
{

    class InterpolationFixedAxis : public NumericalMethod {
    public:
        InterpolationFixedAxis();
        virtual ~InterpolationFixedAxis();
        int getNumberOfStages();
        void doNextPartStep(CalcNode& cur_node, CalcNode& new_node, float time_step, int stage, Mesh* mesh);
        std::string getType();
    protected:
        int prepare_node(CalcNode& cur_node, RheologyMatrixPtr rheologyMatrix,
                float time_step, int stage, Mesh* mesh,
                float* dksi, bool* inner, std::vector<CalcNode>& previous_nodes,
                float* outer_normal, bool debug);
        int prepare_node(CalcNode& cur_node, RheologyMatrixPtr rheologyMatrix,
                float time_step, int stage, Mesh* mesh,
                float* dksi, bool* inner, std::vector<CalcNode>& previous_nodes,
                float* outer_normal);
        int find_nodes_on_previous_time_layer(CalcNode& cur_node, int stage, Mesh* mesh,
                float dksi[], bool inner[], std::vector<CalcNode>& previous_nodes,
                float outer_normal[], bool debug);
        int find_nodes_on_previous_time_layer(CalcNode& cur_node, int stage, Mesh* mesh,
                float dksi[], bool inner[], std::vector<CalcNode>& previous_nodes,
                float outer_normal[]);
        
        void __doNextPartStep(CalcNode& cur_node, CalcNode& new_node, float time_step, int stage, Mesh* mesh);

        USE_LOGGER;
    private:
	void quasi_zerofy();
	void quasi_border(CalcNode& cur_node, CalcNode& new_node, CalcNode& virt_node, Mesh* vm, int dir, std::vector<CalcNode>& previous_nodes, bool* inner, float time_step);
	void quasi_contact(CalcNode& cur_node, CalcNode& new_node, CalcNode& virt_node, CalcNode& virt_node_inner, Mesh* vmi, std::vector<CalcNode>& previous_nodes, bool* inner, std::vector<CalcNode>& virt_previous_nodes, bool* virt_inner, float time_step, int dir);
	void quasi_volume(CalcNode& cur_node, CalcNode& new_node, CalcNode& virt_node_left, CalcNode& virt_node_right, Mesh* vml, Mesh* vmr, std::vector<CalcNode>& previous_nodes, bool* inner, float time_step);
    };
}

#endif    /* GCM_INTERPOLATION_FIXED_AXIS_H */


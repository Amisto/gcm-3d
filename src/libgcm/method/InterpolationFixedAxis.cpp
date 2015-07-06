/*
 * File:   InterpolationFixedAxis.cpp
 * Author: anganar
 *
 * Created on May 3, 2013, 12:00 AM
 */
#include "libgcm/method/InterpolationFixedAxis.hpp"

using namespace gcm;
using std::string;
using std::vector;

string InterpolationFixedAxis::getType()
{
    return "InterpolationFixedAxis";
}

InterpolationFixedAxis::InterpolationFixedAxis()
{
    INIT_LOGGER("gcm.method.InterpolationFixedAxis");
}

InterpolationFixedAxis::~InterpolationFixedAxis()
{
}

int InterpolationFixedAxis::getNumberOfStages()
{
    return 3;
}

void InterpolationFixedAxis::quasi_zerofy()
{
    //Don't do anything or zerofy the fuck of new_node? Let's try don't do anything.
}

void InterpolationFixedAxis::quasi_border(CalcNode& cur_node, CalcNode& new_node, CalcNode& virt_node, Mesh* vm, int dir, vector<CalcNode>& previous_nodes, bool* inner, float time_step)
{
    Engine &e = Engine::getInstance();
    RheologyMatrixPtr curM = cur_node.getRheologyMatrix();
    gcm::real distq = (cur_node.coords[0] - virt_node.coords[0])*(cur_node.coords[0] - virt_node.coords[0]) + (cur_node.coords[1] - virt_node.coords[1])*(cur_node.coords[1] - virt_node.coords[1]) + (cur_node.coords[2] - virt_node.coords[2])*(cur_node.coords[2] - virt_node.coords[2]);
    for (int i=0; i<9; i++)
	if (!inner[i] && curM->getL(i, i)*dir > EQUALITY_TOLERANCE)
	{
	    gcm::real lam_tau = curM->getL(i, i)*time_step;
	    if (lam_tau*lam_tau < distq)
	    {//segment interpolation
	    }
	    else
	    {//tetr interpolation
  	    }	
	    //FIXME Temporarily, hacked with virt values
	    for (int j=0; j<9; j++)
		previous_nodes[i].values[j] = virt_node.values[j];
	    inner[i] = true;
	}
    float outer_normal[3];
    outer_normal[abs(dir)-1] = dir/abs(dir);
    e.getBorderCondition(cur_node.getBorderConditionId())->doCalc(e.getCurrentTime(), cur_node,
                                                                new_node, cur_node.getRheologyMatrix(), previous_nodes, inner, outer_normal);
}

void InterpolationFixedAxis::quasi_contact(CalcNode& cur_node, CalcNode& new_node, CalcNode& virt_node, CalcNode& virt_node_inner, Mesh* vmi, vector<CalcNode>& previous_nodes, bool* inner, vector<CalcNode>& virt_previous_nodes, bool* virt_inner, float time_step, int dir)
{
    Engine &e = Engine::getInstance();
    RheologyMatrixPtr curM = cur_node.getRheologyMatrix();
    gcm::real distq = (cur_node.coords[0] - virt_node_inner.coords[0])*(cur_node.coords[0] - virt_node_inner.coords[0]) + (cur_node.coords[1] - virt_node_inner.coords[1])*(cur_node.coords[1] - virt_node_inner.coords[1]) + (cur_node.coords[2] - virt_node_inner.coords[2])*(cur_node.coords[2] - virt_node_inner.coords[2]);
    for (int i=0; i<9; i++)
	if (!inner[i] && curM->getL(i, i)*dir > EQUALITY_TOLERANCE)
	{
	    gcm::real lam_tau = curM->getL(i, i)*time_step;
            if (lam_tau*lam_tau < distq)
            {//segment interpolation
            }
            else
            {//tetr interpolation
            }
            //FIXME Temporarily, hacked with virt values
            for (int j=0; j<9; j++)
                previous_nodes[i].values[j] = virt_node_inner.values[j];
            inner[i] = true;
	}
    float outer_normal[3];
    outer_normal[abs(dir)-1] = dir/abs(dir);
    e.getContactCondition(cur_node.getContactConditionId())->doCalc(Engine::getInstance().getCurrentTime(), cur_node,
                                                       new_node, virt_node, cur_node.getRheologyMatrix(), previous_nodes, inner,
                                                       virt_node.getRheologyMatrix(), virt_previous_nodes, virt_inner, outer_normal);
}

void InterpolationFixedAxis::quasi_volume(CalcNode& cur_node, CalcNode& new_node, CalcNode& virt_node_left, CalcNode& virt_node_right, Mesh* vml, Mesh* vmr, vector<CalcNode>& previous_nodes, bool* inner, float time_step)
{
    Engine &e = Engine::getInstance();
    RheologyMatrixPtr curM = cur_node.getRheologyMatrix();
    gcm::real distql = (cur_node.coords[0] - virt_node_left.coords[0])*(cur_node.coords[0] - virt_node_left.coords[0]) + (cur_node.coords[1] - virt_node_left.coords[1])*(cur_node.coords[1] - virt_node_left.coords[1]) + (cur_node.coords[2] - virt_node_left.coords[2])*(cur_node.coords[2] - virt_node_left.coords[2]);
    gcm::real distqr = (cur_node.coords[0] - virt_node_right.coords[0])*(cur_node.coords[0] - virt_node_right.coords[0]) + (cur_node.coords[1] - virt_node_right.coords[1])*(cur_node.coords[1] - virt_node_right.coords[1]) + (cur_node.coords[2] - virt_node_right.coords[2])*(cur_node.coords[2] - virt_node_right.coords[2]);
    for (int i=0; i<9; i++)
	if (!inner[i])
	{
	    gcm::real lam_tau = curM->getL(i, i)*time_step;
	    if (curM->getL(i, i) > EQUALITY_TOLERANCE)
	    { //get values from left
		if (lam_tau*lam_tau < distql)
            	{//segment interpolation
            	}
            	else
            	{//tetr interpolation
            	}
		//FIXME Temporarily, hacked with virt values
		for (int j=0; j<9; j++)
		    previous_nodes[i].values[j] = virt_node_left.values[j];
		inner[i] = true;
	    }
	    else if (curM->getL(i, i) < -EQUALITY_TOLERANCE)
	    { //get valies from right
		if (lam_tau*lam_tau < distqr)
            	{//segment interpolation
            	}
            	else
            	{//tetr interpolation
            	}
                //FIXME Temporarily, hacked with virt values
                for (int j=0; j<9; j++)
                    previous_nodes[i].values[j] = virt_node_right.values[j];
		inner[i] = true;
	    }
	    else
	    {
		LOG_INFO("Zero lambda characteristic marked as outer o_O");
	 	return;
	    }
	}
    e.getVolumeCalculator("SimpleVolumeCalculator")->doCalc(cur_node, new_node, cur_node.getRheologyMatrix(), previous_nodes);
}

void InterpolationFixedAxis::__doNextPartStep(CalcNode& cur_node, CalcNode& new_node, float time_step, int stage, Mesh* mesh)
{
    assert_ge(stage, 0);
    assert_le(stage, 2);

    auto& engine = Engine::getInstance();

    LOG_TRACE("Start node prepare for node " << cur_node.number);
    LOG_TRACE("Node: " << cur_node);

    // Variables used in calculations internally

    // Delta x on previous time layer for all the omegas
    //     omega_new_time_layer(ksi) = omega_old_time_layer(ksi+dksi)
    float dksi[9];

    // If the corresponding point on previous time layer is inner or not
    bool inner[9];

    // We will store interpolated nodes on previous time layer here
    // We know that we need five nodes for each direction (corresponding to Lambdas -C1, -C2, 0, C2, C1)
    // TODO  - We can  deal with (lambda == 0) separately
    vector<CalcNode> previous_nodes;
    previous_nodes.resize(9);

    // Outer normal at current point
    float outer_normal[3];

    // Number of outer characteristics
    int outer_count = prepare_node(cur_node, cur_node.getRheologyMatrix(),
                                   time_step, stage, mesh,
                                   dksi, inner, previous_nodes,
                                   outer_normal);

    LOG_TRACE("Done node prepare");

    if(cur_node.number == -777) {
        LOG_INFO("Met: " << cur_node);
        LOG_INFO("Met: " << outer_count);
    }

    // If all the omegas are 'inner'
    // omega = Matrix_OMEGA * u
    // new_u = Matrix_OMEGA^(-1) * omega
    // TODO - to think - if all omegas are 'inner' can we skip matrix calculations and just use new_u = interpolated_u ?
    if (cur_node.isInner()) {
        LOG_TRACE("Start inner node calc");
        if (outer_count == 0)
            // FIXME - hardcoded name
            engine.getVolumeCalculator("SimpleVolumeCalculator")->doCalc(
                                                                          cur_node, new_node, cur_node.getRheologyMatrix(), previous_nodes);
        else
            THROW_BAD_MESH("Outer characteristic for internal node detected");
        LOG_TRACE("Done inner node calc");
    }

    if (cur_node.isBorder())
    {
        LOG_TRACE("Start border node calc");
        // FIXME_ASAP - do smth with this!
        // It is not stable now. See ugly hack below.
        // Think about: (a) cube, (b) rotated cube, (c) sphere.
        //int maxNormDirection = ( fabs(outer_normal[0]) > fabs(outer_normal[1])
        //            ? ( fabs(outer_normal[0]) > fabs(outer_normal[2] ? 0 : 2 ) )
        //            : ( fabs(outer_normal[1]) > fabs(outer_normal[2] ? 1 : 2 ) ) );
        //if( maxNormDirection != stage )
        //    return;
        float val = (outer_normal[stage] >= 0 ? 1.0 : -1.0);
        outer_normal[0] = outer_normal[1] = outer_normal[2] = 0;
        outer_normal[stage] = val;
        // If there is no 'outer' omega - it is ok, border node can be inner for some directions
        if (outer_count == 0)
        {
            // FIXME - hardcoded name
            engine.getVolumeCalculator("SimpleVolumeCalculator")->doCalc(
                                                                          cur_node, new_node, cur_node.getRheologyMatrix(), previous_nodes);
            return;
        }
	//LOG_INFO("01");
        RheologyMatrixPtr curM = cur_node.getRheologyMatrix();
        if (outer_count == 1)
        {   
            float valL = 0;
	    outer_count = 0;
            for (uint i=0; i<9; i++)
            {
                valL = (curM->getL(i, i) > 0 ? 1.0 : (fabs(curM->getL(i, i)) < EQUALITY_TOLERANCE ? 0 : -1.0));
                if (valL == 0) continue;
                if (val != valL) inner[i] = false;
		if (!inner[i]) outer_count++;
            }
        }
	//LOG_INFO("02");
	if (outer_count == 2 || outer_count == 4 || outer_count == 6) //double contact, we have another bodies on both sides
	{
	    //firstly, get virtual nodes and check materials 
            CalcNode virt_node_left, virt_node_right;
	    Mesh *ml, *mr;
	    //RheologyMatrixPtr lM, rM;
	    //LOG_INFO("L");
	    engine.getVirtNode(cur_node, virt_node_left, -stage-1, ml);
	    //LOG_INFO("R");
            engine.getVirtNode(cur_node, virt_node_right, stage+1, mr);
	    //LOG_INFO("D");
            if (virt_node_left.number) virt_node_left.setInContact(true);
	    if (virt_node_right.number) virt_node_right.setInContact(true);
            float virt_left_dksi[9], virt_right_dksi[9];
            bool virt_left_inner[9], virt_right_inner[9];
            vector<CalcNode> virt_left_previous_nodes, virt_right_previous_nodes;
            virt_left_previous_nodes.resize(9);
            virt_right_previous_nodes.resize(9);
            float virt_left_outer_normal[3], virt_right_outer_normal[3];
            int virt_left_outer_count;
	    if (virt_node_left.number) virt_left_outer_count = prepare_node(virt_node_left, virt_node_left.getRheologyMatrix(),
                                                    time_step, stage, ml,
                                                    virt_left_dksi, virt_left_inner, virt_left_previous_nodes,
                                                    virt_left_outer_normal);
            int virt_right_outer_count;
	    if (virt_node_right.number) virt_right_outer_count = prepare_node(virt_node_right, virt_node_right.getRheologyMatrix(),
                                                    time_step, stage, mr,
                                                    virt_right_dksi, virt_right_inner, virt_right_previous_nodes,
                                                    virt_right_outer_normal);
	    //Check Direction for virt nodes
	    RheologyMatrixPtr vlM, vrM;
	    float valL = 0;
	    if (virt_node_left.number)
	    {
		vlM = virt_node_left.getRheologyMatrix();
            	virt_left_outer_count = 0;
            	for (uint i=0; i<9; i++)
            	{
            	    valL = (vlM->getL(i, i) > 0 ? -1.0 : (fabs(vlM->getL(i, i)) < EQUALITY_TOLERANCE ? 0 : 1.0));
            	    if (valL == 0) continue;
            	    if (valL == -1) virt_left_inner[i] = false;
            	    if (!virt_left_inner[i]) virt_left_outer_count++;
            	}
	    }
	    if (virt_node_right.number)
	    {
            	vrM = virt_node_right.getRheologyMatrix();
            	virt_right_outer_count = 0;
            	for (uint i=0; i<9; i++)
            	{
                    valL = (vrM->getL(i, i) > 0 ? -1.0 : (fabs(vrM->getL(i, i)) < EQUALITY_TOLERANCE ? 0 : 1.0));
            	    if (valL == 0) continue;
            	    if (valL == 1) virt_right_inner[i] = false;
            	    if (!virt_right_inner[i]) virt_right_outer_count++;
            	}
	    }

	    //LOG_INFO("03");
	    if (!virt_node_left.number && !virt_node_right.number)
	    {
		quasi_zerofy(); //TODO
		return;
	    }
	    if (virt_node_left.number && !virt_node_right.number)	  
	    {
		cur_node.setMaterialId(virt_node_left.getMaterialId());
                //val = 1;
                //outer_normal[stage] = 1;
		if (virt_left_outer_count == 3 || virt_left_outer_count == 1) //quasi_border();
		    quasi_border(cur_node, new_node, virt_node_left, ml, stage+1, previous_nodes, inner, time_step);
		else quasi_zerofy(); //TODO
		return;
	    }
            if (!virt_node_left.number && virt_node_right.number)
            {
                cur_node.setMaterialId(virt_node_right.getMaterialId());
                //val = -1;
                //outer_normal[stage] = -1;
                if (virt_right_outer_count == 3 || virt_right_outer_count == 1) //quasi_border();
                    quasi_border(cur_node, new_node, virt_node_right, mr, -stage-1, previous_nodes, inner, time_step);
                else quasi_zerofy(); //TODO
		return;
            }
	    //LOG_INFO("04");
	    if (virt_node_left.getMaterialId() == virt_node_right.getMaterialId())
	    {
		if (virt_node_left.getMaterialId() != cur_node.getMaterialId())
		    cur_node.setMaterialId(virt_node_left.getMaterialId());
  	        if ((virt_left_outer_count == 3 || virt_left_outer_count == 1) && (virt_right_outer_count == 3 || virt_right_outer_count == 1)) 
		{
		    quasi_volume(cur_node, new_node, virt_node_left, virt_node_right, ml, mr, previous_nodes, inner, time_step); 
		    return;	
		}
		if (virt_left_outer_count == 3 || virt_left_outer_count == 1) 
		{
		    //val = 1;   
                    //outer_normal[stage] = 1;
		    //quasi_border();//
                    quasi_border(cur_node, new_node, virt_node_left, ml, stage+1, previous_nodes, inner, time_step);
		    return;	
		}
		if (virt_right_outer_count == 3 || virt_right_outer_count == 1) 
		{
		    //val = -1;   
                    //outer_normal[stage] = -1;
		    //quasi_border();
		    quasi_border(cur_node, new_node, virt_node_right, mr, -stage-1, previous_nodes, inner, time_step);
		    return;
		}
                quasi_zerofy(); //TODO
		return;
	    }

	    if (virt_node_left.getMaterialId() != cur_node.getMaterialId() && virt_node_right.getMaterialId() != cur_node.getMaterialId())
		cur_node.setMaterialId(virt_node_left.getMaterialId()); 
	    //"contact calculator"
	    if (virt_node_left.getMaterialId() == cur_node.getMaterialId())
	    {	
                if ((virt_left_outer_count == 3 || virt_left_outer_count == 1) && (virt_right_outer_count == 3 || virt_right_outer_count == 1))
                {  
                    //val = 1;
                    //outer_normal[stage] = 1;	
                    quasi_contact(cur_node, new_node, virt_node_right, virt_node_left, ml, previous_nodes, inner, virt_right_previous_nodes, virt_right_inner, time_step, stage + 1); 
                    return;
                }
                if (virt_left_outer_count == 3 || virt_left_outer_count == 1)
                {  
                    //val = 1;
                    //outer_normal[stage] = 1;
                    //quasi_border();
		    quasi_border(cur_node, new_node, virt_node_left, ml, stage+1, previous_nodes, inner, time_step);
                    return;
                }
                if (virt_right_outer_count == 3 || virt_right_outer_count == 1)
                {  
		    cur_node.setMaterialId(virt_node_right.getMaterialId());	  
                    //val = -1;
                    //outer_normal[stage] = -1;
                    //quasi_border();
		    quasi_border(cur_node, new_node, virt_node_right, mr, -stage-1, previous_nodes, inner, time_step);
                    return;
                }
                quasi_zerofy(); //TODO
		return;
	    }
  	    if (virt_node_right.getMaterialId() == cur_node.getMaterialId())
            {
                if ((virt_left_outer_count == 3 || virt_left_outer_count == 1) && (virt_right_outer_count == 3 || virt_right_outer_count == 1))
                {  
                    //val = -1;
                    //outer_normal[stage] = -1;
                    quasi_contact(cur_node, new_node, virt_node_left, virt_node_right, mr, previous_nodes, inner, virt_left_previous_nodes, virt_left_inner, time_step, -stage-1); 
                    return;
                }
                if (virt_left_outer_count == 3 || virt_left_outer_count == 1)
                {  
		    cur_node.setMaterialId(virt_node_left.getMaterialId());
                    //val = 1;
                    //outer_normal[stage] = 1;
                    //quasi_border();
		    quasi_border(cur_node, new_node, virt_node_left, ml, stage+1, previous_nodes, inner, time_step);
                    return;
                }
                if (virt_right_outer_count == 3 || virt_right_outer_count == 1)
                {  
                    //val = -1;
                    //outer_normal[stage] = -1;
                    //quasi_border();
		    quasi_border(cur_node, new_node, virt_node_right, mr, -stage-1, previous_nodes, inner, time_step);
                    return;
                }
                quasi_zerofy(); //TODO
	 	return;
            }
	    LOG_INFO("OBVIOUSLY, we forgot something 0 " <<outer_count);
	    return;
	}
	//LOG_INFO("06");
	if (outer_count == 3)
	{
	    //LOG_INFO("07");
	    CalcNode virt_node;
            Mesh *m;
            //RheologyMatrixPtr lM, rM;
	    //LOG_INFO("S");
            engine.getVirtNode(cur_node, virt_node, val*(stage+1), m); 
	    //LOG_INFO("D");
            if (virt_node.number) virt_node.setInContact(true);
            float virt_dksi[9];
            bool virt_inner[9];
            vector<CalcNode> virt_previous_nodes;
            virt_previous_nodes.resize(9);
            float virt_outer_normal[3];
            int virt_outer_count;
            if (virt_node.number)
	    { 
	  	//LOG_INFO("08");
		virt_outer_count = prepare_node(virt_node, virt_node.getRheologyMatrix(),
                                                    time_step, stage, m,
                                                    virt_dksi, virt_inner, virt_previous_nodes,
                                                    virt_outer_normal);
            	RheologyMatrixPtr virtM = virt_node.getRheologyMatrix();
            	float valL = 0;
		virt_outer_count = 0;
            	for (uint i=0; i<9; i++)
            	{
                    valL = (virtM->getL(i, i) > 0 ? -1.0 : (fabs(virtM->getL(i, i)) < EQUALITY_TOLERANCE ? 0 : 1.0));
                    if (valL == 0) continue;
                    if (val != valL) virt_inner[i] = false;
		    if (!virt_inner[i]) virt_outer_count++;
            	}
		if (virt_outer_count < 3) LOG_INFO("OBVIOUSLY, we fucked up smt 1 " <<virt_outer_count);
		if (virt_outer_count > 3) //quasi_border(); //FIXME - cut an unlikely branch of algorhythm
		    quasi_border(cur_node, new_node, virt_node, m, (stage+1)*val, previous_nodes, inner, time_step);
		//LOG_INFO("09");
		if (virt_outer_count == 3) engine.getContactCondition(cur_node.getContactConditionId())->doCalc(Engine::getInstance().getCurrentTime(), cur_node,
                                                       new_node, virt_node, cur_node.getRheologyMatrix(), previous_nodes, inner,
                                                       virt_node.getRheologyMatrix(), virt_previous_nodes, virt_inner, outer_normal);  
	 	//LOG_INFO("10");
		return;	
	    }
	    //LOG_INFO("11");
	    engine.getBorderCondition(cur_node.getBorderConditionId())->doCalc(Engine::getInstance().getCurrentTime(), cur_node,
                                                                 new_node, cur_node.getRheologyMatrix(), previous_nodes, inner, outer_normal);
	    //LOG_INFO("12");
	    return;
	}

	LOG_INFO("OBVIOUSLY, we forgot something 3 " <<outer_count);

	return;
        if (outer_count == 3)
        {
            // Border
            if (!cur_node.isInContact() || cur_node.contactDirection != stage) {
                // FIXME
                int borderCondId = cur_node.getBorderConditionId();
                //if(engine.getBorderCondition(borderCondId)->calc->getType() != "FreeBorderCalculator")
                //    LOG_INFO("Node: " << cur_node.number << ". Using calculator: " << engine.getBorderCondition(borderCondId)->calc->getType());
                engine.getBorderCondition(borderCondId)->doCalc(Engine::getInstance().getCurrentTime(), cur_node,
                                                                 new_node, cur_node.getRheologyMatrix(), previous_nodes, inner, outer_normal);
            }
            // Contact
            else
            {
                CalcNode& virt_node = engine.getVirtNode(cur_node.contactNodeNum);

                // FIXME - WA
                Mesh* virtMesh = (Mesh*) engine.getBody(virt_node.contactNodeNum)->getMeshes();
                // Mark virt node as having contact state
                // TODO FIXME - most probably CollisionDetector should do it
                // But we should check it anycase
                virt_node.setInContact(true);
                //virt_node.contactNodeNum = cur_node.contactNodeNum;

                // Variables used in calculations internally

                // Delta x on previous time layer for all the omegas
                //     omega_new_time_layer(ksi) = omega_old_time_layer(ksi+dksi)
                float virt_dksi[9];

                // If the corresponding point on previous time layer is inner or not
                bool virt_inner[9];

                // We will store interpolated nodes on previous time layer here
                // We know that we need five nodes for each direction (corresponding to Lambdas -C1, -C2, 0, C2, C1)
                // TODO  - We can  deal with (lambda == 0) separately
                vector<CalcNode> virt_previous_nodes;
                virt_previous_nodes.resize(9);

                // Outer normal at current point
                float virt_outer_normal[3];

                // Number of outer characteristics
                int virt_outer_count = prepare_node(virt_node, virt_node.getRheologyMatrix(),
                                                    time_step, stage, virtMesh,
                                                    virt_dksi, virt_inner, virt_previous_nodes,
                                                    virt_outer_normal);

                // FIXME_ASAP: WA
                switch (stage) {
                case 0: virt_node.getRheologyMatrix()->decomposeX(virt_node);
                    break;
                case 1: virt_node.getRheologyMatrix()->decomposeY(virt_node);
                    break;
                case 2: virt_node.getRheologyMatrix()->decomposeZ(virt_node);
                    break;
                }
                
                // WA for sharp edges
                if(virt_outer_count == 0) {
                    RheologyMatrixPtr curM = cur_node.getRheologyMatrix();
                    RheologyMatrixPtr virtM = virt_node.getRheologyMatrix();
                    int sign = 0;
                    for(int i = 0; i < 9; i++) {
                        if(!inner[i])
                            sign = (curM->getL(i,i) > 0 ? 1 : -1);
                    }
                    for(int i = 0; i < 9; i++) {
                        if( virtM->getL(i,i) * sign < 0 )
                            virt_inner[i] = false;
                    }
                    virt_outer_count = 3;
                }

                // TODO - merge this condition with the next ones
                if (virt_outer_count != 3) {
                    LOG_DEBUG("EXTENDED DEBUG INFO BEGINS");
                    prepare_node(virt_node, virt_node.getRheologyMatrix(), time_step, stage, virtMesh, 
                                 virt_dksi, virt_inner, virt_previous_nodes, virt_outer_normal, true);
                    LOG_DEBUG("EXTENDED DEBUG INFO ENDS");
                    LOG_DEBUG("Calc contact failed. Mesh: " << mesh->getId()
                          << " Virt mesh: " << virtMesh->getId()
                          << "\nReal node: " << cur_node << "\nVirt node: " << virt_node);
                    LOG_DEBUG("There are " << virt_outer_count << " 'outer' characteristics for virt node.");
                    for (int z = 0; z < 9; z++) {
                        LOG_DEBUG("Dksi[" << z << "]: " << virt_dksi[z]);
                        LOG_DEBUG("Inner[" << z << "]: " << virt_inner[z]);
                        LOG_DEBUG("PrNodes[" << z << "]: " << virt_previous_nodes[z]);
                    }
                    THROW_BAD_METHOD("Illegal number of outer characteristics");
                }

                //                // Check that 'paired node' is in the direction of 'outer' characteristics
                //                // If it is not the case - we have strange situation when
                //                // we replace 'outer' points data with data of 'paired node' from different axis direction.
                //
                //                // For all characteristics of real node and virt node
                //                /*for(int i = 0; i < 9; i++)
                //                {
                //                    float v_x_outer[3];
                //                    float v_x_virt[3];
                //                    // Real node - if characteristic is 'outer'*/
                //    /*                if(!inner[i])
                //                    {
                //                        // Find directions to corresponding 'outer' point and to virt 'paired node'
                //                        for(int j = 0; j < 3; j++) {
                //                            v_x_outer[j] = previous_nodes[ppoint_num[i]].coords[j] - cur_node.coords[j];
                //                            v_x_virt[j] = virt_node.coords[j] - cur_node.coords[j];
                //                        }
                //                        // If directions are different - smth bad happens
                //                        if( (v_x_outer[0] * v_x_virt[0]
                //                             + v_x_outer[1] * v_x_virt[1] + v_x_outer[2] * v_x_virt[2]) < 0 )
                //                        {
                //                            *logger << "MESH " << mesh->zone_num << "REAL NODE " << cur_node.local_num << ": "
                //                                    << "x: " << cur_node.coords[0]
                //                                    << " y: " << cur_node.coords[1]
                //                                    << " z: " < cur_node.coords[2];
                //                            log_node_diagnostics(cur_node, stage, outer_normal, mesh, basis_num, rheologyMatrix, time_step, previous_nodes, ppoint_num, inner, dksi);
                //                            *logger << "'Outer' direction: " << v_x_outer[0] << " "
                //                                << v_x_outer[1] << " " < v_x_outer[2];
                //                            *logger << "'Virt' direction: " << v_x_virt[0] << " "
                //                                << v_x_virt[1] << " " < v_x_virt[2];
                //                            throw GCMException( GCMException::METHOD_EXCEPTION, "Bad contact from real node point of view: 'outer' and 'virt' directions are different");
                //                        }
                //                    }*/
                //    // We switch it off because it conflicts sometimes with 'safe_direction'
                //    /*                // Virt node - if characteristic is 'outer'
                //                    if(!virt_inner[i])
                //                    {
                //                        // Find directions to corresponding 'outer' point and to real 'paired node'
                //                        for(int j = 0; j < 3; j++) {
                //                            v_x_outer[j] = virt_previous_nodes[virt_ppoint_num[i]].coords[j] - virt_node.coords[j];
                //                            v_x_virt[j] = cur_node.coords[j] - virt_node.coords[j];
                //                        }
                //                        // If directions are different - smth bad happens
                //                        if( (v_x_outer[0] * v_x_virt[0]
                //                            + v_x_outer[1] * v_x_virt[1] + v_x_outer[2] * v_x_virt[2]) < 0 )
                //                        {
                //                            *logger << "MESH " << mesh->zone_num << "REAL NODE " << cur_node.local_num << ": "
                //                                    << "x: " << cur_node.coords[0]
                //                                    << " y: " << cur_node.coords[1]
                //                                    << " z: " < cur_node.coords[2];
                //                            log_node_diagnostics(virt_node, stage, virt_outer_normal, virt_node.mesh, basis_num, virt_rheologyMatrix, time_step, virt_previous_nodes, virt_ppoint_num, virt_inner, virt_dksi);
                //                            *logger << "'Outer' direction: " << v_x_outer[0] << " "
                //                                << v_x_outer[1] << " "< v_x_outer[2];
                //                            *logger << "'Virt' direction: " << v_x_virt[0] << " "
                //                                << v_x_virt[1] << " " < v_x_virt[2];
                //                            throw GCMException( GCMException::METHOD_EXCEPTION, "Bad contact from virt node point of view: 'outer' and 'virt' directions are different");
                //                        }
                //                    }*/
                ////                }
                //
                LOG_TRACE("Using calculator: " << engine.getContactCondition(cur_node.getContactConditionId())->calc->getType());
                LOG_TRACE("Outer normal: " << outer_normal[0] << " " << outer_normal[1] << " " << outer_normal[2]);
                engine.getContactCondition(cur_node.getContactConditionId())->doCalc(Engine::getInstance().getCurrentTime(), cur_node,
                                                       new_node, virt_node, cur_node.getRheologyMatrix(), previous_nodes, inner,
                                                       virt_node.getRheologyMatrix(), virt_previous_nodes, virt_inner, outer_normal);
            }
            // It means smth went wrong. Just interpolate the values and report bad node.
        }
        else
        {
            //LOG_WARN("Outer count: " << outer_count);
            //LOG_WARN("Node: " << cur_node);
            //for (int z = 0; z < 9; z++) {
            //    LOG_WARN("Dksi[" << z << "]: " << dksi[z]);
            //    LOG_WARN("Inner[" << z << "]: " << inner[z]);
            //    LOG_WARN("PrNodes[" << z << "]: " << previous_nodes[z]);
            //}
            //THROW_BAD_METHOD("Illegal number of outer characteristics");
            // FIXME - implement border and contact completely
            LOG_TRACE("Using calculator: " << engine.getBorderCondition(0)->calc->getType());
            engine.getBorderCondition(0)->doCalc(Engine::getInstance().getCurrentTime(), cur_node,
                                                  new_node, cur_node.getRheologyMatrix(), previous_nodes, inner, outer_normal);
            //cur_node.setNeighError(stage);
        }
        LOG_TRACE("Done border node calc");
    }
}

void InterpolationFixedAxis::doNextPartStep(CalcNode& cur_node, CalcNode& new_node, float time_step, int stage, Mesh* mesh)
{
    TRACE_ON_EXCEPTION(__doNextPartStep(cur_node, new_node, time_step, stage, mesh));
}

int InterpolationFixedAxis::prepare_node(CalcNode& cur_node, RheologyMatrixPtr rheologyMatrix,
                                              float time_step, int stage, Mesh* mesh,
                                              float* dksi, bool* inner, vector<CalcNode>& previous_nodes,
                                              float* outer_normal)
{
    return prepare_node(cur_node, rheologyMatrix, time_step, stage, mesh, dksi, inner, previous_nodes, outer_normal, false);
}

int InterpolationFixedAxis::prepare_node(CalcNode& cur_node, RheologyMatrixPtr rheologyMatrix,
                                              float time_step, int stage, Mesh* mesh,
                                              float* dksi, bool* inner, vector<CalcNode>& previous_nodes,
                                              float* outer_normal, bool debug)
{
    assert_ge(stage, 0);
    assert_le(stage, 2);

    if (cur_node.isBorder())
        mesh->findBorderNodeNormal(cur_node, &outer_normal[0], &outer_normal[1], &outer_normal[2], false);

    LOG_TRACE("Preparing elastic matrix");
    //  Prepare matrixes  A, Lambda, Omega, Omega^(-1)

    switch (stage) {
    case 0: rheologyMatrix->decomposeX(cur_node);
        break;
    case 1: rheologyMatrix->decomposeY(cur_node);
        break;
    case 2: rheologyMatrix->decomposeZ(cur_node);
        break;
    }
    LOG_TRACE("Preparing elastic matrix done");

    LOG_TRACE("Elastic matrix eigen values:\n" << rheologyMatrix->getL());

    for (int i = 0; i < 9; i++)
        dksi[i] = -rheologyMatrix->getL(i, i) * time_step;

    return find_nodes_on_previous_time_layer(cur_node, stage, mesh, dksi, inner, previous_nodes, outer_normal, debug);
}

int InterpolationFixedAxis::find_nodes_on_previous_time_layer(CalcNode& cur_node, int stage, Mesh* mesh,
                                                                   float dksi[], bool inner[], vector<CalcNode>& previous_nodes,
                                                                   float outer_normal[])
{
    return find_nodes_on_previous_time_layer(cur_node, stage, mesh, dksi, inner, previous_nodes, outer_normal, false);
}

int InterpolationFixedAxis::find_nodes_on_previous_time_layer(CalcNode& cur_node, int stage, Mesh* mesh,
                                                                   float dksi[], bool inner[], vector<CalcNode>& previous_nodes,
                                                                   float outer_normal[], bool debug)
{
    LOG_TRACE("Start looking for nodes on previous time layer");
    
    // For all omegas
    for (int i = 0; i < 9; i++) {
        LOG_TRACE("Looking for characteristic " << i);
        // Check prevoius omegas ...
        bool already_found = false;
        for (int j = 0; j < i; j++) {
            // ... And try to find if we have already worked with the required point
            // on previous time layer (or at least with the point that is close enough)
            if (fabs(dksi[i] - dksi[j]) <= EQUALITY_TOLERANCE * 0.5 * fabs(dksi[i] + dksi[j])) {
                LOG_TRACE("Found old value " << dksi[i] << " - done");
                // If we have already worked with this point - just remember the number
                already_found = true;
                previous_nodes[i] = previous_nodes[j];
                inner[i] = inner[j];
            }
        }

        // If we do not have necessary point in place - ...
        if (!already_found) {
            LOG_TRACE("New value " << dksi[i] << " - preparing vectors");
            // ... Put new number ...
            previous_nodes[i] = cur_node;

            // ... Find vectors ...
            float dx[3];
            dx[0] = dx[1] = dx[2] = 0.0;
            dx[stage] += dksi[i];

            // For dksi = 0 we can skip check and just copy everything
            if (dksi[i] == 0) {
                // no interpolation required - everything is already in place
                inner[i] = true;
                LOG_TRACE("dksi is zero - done");
            }
            else if (cur_node.isInner()) {
                LOG_TRACE("Checking inner node");
                // ... Find owner tetrahedron ...
                bool isInnerPoint;
                mesh->interpolateNode(cur_node, dx[0], dx[1], dx[2], debug,
                                      previous_nodes[i], isInnerPoint);

                if (!isInnerPoint) {
                    LOG_TRACE("Inner node: we need new method here!");
                    LOG_TRACE("Node:\n" << cur_node);
                    LOG_TRACE("Move: " << dx[0] << " " << dx[1] << " " << dx[2]);
                    // TODO: return it back later
                    // Re-run search with debug on
                    //mesh->interpolateNode(origin, dx[0], dx[1], dx[2], true,
                    //                      previous_nodes[i], isInnerPoint);
                }

                inner[i] = true;
                LOG_TRACE("Checking inner node done");
            }
            else if (cur_node.isBorder()) {
                LOG_TRACE("Checking border node");
                // ... Find owner tetrahedron ...
                bool isInnerPoint;
                bool interpolated = mesh->interpolateNode(cur_node, dx[0], dx[1], dx[2], debug,
                                                          previous_nodes[i], isInnerPoint);

                // If we found inner point, it means
                // this direction is inner and everything works as for usual inner point
                if (isInnerPoint) {
                    inner[i] = true;
                    // If we did not find inner point - two cases are possible
                }
                else {
                    inner[i] = false;
                    // We found border cross somehow
                    // It can happen if we work with really thin structures and big time step
                    // We can work as usual in this case
                    if (interpolated) {
                        LOG_TRACE("Border node: we need new method here!");
                        inner[i] = true;
                        // Or we did not find any point at all - it means this characteristic is outer
                    }
                    else {
                        inner[i] = false;
                    }
                }
                LOG_TRACE("Checking border node done");
            }
            else {
                THROW_BAD_MESH("Unsupported case for characteristic location");
            }
        }
        LOG_TRACE("Looking for characteristic " << i << " done");
    }

    int outer_count = 0;
    for (int i = 0; i < 9; i++)
        if (!inner[i])
            outer_count++;

    // assert_true(outer_count == 0 || outer_count == 3);

    LOG_TRACE("Looking for nodes on previous time layer done. Outer count = " << outer_count);

    return outer_count;
}

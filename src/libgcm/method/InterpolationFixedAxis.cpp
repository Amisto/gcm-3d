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
    //LOG_INFO("QB 00");
    Engine &e = Engine::getInstance();
    RheologyMatrixPtr curM = cur_node.getRheologyMatrix();
    //LOG_INFO("QB 01");
    gcm::real distq = (cur_node.coords[0] - virt_node.coords[0])*(cur_node.coords[0] - virt_node.coords[0]) + (cur_node.coords[1] - virt_node.coords[1])*(cur_node.coords[1] - virt_node.coords[1]) + (cur_node.coords[2] - virt_node.coords[2])*(cur_node.coords[2] - virt_node.coords[2]);
    curM->decompose(cur_node, abs(dir)-1);
    int n_outer = 0;
    for (int i=0; i<9; i++)
    {
	//LOG_INFO("QB 02 "<<i);
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
	if (curM->getL(i, i)*dir < -EQUALITY_TOLERANCE)
	    inner[i] = false;
	if (!inner[i]) n_outer++;
    }
    if (n_outer != 3) LOG_INFO("QB 03 " <<n_outer);	
    float outer_normal[3] = {0};
    outer_normal[abs(dir)-1] = dir/abs(dir);
    cur_node.contactDirection = 1;
    e.getBorderCondition(cur_node.getBorderConditionId())->doCalc(e.getCurrentTime(), cur_node,
                                                                new_node, cur_node.getRheologyMatrix(), previous_nodes, inner, outer_normal);
    if (!cur_node.contactDirection)
    {
	LOG_INFO("QB \n cur \n"<<cur_node.getRheologyMatrix()->getU());
	LOG_INFO("o_n "<<outer_normal[0]<<" "<<outer_normal[1]<<" "<<outer_normal[2]<<" ");
    }
    //LOG_INFO("QB 04");
}

void InterpolationFixedAxis::quasi_contact(CalcNode& cur_node, CalcNode& new_node, CalcNode& virt_node, CalcNode& virt_node_inner, Mesh* vmi, vector<CalcNode>& previous_nodes, bool* inner, vector<CalcNode>& virt_previous_nodes, bool* virt_inner, float time_step, int dir)
{
    Engine &e = Engine::getInstance();
    RheologyMatrixPtr curM = cur_node.getRheologyMatrix();
    gcm::real distq = (cur_node.coords[0] - virt_node_inner.coords[0])*(cur_node.coords[0] - virt_node_inner.coords[0]) + (cur_node.coords[1] - virt_node_inner.coords[1])*(cur_node.coords[1] - virt_node_inner.coords[1]) + (cur_node.coords[2] - virt_node_inner.coords[2])*(cur_node.coords[2] - virt_node_inner.coords[2]);
    int n_outer = 0;
    for (int i=0; i<9; i++)
    {
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
	if (curM->getL(i, i)*dir < -EQUALITY_TOLERANCE)
            inner[i] = false;
        if (!inner[i]) n_outer++;
    }
    if (n_outer != 3) LOG_INFO("QC " <<n_outer);
    float outer_normal[3] = {0};
    outer_normal[abs(dir)-1] = dir/abs(dir);
    //LOG_INFO("QC \n cur \n"<<cur_node.getRheologyMatrix()->getU());
    //LOG_INFO("virt \n"<<virt_node.getRheologyMatrix()->getU());
    //LOG_INFO("o_n "<<outer_normal[0]<<" "<<outer_normal[1]<<" "<<outer_normal[2]<<" ");
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
	//LOG_INFO("outers: "<<outer_count);
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
	//LOG_INFO("outers: "<<outer_count);
        //if (outer_count == 1)
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
	if (outer_count > 3)
        {
            quasi_zerofy();
            return;
        }
	if (outer_count != 3) LOG_INFO("outers: "<<outer_count);
	//LOG_INFO("02");
/*	if (outer_count == 2 || outer_count == 4 || outer_count == 6) //double contact, we have another bodies on both sides
	{
	    //firstly, get virtual nodes and check materials 
	    CalcNode v_n_l, v_n_r;
            CalcNode &virt_node_left = v_n_l, &virt_node_right = v_n_r;
	    Mesh *ml = NULL, *mr = NULL;
	    virt_node_left.number = 0;
	    virt_node_right.number = 0;
	    //RheologyMatrixPtr lM, rM;
//	    LOG_INFO("L");
	    engine.getVirtNode(cur_node, virt_node_left, -stage-1, ml);
//	    LOG_INFO("R");
            engine.getVirtNode(cur_node, virt_node_right, stage+1, mr);
//	    LOG_INFO("D");
            if (virt_node_left.number) virt_node_left.setInContact(true);
	    if (virt_node_right.number) virt_node_right.setInContact(true);
            float virt_left_dksi[9], virt_right_dksi[9];
            bool virt_left_inner[9], virt_right_inner[9];
            vector<CalcNode> virt_left_previous_nodes, virt_right_previous_nodes;
            virt_left_previous_nodes.resize(9);
            virt_right_previous_nodes.resize(9);
            float virt_left_outer_normal[3] = {0}, virt_right_outer_normal[3] = {0};
	    virt_left_outer_normal[stage] = 1;
	    virt_right_outer_normal[stage] = -1;
            int virt_left_outer_count;
	    RheologyMatrixPtr vlM, vrM;
	    if (virt_node_left.number) 
	    {
		ml = engine.getBody(virt_node_left.contactDirection)->getMeshes();
		vlM = virt_node_left.getRheologyMatrix();
		virt_node_left.setIsBorder(false);
		virt_left_outer_count = prepare_node(virt_node_left, vlM,
                                                    time_step, stage, ml,
                                                    virt_left_dksi, virt_left_inner, virt_left_previous_nodes,
                                                    virt_left_outer_normal);
		virt_node_left.setIsBorder(true);
	    }
            int virt_right_outer_count;
	    if (virt_node_right.number) 
	    {
		mr = engine.getBody(virt_node_right.contactDirection)->getMeshes();
		vrM = virt_node_right.getRheologyMatrix();
		virt_node_right.setIsBorder(false);
		virt_right_outer_count = prepare_node(virt_node_right, vrM,
                                                    time_step, stage, mr,
                                                    virt_right_dksi, virt_right_inner, virt_right_previous_nodes,
                                                    virt_right_outer_normal);
		virt_node_right.setIsBorder(true);
	    }
	    //LOG_INFO("Prepared");
	    //Check Direction for virt nodes
//	    LOG_INFO("P");
	    float valL = 0;
	    if (virt_node_left.number)
	    {
		//vlM = virt_node_left.getRheologyMatrix();
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
            	//vrM = virt_node_right.getRheologyMatrix();
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
		    //LOG_INFO("LR 0");
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
		    //LOG_INFO("LR 1");
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
	} */
	//LOG_INFO("06");
	if (outer_count == 3)
	{
	    //LOG_INFO("07");
	    CalcNode v_n;
	    CalcNode &virt_node = v_n;
	    virt_node.number = 0;
            Mesh *m = NULL;
            //RheologyMatrixPtr lM, rM;
//	    LOG_INFO("S");
            engine.getVirtNode(cur_node, virt_node, val*(stage+1), m); 
	    //LOG_INFO("Got virt node "<<virt_node.number <<" body "<<(int)virt_node.contactDirection <<" mID " <<(int)virt_node.getMaterialId());
	    if (virt_node.number) 
	    {
		m = engine.getBody(virt_node.contactDirection)->getMeshes();
	    	if (!m) LOG_INFO("No mesh for virt node");
	    }
//	    LOG_INFO("D");
            if (virt_node.number) virt_node.setInContact(true);
            float virt_dksi[9];
            bool virt_inner[9];
            vector<CalcNode> virt_previous_nodes;
            virt_previous_nodes.resize(9);
            float virt_outer_normal[3] = {-outer_normal[0], -outer_normal[1], -outer_normal[2]};
            int virt_outer_count;
            if (virt_node.number)
	    { 
//	  	LOG_INFO("08 " <<(int)virt_node.getMaterialId());
		//cur_node.setIsBorder(false);
		RheologyMatrixPtr vrm = virt_node.getRheologyMatrix();
		//vrm->decomposeX(virt_node);
		virt_node.setIsBorder(false);
		virt_outer_count = prepare_node(virt_node, vrm,
                                                    time_step, stage, m,
                                                    virt_dksi, virt_inner, virt_previous_nodes,
                                                    virt_outer_normal);
		virt_node.setIsBorder(true);
//		LOG_INFO("08 between");
            	RheologyMatrixPtr virtM = virt_node.getRheologyMatrix();
            	float valL = 0;
		virt_outer_count = 0;
//	 	LOG_INFO("080");
            	for (uint i=0; i<9; i++)
            	{
//		    LOG_INFO("081 " <<i);
                    valL = (virtM->getL(i, i) > 0 ? -1.0 : (fabs(virtM->getL(i, i)) < EQUALITY_TOLERANCE ? 0 : 1.0));
                    if (valL == 0) continue;
                    if (val != valL) virt_inner[i] = false;
		    if (!virt_inner[i]) virt_outer_count++;
            	}
//	 	LOG_INFO("082 " <<virt_outer_count);
		if (virt_outer_count < 3) LOG_INFO("OBVIOUSLY, we fucked up smt 1 " <<virt_outer_count);
		if (virt_outer_count > 3) //quasi_border(); //FIXME - cut an unlikely branch of algorhythm
		    quasi_border(cur_node, new_node, virt_node, m, (stage+1)*val, previous_nodes, inner, time_step);
		//LOG_INFO("09");
    		//LOG_INFO("cur \n"<<cur_node.getRheologyMatrix()->getU());
    		//LOG_INFO("virt \n"<<virt_node.getRheologyMatrix()->getU());
    		//LOG_INFO("o_n "<<outer_normal[0]<<" "<<outer_normal[1]<<" "<<outer_normal[2]<<" ");
		//cur_node.getRheologyMatrix()->decompose(cur_node, stage);
		//virt_node.getRheologyMatrix()->decompose(virt_node, stage);
		//LOG_INFO("TC "<<outer_normal[0] <<" "<<outer_normal[1] <<" "<<outer_normal[2] <<" ");
		virt_node.number = 1;
		if (virt_outer_count == 3) engine.getContactCondition(cur_node.getContactConditionId())->doCalc(Engine::getInstance().getCurrentTime(), cur_node,
                                                       new_node, virt_node, cur_node.getRheologyMatrix(), previous_nodes, inner,
                                                       virt_node.getRheologyMatrix(), virt_previous_nodes, virt_inner, outer_normal);  
		if (!virt_node.number)
		{
		    LOG_INFO("09");
                    LOG_INFO("cur \n"<<cur_node.getRheologyMatrix()->getU());
                    LOG_INFO("virt \n"<<virt_node.getRheologyMatrix()->getU());
                    LOG_INFO("o_n "<<outer_normal[0]<<" "<<outer_normal[1]<<" "<<outer_normal[2]<<" ");
		}
	 	//LOG_INFO("10");
		return;
	    }
	    //LOG_INFO("11 " <<(int)cur_node.getBorderConditionId());
	    if (cur_node.getBorderConditionId())
	    	engine.getBorderCondition(cur_node.getBorderConditionId())->doCalc(Engine::getInstance().getCurrentTime(), cur_node,
                                                                 new_node, cur_node.getRheologyMatrix(), previous_nodes, inner, outer_normal);
	    else
	    	engine.getBorderCalculator("FreeBorderCalculator")->doCalc(cur_node, new_node, cur_node.getRheologyMatrix(), previous_nodes, inner, outer_normal, 1.0);
	    //LOG_INFO("12");
	    return;
	}

	//LOG_INFO("OBVIOUSLY, we forgot something 3 " <<outer_count);

	return;
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
//    LOG_INFO("PN 00");
    if (cur_node.isBorder())
        mesh->findBorderNodeNormal(cur_node, &outer_normal[0], &outer_normal[1], &outer_normal[2], false);
//    LOG_INFO("PN 01");
    LOG_TRACE("Preparing elastic matrix");
    //  Prepare matrixes  A, Lambda, Omega, Omega^(-1)
//    LOG_INFO("Elastic matrix:\n" << rheologyMatrix);
//    LOG_INFO("Elastic matrix eigen values: before \n" << rheologyMatrix->getL());
    switch (stage) {
    case 0: 
	//LOG_INFO("PN 01 0 " <<stage);
	rheologyMatrix->decomposeX(cur_node);
        break;
    case 1: 
	//LOG_INFO("PN 01 1 " <<stage);
	rheologyMatrix->decomposeY(cur_node);
        break;
    case 2: 
	//LOG_INFO("PN 01 2 " <<stage);
	rheologyMatrix->decomposeZ(cur_node);
        break;
    }
//    LOG_INFO("PN 02");
//    LOG_INFO("Preparing elastic matrix done");

//    LOG_INFO("Elastic matrix eigen values: after \n" << rheologyMatrix->getL());

    for (int i = 0; i < 9; i++)
        dksi[i] = -rheologyMatrix->getL(i, i) * time_step;
//    LOG_INFO("PN 03");
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
            else //if (cur_node.isBorder()) 
	    {
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
                else 
		{
		    Engine& e = Engine::getInstance();
		    inner[i] = false;
		    Mesh* m = NULL;
		    for (int j=0; j<e.getNumberOfBodies(); j++)	
		    {
			if (inner[i]) continue;
			m = e.getBody(j)->getMeshes();
			if (!m) continue;
 			bool found = m->interpolateBorderNode(cur_node.coords[0], cur_node.coords[1], cur_node.coords[2], dx[0], dx[1], dx[2], previous_nodes[i]);
			if (found) 
			{
			    inner[i] = true;
			    //LOG_INFO("FOUND A NODE");
			}
	   	    }
                    //inner[i] = false;
                    // We found border cross somehow
                    // It can happen if we work with really thin structures and big time step
                    // We can work as usual in this case
                    //if (interpolated) {
                    //    LOG_TRACE("Border node: we need new method here!");
                    //    inner[i] = true;
                        // Or we did not find any point at all - it means this characteristic is outer
                    //}
                    //else {
                    //    inner[i] = false;
                    //}
                }
                LOG_TRACE("Checking border node done");
            }
            //else {
            //    THROW_BAD_MESH("Unsupported case for characteristic location");
            //}
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

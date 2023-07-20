#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "VaultMacroElement_2D.h"
#include <Channel.h>
#include <CrdTransf.h>
#include <CompositeResponse.h>
#include <Domain.h>
#include <ID.h>
#include <Information.h>
#include <Node.h>
#include <Parameter.h>
#include <Renderer.h>
#include <TransientIntegrator.h>
#include <Vector.h>
#include <elementAPI.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <FEM_ObjectBroker.h>
#include <Matrix.h>
#include <SectionForceDeformation.h>
#include <ElasticSection2d.h> 
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <limits>
#include <string>
#include <cmath>
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
#ifndef DBL_EPSILON
#define DBL_EPSILON (std::numeric_limits<double>::epsilon())
#endif
#include <typeinfo>
#include <string>

// -----This from ME: -------------------
#ifdef _USRDLL
#include <windows.h>
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif
//------------------------------------------

Matrix VaultMacroElement_2d::K(6,6); //Iz ME je 18 x 18, u 2d bi trebalo ba bude 6x6
Vector VaultMacroElement_2d::P(6);    // Iz ME je 18 - ti ovde imas 2 node = 6
double VaultMacroElement_2d::workArea[100];

// TODO: Is it needed? Da, za neku error -------------------
static int numMyMacroelement = 0;
//------------------------------------------

// ------This from VaultMacroElement: -----------------------------------
// const std::string VaultMacroElement_2d::acc_file_name = "acc.csv"; //TODO: promeni u accel koji dolaze iz ME?
const std::string VaultMacroElement_2d::result_file_name = "result.csv";
const std::string VaultMacroElement_2d::python_script_name = "VaultMacroElement.py";


// Export acc into a csv
// void VaultMacroElement_2d::setAcc(std::vector<std::vector<double>> acc) {
//	std::ofstream acc_file(acc_file_name);

//	if (!acc_file.is_open()) {
//		std::cerr << "Could not open " << acc_file_name << std::endl;
//		return;
//	}

//	for (int i = 0; i < acc.size(); i++) {
//		for (int j = 0; j < acc[i].size(); j++) {
//			acc_file << acc[i][j] << ",";
//		}
//		acc_file << std::endl;
//	}

//	acc_file.close();
//}
//--------------------------------------------------------------------------


OPS_Export void *
OPS_VaultMacroElement_2d()
{
	  //Vector intLength(2);
    //intLength(0) = 1.0/2;
    //intLength(1) = 1.0/2;
    //Vector intLengthMasses(2);
    //intLengthMasses(0) = 0.50;
    //intLengthMasses(1) = 0.50;
        
    if (numMyMacroelement == 0) {
          opserr << "VaultMacroElement_2d loaded from external library - Written by IB, EPFL, 2023" << endln;
          numMyMacroelement++;
    }

    int remaining = OPS_GetNumRemainingInputArgs();

    if (remaining < 7) { 
		    opserr << "WARNING: insufficient arguments for element geometry." << endln;
		    opserr << "Required: eleTag, iNode, jNode, transfTag, integrationTag, dt, theta <-flags>" << endln;
	  return 0;
    }

    // inputs: 
    int iData[5]; //eletag, iNode, jNode, transfTag, integrationTag
    int numData = 5;
    if(OPS_GetIntInput(&numData,&iData[0]) < 0) {
	  opserr<<"WARNING: invalid integer inputs, element geometry definition." << endln;
	  return 0;
    }

    double dData[2]; //dt, theta
    numData = 2;
    if(OPS_GetDoubleInput(&numData,&dData[0]) < 0) {
	  opserr<<"WARNING: invalid double inputs, element geometry definition (tag: "<< iData[0] << ")" << endln;
	  return 0;
    }

    double dt = dData[0];
    doubletheta = dData[1];


    // options
    double mass = 0.0;
    int cmass = 0;
    numData = 1;
    while(OPS_GetNumRemainingInputArgs() > 0) {
	const char* type = OPS_GetString();
	if(strcmp(type, "-cMass") == 0) {
	    cmass = 1;
	} else if(strcmp(type,"-mass") == 0) {
	    if(OPS_GetNumRemainingInputArgs() > 0) {
		if(OPS_GetDoubleInput(&numData,&mass) < 0) {
		    opserr<<"WARNING: invalid mass\n";
		    return 0;
		}
	    }
	}
    }

    int dampingModel = 0;

    // check transf
    CrdTransf* theTransf = OPS_getCrdTransf(iData[3]);
    if(theTransf == 0) {
	opserr<<"coord transfomration not found\n";
	return 0;
    }

    // check beam integrataion
    BeamIntegrationRule* theRule = OPS_getBeamIntegrationRule(iData[4]);
    if(theRule == 0) {
	opserr<<"beam integration not found\n";
	return 0;
    }
    BeamIntegration* bi = theRule->getBeamIntegration();
    if(bi == 0) {
	opserr<<"beam integration is null\n";
	return 0;
    }

    
    // check sections
    const ID& secTags = theRule->getSectionTags();
    SectionForceDeformation** sections = new SectionForceDeformation *[secTags.Size()];
    for(int i=0; i<secTags.Size(); i++) {
	sections[i] = OPS_getSectionForceDeformation(secTags(i));
	if(sections[i] == 0) {
	    opserr<<"section "<<secTags(i)<<"not found\n";
		delete [] sections;
	    return 0;
	}
    }


    //TODO: Da li dt i theta napisati imenima ili dData[0], dData[1]
    Element *theEle =  new VaultMacroElement_2d(iData[0],iData[1],iData[2],secTags.Size(),sections,
					    *bi,*theTransf, dt, theta, mass, cmass, dampingModel);
    delete [] sections;
    return theEle;


}



VaultMacroElement_2d::VaultMacroElement_2d(int tag, int nd1, int nd2,
				   int numSec, SectionForceDeformation **s,
				   BeamIntegration& bi,
				   CrdTransf &coordTransf, double dt, double theta, double r, int cm, int dampingModel)
:Element (tag, ELE_TAG_VaultMacroElement_2d), 
 numSections(numSec), theSections(0), crdTransf(0), beamInt(0),
  connectedExternalNodes(2),
  Q(6), q(3),  dt(dt), theta(theta), rho(r), cMass(cm), parameterID(0), dampingModel(dampingModel)
{


  // Allocate arrays of pointers to SectionForceDeformations
  theSections = new SectionForceDeformation *[numSections];
    
  if (theSections == 0) {
    opserr << VaultMacroElement_2d::VaultMacroElement_2d- failed to allocate section model pointer\n";
    exit(-1);
  }

  for (int i = 0; i < numSections; i++) {
    
    // Get copies of the material model for each integration point
    theSections[i] = s[i]->getCopy();
    
    // Check allocation
    if (theSections[i] == 0) {
      opserr << "VaultMacroElement_2d::VaultMacroElement_2d -- failed to get a copy of section model\n";
      exit(-1);
    }
  }
  
  beamInt = bi.getCopy();
  
  if (beamInt == 0) {
    opserr << "VaultMacroElement_2d::VaultMacroElement_2d - failed to copy beam integration\n";
    exit(-1);
  }

  crdTransf = coordTransf.getCopy2d();
  
  if (crdTransf == 0) {
    opserr << "VaultMacroElement_2d::VaultMacroElement_2d- failed to copy coordinate transformation\n";
    exit(-1);
  }
  
  // Set connected external node IDs
  connectedExternalNodes(0) = nd1;
  connectedExternalNodes(1) = nd2;
  
  theNodes[0] = 0;
  theNodes[1] = 0;
  
  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  
  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
   
    
}

VaultMacroElement_2d::VaultMacroElement_2d()
:Element (tag, ELE_TAG_VaultMacroElement_2d), 
 numSections(numSec), theSections(0), crdTransf(0), beamInt(0),
  connectedExternalNodes(2),
  Q(6), q(3),  dt(dt), theta(theta), rho(r), cMass(cm), parameterID(0), dampingModel(dampingModel)
{
    q0[0] = 0.0;
    q0[1] = 0.0;
    q0[2] = 0.0;

    p0[0] = 0.0;
    p0[1] = 0.0;
    p0[2] = 0.0;

    theNodes[0] = 0;
    theNodes[1] = 0;
}




VaultMacroElement_2d::~VaultMacroElement_2d()
{    
  for (int i = 0; i < numSections; i++) {
    if (theSections[i])
      delete theSections[i];
  }
  
  // Delete the array of pointers to SectionForceDeformation pointer arrays
  if (theSections)
    delete [] theSections;
  
  if (crdTransf)
    delete crdTransf;

  if (beamInt != 0)
    delete beamInt;
}


int
VaultMacroElement_2d::getNumExternalNodes()const {
    return 2;
}

const ID&
VaultMacroElement_2d::getExternalNodes() {
    return connectedExternalNodes;
} 

Node **
VaultMacroElement_2d::getNodePtrs() {
    return theNodes;
}

int
VaultMacroElement_2d::getNumDOF() {
    return 6;
}



void
VaultMacroElement_2d::setDomain(Domain *theDomain) {
	// Check Domain is not null - invoked when object removed from a domain

	  if (theDomain == 0) {
  theNodes[0] = 0;
  theNodes[1] = 0;
  return;
    }

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    

    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);
    
    if (theNodes[0] == 0 || theNodes[1] == 0) {
      opserr << "WARNING VaultMacroElement_2d (tag: %d), node not found in domain" << this->getTag() << endln;;
      return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();

	  if (dofNd1 != 3 || dofNd2 != 3) {
	return;
    }

  if (crdTransf->initialize(theNodes[0], theNodes[1])) {
		// Add some error check
	}

	double L = crdTransf->getInitialLength();

	if (L == 0.0) {
		// Add some error check
	}

    this->DomainComponent::setDomain(theDomain);

	this->update();
}



int
VaultMacroElement_2d::commitState() {

    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
		opserr << "VaultMacroElement_2d:commitState () - failed in base class";
    }    
    

    
	// Loop over the integration points and commit the material states
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->commitState();

    retVal += crdTransf->commitState();

    return retVal;
}



int
VaultMacroElement_2d::update(void)
{
  int err = 0;

  // Update the transformation
  crdTransf->update();
  
  /// Get basic deformations
  const Vector &v = crdTransf->getBasicTrialDisp();
  const Vector &v = crdTransf->getBasicTrialAccel();
  

   //FIXME: IT4R: Pass to my python getBasicTrialAcc/Displ 



  //TODO: Lines below not needed.
  //double L = crdTransf->getInitialLength();
  //double oneOverL = 1.0/L;

  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  //double xi[maxNumSections];
  //beamInt->getSectionLocations(numSections, L, xi);
  
  // Loop over the integration points
  //for (int i = 0; i < numSections; i++) {
    
   // int order = theSections[i]->getOrder();
   // const ID &code = theSections[i]->getType();
    
   // Vector e(workArea, order);
    
    //double xi6 = 6.0*pts(i,0);
   // double xi6 = 6.0*xi[i];
    
   // int j;
   // for (j = 0; j < order; j++) {
   //   switch(code(j)) {
  //    case SECTION_RESPONSE_P:
	//e(j) = oneOverL*v(0); break;
   //   case SECTION_RESPONSE_MZ:
	//e(j) = oneOverL*((xi6-4.0)*v(1) + (xi6-2.0)*v(2)); break;
   //   default:
	//e(j) = 0.0; break;
    //  }
   // }
    
    // Set the section deformations
    //err += theSections[i]->setTrialSectionDeformation(e);
 // }
  
 // Run python optimization
//void VaultMacroElement_2d::compute() {
// Call python script
//	std::string command = "python " + python_script_name + " " + std::to_string(dt) + " " + std::to_string(theta);

//    int result = std::system(command.c_str());

//    if (result != 0) {
//        std::cerr << "Python command failed with error code " << result << std::endl;
//    }

//	loadResult();
//}

// Load optimization result (csv) into memory
//void VaultMacroElement_2d::loadResult() {
//	std::ifstream result_file(result_file_name);

//	if (!result_file.is_open()) {
//		std::cerr << "Could not open " << result_file_name << std::endl;
//		return;
//	}

//	result = std::vector<std::vector<double>>(); // clear result currently in memory
//	std::string line;

//	while (std::getline(result_file, line)) {
//		std::stringstream line_stream(line);
//		std::vector<double> line_values;
//		std::string value;

//		while (std::getline(line_stream, value, ',')) {
//			line_values.push_back(std::stod(value));
//		}

//		result.push_back(line_values);
//	}

//	result_file.close();
//}



  if (err != 0) {
    opserr << "VaultMacroElement_2d::update() - failed setTrialSectionDeformations()\n";
    return err;
  }

  return 0;
}


void
VaultMacroElement_2d::getBasicStiff(Matrix &kb, int initial)
{
  // Zero for integral
  kb.Zero();
  
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

    Matrix ka(workArea, order, 3);
    ka.Zero();

    double xi6 = 6.0*xi[i];

    // Get the section tangent stiffness and stress resultant
    const Matrix &ks = theSections[i]->getSectionTangent();
        
    // Perform numerical integration
    //kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);
    //double wti = wts(i)*oneOverL;
    double wti = wt[i]*oneOverL;
    double tmp;
    int j, k;
    for (j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	for (k = 0; k < order; k++)
	  ka(k,0) += ks(k,j)*wti;
	break;
      case SECTION_RESPONSE_MZ:
	for (k = 0; k < order; k++) {
	  tmp = ks(k,j)*wti;
	  ka(k,1) += (xi6-4.0)*tmp;
	  ka(k,2) += (xi6-2.0)*tmp;
	}
	break;
      default:
	break;
      }
    }
    for (j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_P:
	for (k = 0; k < 3; k++)
	  kb(0,k) += ka(j,k);
	break;
      case SECTION_RESPONSE_MZ:
	for (k = 0; k < 3; k++) {
	  tmp = ka(j,k);
	  kb(1,k) += (xi6-4.0)*tmp;
	  kb(2,k) += (xi6-2.0)*tmp;
	}
	break;
      default:
	break;
      }
    }
  }
}



const Matrix&
VaultMacroElement_2d::getTangentStiff()
{
  static Matrix kb(3,3);

  this->getBasicStiff(kb);

  // Zero for integral
  q.Zero();
  
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  
  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

    //double xi6 = 6.0*pts(i,0);
    double xi6 = 6.0*xi[i];

    // Get the section tangent stiffness and stress resultant
    const Vector &s = theSections[i]->getStressResultant();
        
    // Perform numerical integration
    //q.addMatrixTransposeVector(1.0, *B, s, wts(i));
    double si;
    for (int j = 0; j < order; j++) {
      //si = s(j)*wts(i);
      si = s(j)*wt[i];
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	q(0) += si; break;
      case SECTION_RESPONSE_MZ:
	q(1) += (xi6-4.0)*si; q(2) += (xi6-2.0)*si; break;
      default:
	break;
      }
    }
    
  }
  
  // Add effects of element loads, q = q(v) + q0 //TODO:
  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];

  // Transform to global stiffness
  K = crdTransf->getGlobalStiffMatrix(kb, q);

  return K;
}



const Matrix&
VaultMacroElement_2d::getInitialBasicStiff()
{
  static Matrix kb(3,3);

  // Zero for integral
  kb.Zero();
  
  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  
  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();
  
    Matrix ka(workArea, order, 3);
    ka.Zero();

    //double xi6 = 6.0*pts(i,0);
    double xi6 = 6.0*xi[i];
    
    // Get the section tangent stiffness and stress resultant
    const Matrix &ks = theSections[i]->getInitialTangent();
    
    // Perform numerical integration
    //kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);
    //double wti = wts(i)*oneOverL;
    double wti = wt[i]*oneOverL;
    double tmp;
    int j, k;
    for (j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	for (k = 0; k < order; k++)
	  ka(k,0) += ks(k,j)*wti;
	break;
      case SECTION_RESPONSE_MZ:
	for (k = 0; k < order; k++) {
	  tmp = ks(k,j)*wti;
	  ka(k,1) += (xi6-4.0)*tmp;
	  ka(k,2) += (xi6-2.0)*tmp;
	}
	break;
      default:
	break;
      }
    }
    for (j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_P:
	for (k = 0; k < 3; k++)
	  kb(0,k) += ka(j,k);
	break;
      case SECTION_RESPONSE_MZ:
	for (k = 0; k < 3; k++) {
	  tmp = ka(j,k);
	  kb(1,k) += (xi6-4.0)*tmp;
	  kb(2,k) += (xi6-2.0)*tmp;
	}
	break;
      default:
	break;
      }
    }
    
  }

  return kb;
}



const Matrix&
VaultMacroElement_2d::getInitialStiff()
{
  //const Matrix &kb = this->getInitialBasicStiff();
  //FIXME: IT4R: input of a matrix from python. 
  // Transform to global stiffness
  K = crdTransf->getInitialGlobalStiffMatrix(kb);

  return K;
}




const Matrix&
VaultMacroElement_2d::getMass()
{
  K.Zero();

  if (rho == 0.0)
    return K;
  
  double L = crdTransf->getInitialLength();
  if (cMass == 0)  {
    // lumped mass matrix
    // TODO: Do I need to add also mass moments of inertia for rotations? 
    double m = 0.5*rho*L; 
    K(0,0) = K(1,1) = K(3,3) = K(4,4) = m;

  }
  return K;

}


void
VaultMacroElement_2d::zeroLoad(void)
{
  Q.Zero();

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  
  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  
  return;
}



int 
VaultMacroElement_2d::addInertiaLoadToUnbalance(const Vector &accel)  {
  if (rho == 0.0) {
    return 0;
  }

  // Get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  
   if (3 != Raccel1.Size() || 3 != Raccel2.Size()) {
    opserr << "VaultMacroElement_2d::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
    return -1;
  }

//TODO: Proveri jer si ovo uzela iz DispBasedBeam2d; nemas u rotacijama nista ovako. Da li treba isto i u 3d? 
  if (cMass == 0)  {
     // take advantage of lumped mass matrix
    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;
    
    Q(0) -= m*Raccel1(0);
    Q(1) -= m*Raccel1(1);
    Q(3) -= m*Raccel2(0);
    Q(4) -= m*Raccel2(1);
    

  }

   return 0;
}


//TODO: Here, probably I need to return my reaction force from py.
const Vector&
VaultMacroElement_2d::getResistingForce()
{

  //FIXME: IT4R: Here return reactions from python.
  /* double L = crdTransf->getInitialLength();
  
  double oneOverL = 1.0/L;
  
  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  //const Vector &wts = quadRule.getIntegrPointWeights(numSections);  
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  // Zero for integration
  q.Zero();
  
  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {
    
    int order = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();
  
    //double xi6 = 6.0*pts(i,0);
    double xi6 = 6.0*xi[i];
    
    // Get section stress resultant
    const Vector &s = theSections[i]->getStressResultant();
    
    // Perform numerical integration on internal force
    //q.addMatrixTransposeVector(1.0, *B, s, wts(i));
    
    double si;
    for (int j = 0; j < order; j++) {
      //si = s(j)*wts(i);
      si = s(j)*wt[i];
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	q(0) += si; break;
      case SECTION_RESPONSE_MZ:
	q(1) += (xi6-4.0)*si; q(2) += (xi6-2.0)*si; break;
      default:
	break;
      }
    }    
  }
  
  // Add effects of element loads, q = q(v) + q0 //TODO: this can be set as zero?  
  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
 */

  // Vector for reactions in basic system
  Vector p0Vec(p0, 3);

  P = crdTransf->getGlobalResistingForce(q, p0Vec);  //TODO: Here I should return my reaction from .py? 

  // Subtract other external nodal loads ... P_res = P_int - P_ext
  if (rho != 0)
    P.addVector(1.0, Q, -1.0);
  
  return P;
}


const Vector&
VaultMacroElement_2d::getResistingForceIncInertia()
{
  P = this->getResistingForce();
  
  if (rho != 0.0) {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();
    
  if (cMass == 0)  {
    // take advantage of lumped mass matrix
    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;
    
    P(0) += m*accel1(0);
    P(1) += m*accel1(1);
    P(3) += m*accel2(0);
    P(4) += m*accel2(1);
  } 
    
    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      P.addVector(1.0, this->getRayleighDampingForces(), 1.0);

  } else {
    
    // add the damping forces if rayleigh damping
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
  }

  return P;
}



// Missing methods: Print; displaySelf; setResponse; getResponse



//Copied from ME, just arguments set as zero. 
int
VaultMacroElement_2d::setRayleighDampingFactors(double alpham, double betak, double betak0, double betakc)
{

  alpham = 0.0;
	betak = 0.0;
	betak0 = 0.0;
	betakc = 0.0;

	alphaM = alpham;
	betaK = betak;
	betaK0 = betak0;
	betaKc = betakc;

	// check that memory has been allocated to store compute/return
	// damping matrix & residual force calculations
	if (index == -1) {
		int numDOF = this->getNumDOF();

		for (int i = 0; i<numMatrices; i++) {
			Matrix *aMatrix = theMatrices[i];
			if (aMatrix->noRows() == numDOF) {
				index = i;
				i = numMatrices;
			}
		}
		if (index == -1) {
			Matrix **nextMatrices = new Matrix *[numMatrices + 1];
			if (nextMatrices == 0) {
				opserr << "Element::getTheMatrix - out of memory\n";
			}
			int j;
			for (j = 0; j<numMatrices; j++)
				nextMatrices[j] = theMatrices[j];
			Matrix *theMatrix = new Matrix(numDOF, numDOF);
			if (theMatrix == 0) {
				opserr << "Element::getTheMatrix - out of memory\n";
				exit(-1);
			}
			nextMatrices[numMatrices] = theMatrix;

			Vector **nextVectors1 = new Vector *[numMatrices + 1];
			Vector **nextVectors2 = new Vector *[numMatrices + 1];
			if (nextVectors1 == 0 || nextVectors2 == 0) {
				opserr << "Element::getTheVector - out of memory\n";
				exit(-1);
			}

			for (j = 0; j<numMatrices; j++) {
				nextVectors1[j] = theVectors1[j];
				nextVectors2[j] = theVectors2[j];
			}

			Vector *theVector1 = new Vector(numDOF);
			Vector *theVector2 = new Vector(numDOF);
			if (theVector1 == 0 || theVector2 == 0) {
				opserr << "Element::getTheVector - out of memory\n";
				exit(-1);
			}

			nextVectors1[numMatrices] = theVector1;
			nextVectors2[numMatrices] = theVector2;

			if (numMatrices != 0) {
				delete[] theMatrices;
				delete[] theVectors1;
				delete[] theVectors2;
			}
			index = numMatrices;
			numMatrices++;
			theMatrices = nextMatrices;
			theVectors1 = nextVectors1;
			theVectors2 = nextVectors2;
		}
	}

	// if need storage for Kc go get it
	if (betaKc != 0.0) {
		if (Kc == 0)
			Kc = new Matrix(this->getTangentStiff());
		if (Kc == 0) {
			opserr << "WARNING - ELEMENT::setRayleighDampingFactors - out of memory\n";
			betaKc = 0.0;
		}

		// if don't need storage for Kc & have allocated some for it, free the memory
	}
	else if (Kc != 0) {
		delete Kc;
		Kc = 0;
	}

	return 0;
}

//TODO: DA li i ovde treba staviti parametre da su jednaki nuli? 
const Matrix &
VaultMacroElement_2d::getDamp(void)
{
  alphaM = 0.0;
	betaK = 0.0;
	betaK0 = 0.0;
	betaKc = 0.0;

	if (index == -1) {
		this->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
	}

	// now compute the damping matrix
	Matrix *theMatrix = theMatrices[index];
	theMatrix->Zero();
	if (alphaM != 0.0)
		theMatrix->addMatrix(0.0, this->getMass(), alphaM);
	if (betaK != 0.0)
		theMatrix->addMatrix(1.0, this->getTangentStiff(), betaK);
	if (betaK0 != 0.0)
		if (dampingModel == 0) {
			theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0);
		}
		else {
			if (ductilityDemandCommitted < 1) {
				theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0);
			}
			else {
				if (dampingModel == 1) {
					theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0 / ductilityDemandCommitted);
				}
				else {
					theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0 / ductilityDemandCycle);
				}
			}
		}
	if (betaKc != 0.0)
		theMatrix->addMatrix(1.0, this->getSecantStiff(), betaKc);

	// return the computed matrix
	return *theMatrix;
}


const Vector &
VaultMacroElement_2d::getRayleighDampingForces(void)
{
  
  alphaM = 0.0;
	betaK = 0.0;
	betaK0 = 0.0;
	betaKc = 0.0;

	if (index == -1) {
		this->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
	}

	Matrix *theMatrix = theMatrices[index];
	Vector *theVector = theVectors2[index];
	Vector *theVector2 = theVectors1[index];

	//
	// perform: R = (alphaM * M + betaK0 * K0 + betaK * K) * v
	//            = D * v
	//

	// determine the vel vector from ele nodes
	Node **theNodes = this->getNodePtrs();
	int numNodes = this->getNumExternalNodes();
	int loc = 0;
	for (int i = 0; i<numNodes; i++) {
		const Vector &vel = theNodes[i]->getTrialVel();
		for (int i = 0; i<vel.Size(); i++) {
			(*theVector2)(loc++) = vel[i];
		}
	}

	// now compute the damping matrix
	theMatrix->Zero();
	if (alphaM != 0.0)
		theMatrix->addMatrix(0.0, this->getMass(), alphaM);
	if (betaK != 0.0)
		theMatrix->addMatrix(1.0, this->getTangentStiff(), betaK);
	if (betaK0 != 0.0)
		if (dampingModel == 0) {
			theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0);
		}
		else {
			if (ductilityDemandCommitted < 1.0) {
				theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0);
			}
			else {
				if (dampingModel == 1) {
					theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0 / ductilityDemandCommitted);
				}
				else {
					theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0 / ductilityDemandCycle);
				}
			}
		}
	if (betaKc != 0.0)
		theMatrix->addMatrix(1.0, this->getSecantStiff(), betaKc);

	// finally the D * v
	theVector->addMatrixVector(0.0, *theMatrix, *theVector2, 1.0);

	return *theVector;
}


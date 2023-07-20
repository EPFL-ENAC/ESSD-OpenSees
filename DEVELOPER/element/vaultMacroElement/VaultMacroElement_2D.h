#ifndef VaultMacroelement2d_h // Proveri da li treba da bude istog naziva kao i file? Ne treba, tu definises macro. 
#define VaultMacroelement2d_h

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <BeamIntegration.h>
#include <vector>
#include <string>

class Node;
class SectionForceDeformation;
class CrdTransf; // Da li ovo treba CrdTransf, ili moja matrica za geom transformaciju? Posto imam normalne dof u cvorovima, moze ova default CrdTransf
class Response;
class Channel; // ovo je nesto za output
//class UniaxialMaterial; // This is added to bet he same as Truss2.h

class VaultMacroElement_2d : public Element
{
  public:
    VaultMacroElement_2d(int tag,
                        int nd1, int nd2,
                        int numSections, SectionForceDeformation **s, 
                       // Vector intLength, Vector intLengthMasses,
                        BeamIntegration &bi, CrdTransf &coordTransf,
                        double dt, double theta,
                        int dampingModel = 0,
                        double rho=0.0,
                        int cMass = 0); 
    VaultMacroElement_2d(); //FIXME: Element (as in ENACIT), or VaultMacroElement (as in FV) in the beg.of this line(in green) ?
    ~VaultMacroElement_2d(); //ovo  ~ je destructor

    const char *getClassType(void) const {return "VaultMacroElement_2d";}; // IB: Is it needed? Seems that is for debugging

    // public methods to obtain information about dof & connectivity    
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);
    int getNumDOF(void);
    void setDomain(Domain *theDomain); // initialization

    //  ---- methods dealing with committed state and update - u ME nisu jednake nuli! DA LI TREBA = 0??
    int commitState(void) = 0;  // called when a converged solution has been obtained for a time step
    // public methods to obtain stiffness, mass, damping and residual information
    int update(void) = 0;  // called when a new trial step has been set at the nodes

    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getMass(void);

    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInfo);

    // damping methods
    int setRayleighDampingFactors(double alphaM, double betaK, double betaK0, double betaKc);
    const Matrix &getDamp(void);


    //--- This is from VaultMacroElement: --------------
    // Set input acc
    void setAcc(std::vector<std::vector<double>> acc);
    // Perform optimization
    void compute();
    // Get output for usage outside of the element
    std::vector<std::vector<double>> getResult();
    //---------------------------------------------------

  protected: 
    const Vector &getRayleighDampingForces(void); // set as zero. Here or in cpp? Probably cpp

  private:

    const int numSections;                      // number of sections = 1
    SectionForceDeformation **theSections;      // pointers to sectional models
    BeamIntegration* beamInt;
    CrdTransf* crdTransf;  

    ID  connectedExternalNodes;                 // tags of element nodes (i, j)
    Node *theNodes[2];                          // pointers to element nodes (i,j) 

    static Matrix K;                            // stores element stiffness, damping, or mass Matrix
    static Vector P;                            // Element resisting force vector - Treba za getResistingForce - u getResponse
    
    Vector Q;                                   // Applied nodal loads
    Vector q;                                   //  Basic force
    double q0[3];                                // Fixed end forces in basic system
    double p0[3];                                // Reactions in basic system of distributed loads 

    double rho;                                 // Mass density per unit length
    
    // DO I need this? I think yes, for connectedExternalNode
    int parameterID;   //ovo je isto nesto za sensitivity analysis
    static double workArea[];  // not used here

    //TODO: da li mi ovo treba? sto sam komentarisala?
    //double* nodeIInitialDisp;                    // initial displacements of the end nodes 
    //double* nodeJInitialDisp; 
 	  //Vector intLength, intLengthMasses;           // integration lengths for sectional responses and lumped masses
    double dt;
    double theta;
	
    
    double rho;	   // Mass density per unit length
    int cMass;     // consistent mass flag

    int dampingModel;


    //--------This from VaultMacroELement: --------------
    // static const std::string acc_file_name;
    // static const std::string result_file_name;
    // static const std::string python_script_name;

    // std::vector<std::vector<double>> result;

    // Load output in ::result
    void loadResult();
    //---------------------------------------------------    
	

	// Link with Python
	int initializePython();
	PyObject* pModule;
	PyObject* pUpdateFunc;
	PyObject* pGetResistingForceFunc;
	PyObject* pGetInitialStiffFunc;

};

#endif

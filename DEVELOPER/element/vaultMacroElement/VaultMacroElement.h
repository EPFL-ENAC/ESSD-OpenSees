#ifndef VaultMacroElement_h
#define VaultMacroElement_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <NDMaterial.h>
#include <UniaxialMaterial.h>
#include <vector>
#include <string>


class VaultMacroElement : public Element
{
public:
	VaultMacroElement(double dt, double theta):
		Element(0, 0), dt(dt), theta(theta) {};
	~VaultMacroElement();

	// Set input acc
	void setAcc(std::vector<std::vector<double>> acc);
	// Perform optimization
	void compute();
	// Get output for usage outside of the element
	std::vector<std::vector<double>> getResult();


private:
	static const std::string acc_file_name;
	static const std::string result_file_name;
	static const std::string python_script_name;

	double dt;
	double theta;
	std::vector<std::vector<double>> result;

	// Load output in ::result
	void loadResult();


	// Mandatory methods for Element ------------------------------------------

public:
	// initialization
	void setDomain(Domain *theDomain);

	// methods dealing with nodes and number of external dof
	int getNumExternalNodes(void) const = 0;
	const ID &getExternalNodes(void) = 0;
	Node **getNodePtrs(void) = 0;
	int getNumDOF(void) = 0;

	// methods dealing with committed state and update
	int commitState(void) = 0;  // called when a converged solution has been obtained for a time step
	int revertToLastCommit(void) = 0; // called when the soln algorithm has failed to converge to a solution at a time step
	int revertToStart(void) = 0; // called when model is rest to initial conditions
	int update(void) = 0;  // called when a new trial step has been set at the nodes

	// methods dealing with element stiffness
	const Matrix &getTangentStiff(void) = 0;
	const Matrix &getInitialStiff(void) = 0;

	// methods dealing with element forces 
	void zeroLoad(void) = 0;
	int addLoad(ElementalLoad *theLoad, double loadFactor) = 0;
	const Vector &getResistingForce(void) = 0;

	// public methods for output
	void Print(OPS_Stream &s, int flag = 0) = 0;
	Response *setResponse(const char **argv, int argc, OPS_Stream &theHandler) = 0;
	int getResponse(int responseID, Information &eleInformation) = 0;

	// method for database/parallel processing
	int sendSelf(int commitTag, Channel &theChannel) = 0;
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) = 0;

	// ------------------------------------------------------------------------
};


#endif // VaultMacroElement_h


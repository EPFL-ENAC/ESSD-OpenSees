#include "VaultMacroElement.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>


const std::string VaultMacroElement::acc_file_name = "acc.csv";
const std::string VaultMacroElement::result_file_name = "result.csv";
const std::string VaultMacroElement::python_script_name = "VaultMacroElement.py";


// Export acc into a csv
void VaultMacroElement::setAcc(std::vector<std::vector<double>> acc) {
	std::ofstream acc_file(acc_file_name);

	if (!acc_file.is_open()) {
		std::cerr << "Could not open " << acc_file_name << std::endl;
		return;
	}

	for (int i = 0; i < acc.size(); i++) {
		for (int j = 0; j < acc[i].size(); j++) {
			acc_file << acc[i][j] << ",";
		}
		acc_file << std::endl;
	}

	acc_file.close();
}


// Run python optimization
void VaultMacroElement::compute() {
	// Call python script
	std::string command = "python " + python_script_name + " " + std::to_string(dt) + " " + std::to_string(theta);

    int result = std::system(command.c_str());

    if (result != 0) {
        std::cerr << "Python command failed with error code " << result << std::endl;
    }

	loadResult();
}


// Load optimization result (csv) into memory
void VaultMacroElement::loadResult() {
	std::ifstream result_file(result_file_name);

	if (!result_file.is_open()) {
		std::cerr << "Could not open " << result_file_name << std::endl;
		return;
	}

	result = std::vector<std::vector<double>>(); // clear result currently in memory
	std::string line;

	while (std::getline(result_file, line)) {
		std::stringstream line_stream(line);
		std::vector<double> line_values;
		std::string value;

		while (std::getline(line_stream, value, ',')) {
			line_values.push_back(std::stod(value));
		}

		result.push_back(line_values);
	}

	result_file.close();
}


// Get result in memory
std::vector<std::vector<double>> VaultMacroElement::getResult() {
	return result;
}

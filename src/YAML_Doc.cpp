#include <ctime>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#ifdef REDSTORM
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif
#include "YAML_Doc.hpp"

#if MINIFE_INFO != 0
#include <miniFE_info.hpp>
#else
#include <miniFE_no_info.hpp>
#endif


using namespace std;


//set the microapp_name and version which will become part of the YAML
YAML_Doc::YAML_Doc(const std::string& miniApp_Name,
                   const std::string& miniApp_Version,
                   const std::string& destination_Directory,
                   const std::string& destination_FileName)
{
	miniAppName = miniApp_Name;
	miniAppVersion = miniApp_Version;
	destinationDirectory = destination_Directory;
	destinationFileName = destination_FileName;
}

//inherits the destructor from YAML_Element
YAML_Doc::~YAML_Doc(void)
{}

/*
 * generates YAML from the elements of the document and saves it
 * to a file
 */
string YAML_Doc::generateYAML()
{
	string yaml;
	yaml =  yaml + "Mini-Application Name: " + miniAppName + "\n";
	yaml =  yaml + "Mini-Application Version: " + miniAppVersion + "\n";
	for(size_t i=0; i<children.size(); i++){
		yaml = yaml + children[i]->printYAML("");
	}

	time_t rawtime;
	tm * ptm;
	time ( &rawtime );
	ptm = localtime(&rawtime);
	char sdate[25];
	//use tm_mon+1 because tm_mon is 0 .. 11 instead of 1 .. 12
	sprintf (sdate,"%04d:%02d:%02d-%02d:%02d:%02d",ptm->tm_year + 1900, ptm->tm_mon+1,
	         ptm->tm_mday, ptm->tm_hour, ptm->tm_min,ptm->tm_sec);

	string filename;
	if (destinationFileName=="")
		filename = miniAppName + "-" + miniAppVersion + "_";
	else
		filename = destinationFileName;
	filename = filename + string(sdate) + ".yaml";
	if (destinationDirectory!="" && destinationDirectory!=".") {
		string mkdir_cmd = "mkdir " + destinationDirectory;
		#ifdef REDSTORM
		mkdir(destinationDirectory.c_str(),0755);
		#else
		system(mkdir_cmd.c_str());
		#endif
		filename = destinationDirectory + "/" + destinationFileName;
	}
	else
		filename = "./" + filename;

	ofstream myfile;
	myfile.open(filename.c_str());
	myfile << yaml;
	myfile.close();
	return yaml;
}

void YAML_Doc::add_params_to_yaml(miniFE::Parameters& params)
{
	add("Global Run Parameters","");
	get("Global Run Parameters")->add("dimensions","");
	get("Global Run Parameters")->get("dimensions")->add("nx", params.nx);
	get("Global Run Parameters")->get("dimensions")->add("ny", params.ny);
	get("Global Run Parameters")->get("dimensions")->add("nz", params.nz);
	get("Global Run Parameters")->add("load_imbalance", params.load_imbalance);
	get("Global Run Parameters")->add("mv_overlap_comm_comp", "0 (no)");
}

void YAML_Doc::add_configuration_to_yaml(int numprocs)
{
	get("Global Run Parameters")->add("number of processors", numprocs);

	add("Platform","");
	get("Platform")->add("hostname", MINIFE_HOSTNAME);
	get("Platform")->add("kernel name", MINIFE_KERNEL_NAME);
	get("Platform")->add("kernel release", MINIFE_KERNEL_RELEASE);
	get("Platform")->add("processor", MINIFE_PROCESSOR);

	add("Build","");
	get("Build")->add("CXX", MINIFE_CXX);
	#if MINIFE_INFO != 0
	get("Build")->add("compiler version",MINIFE_CXX_VERSION);
	#endif
	get("Build")->add("CXXFLAGS", MINIFE_CXXFLAGS);
}

void YAML_Doc::add_timestring_to_yaml()
{
	std::time_t rawtime;
	struct tm * timeinfo;
	std::time(&rawtime);
	timeinfo = std::localtime(&rawtime);
	std::ostringstream osstr;
	osstr.fill('0');
	osstr << timeinfo->tm_year+1900 << "-";
	osstr.width(2);
	osstr << timeinfo->tm_mon + 1 << "-";
	osstr.width(2);
	osstr << timeinfo->tm_mday << ", ";
	osstr.width(2);
	osstr << timeinfo->tm_hour << "-";
	osstr.width(2);
	osstr << timeinfo->tm_min << "-";
	osstr.width(2);
	osstr << timeinfo->tm_sec;
	std::string timestring = osstr.str();
	add("Run Date/Time",timestring);
}


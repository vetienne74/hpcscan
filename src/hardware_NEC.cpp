
//-------------------------------------------------------------------------------------------------------
// This class handles all queries related to the hardware
// Specialized class that targets NEC SX-Aurora
// Derived class from Hardware
// Associated to all test modes when hpcscan is compiled with NEC compiler
//-------------------------------------------------------------------------------------------------------

#include "hardware_NEC.h"

extern "C"
{
#include <libsysve.h> // to get NEC device type
}

#include "config.h"
#include "constant.h"
#include "global.h"
#include "output_report.h"

using namespace std;

namespace hpcscan {

typedef struct {
	int type;
	string name;
} s_device;

static s_device device_table[] = {
		{18, "10C-A"},
		{131, "10B-P"},
		{163, "10B-L"},
		{165, "10A-L"},
		{167, "10AE-L"},
		{135, "10BE-P"},
		{147, "10BE-A"},
		{16, "10CE-A"},
		{235, "20A-L"},
		{233, "20B-L"},
		{203, "20B-P"},
		{218, "20B-A"},
		{133, "20B-P(*)"}, // Internal
		{136, "10B-P(*)"} // Internal
};

static int device_table_size = sizeof(device_table) / sizeof(s_device);

//-------------------------------------------------------------------------------------------------------

Hardware_NEC::Hardware_NEC(string gridMode) : Hardware(gridMode)
{
	printDebug(LIGHT_DEBUG, "IN Hardware_NEC::Hardware_NEC");

	supportGetPowerUsage = checkSupportGetPowerUsage() ;

	printDebug(LIGHT_DEBUG, "OUT Hardware_NEC::Hardware_NEC");
}

//-------------------------------------------------------------------------------------------------------

Hardware_NEC::~Hardware_NEC(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware_NEC::~Hardware_NEC");

	// TODO ~Hardware_NEC

	printDebug(LIGHT_DEBUG, "OUT Hardware_NEC::~Hardware_NEC");
}

//-------------------------------------------------------------------------------------------------------

void Hardware_NEC::info(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware_NEC::info");

	print_blank() ;
	printInfo(MASTER, " Hardware information") ;

	printInfo(MASTER, " Target hardware:") ;

	// type of Vector Engine
	std::string s = "type";
	std::string dev;
	char buffer[16];
	int t, i, v;

	ve_get_ve_info(const_cast<char*>(s.c_str()), buffer, 16);
	t = atoi(buffer);
	v = -1;
	for (i = 0; i < device_table_size; i++) {
		if (t == device_table[i].type) {v = i; break;}
	}
	if (v == -1) {
		dev = "NEC SX-Aurora TSUBASA Unknown (" + std::to_string(t) + ")";
	} else {
		dev = "NEC SX-Aurora TSUBASA " + device_table[v].name;
	}

	printInfo(MASTER, " NEC Vector Engine in native mode") ;
	printInfo(MASTER, " Device type:\t", dev );


	// support for power usage
	if (supportGetPowerUsage)
	{
		printInfo(MASTER, " Read power usage", "SUPPORTED") ;
	}
	else
	{
		printInfo(MASTER, " Read power usage", "NOT SUPPORTED") ;
	}
	print_line2() ;

	printDebug(LIGHT_DEBUG, "OUT Hardware_NEC::info");
}

//-------------------------------------------------------------------------------------------------------

bool Hardware_NEC::checkSupportGetPowerUsage(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware_NEC::checkSupportGetPowerUsage");

	bool retVal = false ;

	char *str_ve_node_number;
	char dev_veslot[32];
	char dev_ve[32];

	// TODO adapt to multiple VEs
	int venode = 0 ;

	if (venode == -1) {
		str_ve_node_number = getenv("VE_NODE_NUMBER");
		if (!str_ve_node_number) {
			venode = 0;
		} else {
			venode = atoi(str_ve_node_number);
		}
	}

	sprintf(dev_veslot, "/dev/veslot%d", venode);
	if (readlink(dev_veslot, dev_ve, 32) == -1) {
		fprintf(stderr, "readlink failed\n");
	}
	else
	{
		deviceId = dev_ve[2] - '0';
		retVal = true ;
	}

	printDebug(LIGHT_DEBUG, "OUT Hardware_NEC::checkSupportGetPowerUsage");
	return(retVal) ;
}

//-------------------------------------------------------------------------------------------------------

Myfloat Hardware_NEC::measureCurrentPower(void)
{
	printDebug(LIGHT_DEBUG, "IN Hardware_NEC::measureCurrentPower");

	Myfloat retVal = UNSPECIFIED ;
	if (supportGetPowerUsage)
	{
		FILE *stream;
		char filename[128];
		int val = 0;
		float aux_12v_v, aux_12v_c;
		float edge_12v_v, edge_12v_c;
		float p;

		/* power consumption  = (AUX 12V V) * (AUX 12V C) + (Edge 12V V) * (Edge 12V C) + 5W
		 * AUX  12V V: sensor_8   micro(10^-6) Volt
		 * AUX  12V C: sensor_12  milli(10^-3) Ampere
		 * Edge 12V V: sensor_9   micro(10^-6) Volt
		 * Edge 12V C: sensor_13  milli(10^-3) Ampere
		 */

		/* AUX 12V V, mV */
		sprintf(filename, "/sys/class/ve/ve%d/sensor_8", deviceId);
		stream = fopen(filename, "r");
		fscanf(stream, "%d", &val);
		fclose(stream);
		aux_12v_v = ((float) val) / 1000000.0;

		/* AUX 12V C, mA */
		sprintf(filename, "/sys/class/ve/ve%d/sensor_12", deviceId);
		stream = fopen(filename, "r");
		fscanf(stream, "%d", &val);
		fclose(stream);
		aux_12v_c = (float) val;

		/* Edge 12V V, mV */
		sprintf(filename, "/sys/class/ve/ve%d/sensor_9", deviceId);
		stream = fopen(filename, "r");
		fscanf(stream, "%d", &val);
		fclose(stream);
		edge_12v_v = ((float) val) / 1000000.0;

		/* Edge 12V C, mA */
		sprintf(filename, "/sys/class/ve/ve%d/sensor_13", deviceId);
		stream = fopen(filename, "r");
		fscanf(stream, "%d", &val);
		fclose(stream);
		edge_12v_c = (float) val;

		retVal = (aux_12v_v * aux_12v_c + edge_12v_v * edge_12v_c) / 1000.0 + 5.0;
	}

	printDebug(LIGHT_DEBUG, "OUT Hardware_NEC::measureCurrentPower");
	return(retVal) ;
}


} // namespace hpcscan

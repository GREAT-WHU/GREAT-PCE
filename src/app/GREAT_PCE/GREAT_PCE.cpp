#include <chrono>
#include <thread>

#include "gcfg_pce.h"


using namespace std;
using namespace gnut;
using namespace std::chrono;

void catch_signal(int) { cout << "Program interrupted by Ctrl-C [SIGINT,2]\n"; exit(1);}

// MAIN
// ----------
int main(int argc, char** argv)
{
	// Only to cout the Reminder here
	signal(SIGINT, catch_signal);

	bool thrd = false;

	// Creat and set the log file : clk.log
	t_glog glog;
	glog.mask("great_pcelsq.app_log");
	glog.append(false);
	glog.cache_size(99);
	glog.tsys(t_gtime::GPS);
	glog.time_stamp(true);

	// Construct the gset class and init some values in the class
	// Get the arguments from the command line
	t_gcfg_pce gset;
	gset.app("GREAT/CLK-LSQ", "0.9.0", "$Rev: 2448 $", "(WHU-SGG GREAT)", __DATE__, __TIME__);
	gset.arg(argc, argv, true, false);
	glog.verb(dynamic_cast<t_gsetout*>(&gset)->verb());

	// Prepare site list from gset
	// Prepare input files list form gset
	// Get sample intval from gset. if not, init with the default value
	set<string> sites = dynamic_cast<t_gsetgen*>(&gset)->recs();
	multimap<IFMT, string> inp = gset.inputs_all();
	int sample = int(dynamic_cast<t_gsetgen*>(&gset)->sampling());
	if (!sample) sample = int(dynamic_cast<t_gsetgen*>(&gset)->sampling_default());

	// DECLARATIONS/INITIALIZATIONS
	t_gallobs*    gobs = new t_gallobs();  gobs->glog(&glog); gobs->gset(&gset);
	t_gallnav*    gorb = new t_gallprec(); gorb->glog(&glog);
	dynamic_cast<t_gallprec*>(gorb)->use_clknav(true);
	t_gdata*      gdata = nullptr;
	t_gallpcv*    gpcv = nullptr; if (gset.input_size("atx") > 0) { gpcv = new t_gallpcv;  gpcv->glog(&glog); }
	t_gallotl*    gotl = nullptr; if (gset.input_size("blq") > 0) { gotl = new t_gallotl;  gotl->glog(&glog); }
	t_gallrecover*    grcv = new t_gallrecover(); grcv->glog(&glog);
	t_gallobj*    gobj = new t_gallobj(gpcv, gotl); gobj->glog(&glog);
	t_gallbias*   gbia = nullptr; if (gset.input_size("biasinex") > 0 || gset.input_size("biabern") > 0) {
		gbia = new t_gallbias; gbia->glog(&glog);
	}
	t_gnavde*     gde = new t_gnavde;
	t_gpoleut1*	   gerp = new t_gpoleut1;
	t_gleapsecond*   gleap = new t_gleapsecond;

	//vgclk for clock estimation
	shared_ptr<t_glsqproc> vgclk = nullptr;

	// runepoch for the time costed each epoch 
	t_gtime runepoch(t_gtime::GPS);
	// lstepoch for the time of all epoches 
	t_gtime lstepoch(t_gtime::GPS);

	// SET OBJECTS
	set<string>::const_iterator itOBJ;
	set<string> obj = dynamic_cast<t_gsetrec*>(&gset)->objects();
	for (const auto& name : obj) {
		shared_ptr<t_grec> rec = dynamic_cast<t_gsetrec*>(&gset)->grec(name, &glog);
		gobj->add(rec);
	}

	// CHECK INPUTS, sp3+rinexc+rinexo, Necessary data
	if (gset.input_size("sp3") == 0 && gset.input_size("rinexc") == 0 && gset.input_size("rinexo") == 0) {
		glog.comment(0, "main", "Error: incomplete input: rinexo + rinexc + sp3");
		gset.usage();
	}

	// DATA READING
	// for multiply thread
	vector<t_gcoder*> gcoder_thrd;
	vector<thread> gthread;
	clock_t beg_time = clock();
	clock_t end_time = clock();

	multimap<IFMT, string>::const_iterator itINP = inp.begin();
	for (size_t i = 0; i < inp.size() && itINP != inp.end(); ++i, ++itINP)
	{
		// Get the file format/path, which will be used in decoder
		IFMT   ifmt(itINP->first);
		string path(itINP->second);
		string id("ID" + int2str(i));

		t_gio*    gio = nullptr;
		t_gcoder* gcoder = nullptr;

		// For different file format, we prepare different data container and decoder for them.
		if (ifmt == SP3_INP) { gdata = gorb; gcoder = new t_sp3(&gset, "", 8172); }
		else if (ifmt == RINEXO_INP) { gdata = gobs; gcoder = new t_rinexo(&gset, "", 4096); }
		else if (ifmt == RINEXC_INP) { gdata = gorb; gcoder = new t_rinexc(&gset, "", 4096); }
		else if (ifmt == RINEXN_INP) { gdata = gorb; gcoder = new t_rinexn(&gset, "", 4096); }
		else if (ifmt == ATX_INP) { gdata = gpcv; gcoder = new t_atx(&gset, "", 4096); }
		else if (ifmt == BLQ_INP) { gdata = gotl; gcoder = new t_blq(&gset, "", 4096); }
		else if (ifmt == BIASINEX_INP) { gdata = gbia; gcoder = new t_biasinex(&gset, "", 20480); }
		else if (ifmt == BIABERN_INP) { gdata = gbia; gcoder = new t_biabernese(&gset, "", 20480); }
		else if (ifmt == DE_INP) { gdata = gde; gcoder = new t_dvpteph405(&gset, "", 4096); }
		else if (ifmt == SINEX_INP) { gcoder = new t_sinex(&gset, "", 20480); }
		else if (ifmt == POLEUT1_INP) { gdata = gerp; gcoder = new t_poleut1(&gset, "", 4096); }
		else if (ifmt == LEAPSECOND_INP) { gdata = gleap; gcoder = new t_leapsecond(&gset, "", 4096); }
		else {
			glog.comment(0, "main", "Error: unrecognized format " + int2str(ifmt));
			gdata = nullptr;
		}

		// Check the file path
		if (path.substr(0, 7) == "file://")
		{
			glog.comment(0, "main", "path is file!");
		}

		// READ DATA FROM FILE
		if (gcoder)
		{
			// get read time 
			beg_time = clock();

			gio = new t_gfile;
			gio->glog(&glog);
			gio->path(path);

			// Put the file into gcoder
			gcoder->clear();
			gcoder->path(path);
			gcoder->glog(&glog);
			// Put the data container into gcoder
			gcoder->add_data(id, gdata);
			gcoder->add_data("OBJ", gobj);
			// Put the gcoder into the gio
			// Note, gcoder contain the gdata and gio contain the gcoder
			gio->coder(gcoder);

			runepoch = t_gtime::current_time(t_gtime::GPS);
			// Read the data from file here
			gio->run_read();
			lstepoch = t_gtime::current_time(t_gtime::GPS);
			// Write the information of reading process to log file
			glog.comment(0, "main", "READ: " + path + " time: "
				+ dbl2str(lstepoch.diff(runepoch)) + " sec");
			if (gio) { delete gio; gio = nullptr; }
			if (gcoder) { delete gcoder; gcoder = nullptr; }

			// get read time 
			end_time = clock();
			cout << "Time for " + path + " is " << (end_time - beg_time) / CLOCKS_PER_SEC << " sec " << endl;
		}
	}

	// set antennas for satllites (must be before PCV assigning)
	// assigning PCV pointers to objects
	t_gtime beg = dynamic_cast<t_gsetgen*>(&gset)->beg();
	t_gtime end = dynamic_cast<t_gsetgen*>(&gset)->end();
	int frequency = dynamic_cast<t_gsetproc*>(&gset)->frequency();
	set<string> system = dynamic_cast<t_gsetgen*>(&gset)->sys();
	IFCB_MODEL ifcb_mode = dynamic_cast<t_gsetproc*>(&gset)->ifcb_model();
	gobj->read_satinfo(beg);
	gobj->sync_pcvs();
	
	// ADD DATA
	t_gallproc* data = new t_gallproc();
	if (!data->Add_Data(gobj) || !data->Add_Data(gobs) || !data->Add_Data(gorb) ||
		!data->Add_Data(gbia) || !data->Add_Data(gotl) || !data->Add_Data(grcv) ||
		!data->Add_Data(gde) || !data->Add_Data(gerp) || !data->Add_Data(gleap))
		return -1;

	// PROCESSING
	vgclk = make_shared<t_gpcelsqIF>(&gset, data, &glog);
	
	vgclk->add_coder(gcoder_thrd);
	t_gtime epo(t_gtime::GPS);

	// Start write log for the processing 
	glog.comment(0, "main", " processing started [ Satellite Clock Estimation ]");
	glog.comment(0, "main", beg.str_ymdhms("  beg: ") + end.str_ymdhms("  end: "));
	// Mark for the begin time of processing
	runepoch = t_gtime::current_time(t_gtime::GPS);
	// Main processing function
	cout << beg.str_ymdhms("  beg: ") << endl;
	vgclk->ProcessBatch(data, beg, end);

	vgclk->GenerateProduct();
	// Get the time of duration 
	lstepoch = t_gtime::current_time(t_gtime::GPS);
	// Write finished log
	glog.comment(0, "main", " processing finished : duration  " + dbl2str(lstepoch.diff(runepoch)) + " sec");

	// Normal End
	glog.comment(0, "main", " Normal End! ");
	
	// Delete pointer
	for (size_t i = 0; i < gcoder_thrd.size(); ++i) { delete gcoder_thrd[i]; }; gcoder_thrd.clear();

	if (gobs) { delete gobs; gobs = nullptr; }
	if (gerp) { delete gerp; gerp = nullptr; }
	if (gde) { delete gde; gde = nullptr; }
	if (gpcv) { delete gpcv; gpcv = nullptr; }
	if (grcv) { delete grcv; grcv = nullptr; }
	if (gotl) { delete gotl; gotl = nullptr; }
	if (gleap) { delete gleap; gleap = nullptr; }
	if (gobj) { delete gobj; gobj = nullptr; }
	if (gorb) { delete gorb; gorb = nullptr; }
	if (gbia) { delete gbia; gbia = nullptr; }
	if (data) { delete data; data = nullptr; }
}
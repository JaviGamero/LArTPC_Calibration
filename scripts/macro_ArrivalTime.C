#include <string>

using std::filesystem::directory_iterator;

void macro_ArrivalTime(){

	char dir[256];
	char skip[256] = "/Users/javigamero/MyMac/DS_Master/TFM/data/sample_particles_v2/.DS_Store";

  	getcwd(dir, 256); 
  	strcat(dir, "/data/sample_particles_v2"); // you must be in main folder of the project

  	// PMT IDs in the TPC with X<0
  	int realisticPMT_IDs_TPC0[60] = {6, 8, 10, 12, 14, 16, 36, 38, 40, 60, 62, 64, 66, 68, 70, 84, 86, 88, 90, 92, 94, 114, 116, 118, 138, 140, 142, 144, 146, 148, 162, 164, 166, 168, 170, 172, 
  	192, 194, 196, 216, 218, 220, 222, 224, 226, 240, 242, 244, 246, 248, 250, 270, 272, 274, 294, 296, 298, 300, 302, 304};
	  
  	// PMT IDs in the TPC with X>0
  	int realisticPMT_IDs_TPC1[60] = {7, 9, 11, 13, 15, 17, 37, 39, 41, 61, 63, 65, 67, 69, 71, 85, 87, 89, 91, 93, 95, 115, 117, 119, 139, 141, 143, 145, 147, 149, 163, 165, 167, 169, 171, 173, 
  	193, 195, 197, 217, 219, 221, 223, 225, 227, 241, 243, 245, 247, 249, 251, 271, 273, 275, 295, 297, 299, 301, 303, 305};
  
	// All PMTs IDs
  	int realisticPMT_IDs[120] = {6, 8, 10, 12, 14, 16, 36, 38, 40, 60, 62, 64, 66, 68, 70, 84, 86, 88, 90, 92, 94, 114, 116, 118, 138, 140, 142, 144, 146, 148, 162, 164, 166, 168, 170, 172, 192, 
  	194, 196, 216, 218, 220, 222, 224, 226, 240, 242, 244, 246, 248, 250, 270, 272, 274, 294, 296, 298, 300, 302, 304, 7, 9, 11, 13, 15, 17, 37, 39, 41, 61, 63, 65, 67, 69, 71, 85, 87, 89, 91, 
  	93, 95, 115, 117, 119, 139, 141, 143, 145, 147, 149, 163, 165, 167, 169, 171, 173, 193, 195, 197, 217, 219, 221, 223, 225, 227, 241, 243, 245, 247, 249, 251, 271, 273, 275, 295, 297, 299, 
  	301, 303, 305};

  	vector<int> selected_pmtid;
  	vector<int> selected_pmtid_2;
  	selected_pmtid_2.assign(realisticPMT_IDs, realisticPMT_IDs + 120);
  	vector<int>::iterator it;
  	vector<int>::iterator it2;
  
  	const double shiftStamp = 135.;//ns
  	const double fSamplingTime = 2.;//ns
  	const int fBaseline=8000;//ADC

  	for (const auto &folder: directory_iterator(dir)) // iterate above folders
  	{  
	  	if (folder.path() == skip) continue; // skip .DS_Store/ folder
	  	for (const auto &treePath: directory_iterator(folder)) // above trees in folder
	  	{

		    // =======================================================================
		    // Read tree with root
		    // =======================================================================
		    char *path = new char[treePath.path().string().length()+1];
		    strcpy(path, treePath.path().string().c_str());

		    cout << path << endl;
		    TFile *fh = new TFile(path);
			TTree *tree = fh -> Get<TTree>("opanatree/OpAnaTree");
		   
		    // tree->Print(); 

		   	// =======================================================================
		    // Taking variables (branches of interest) from the tree
		    // =======================================================================
		    unsigned int event; tree->SetBranchAddress("eventID",&event);

		    vector<vector<double> >* fSimPhotonsLiteVIS=new vector<vector<double> >();
		    vector<vector<double> >* fSimPhotonsLiteVUV=new vector<vector<double> >();
		    tree->SetBranchAddress("SimPhotonsLiteVIS",&fSimPhotonsLiteVIS);
		    tree->SetBranchAddress("SimPhotonsLiteVUV",&fSimPhotonsLiteVUV);
		    
		    vector<vector<double>> * fstepX=new vector<vector<double>>;
		    tree->SetBranchAddress("stepX",&fstepX);

		    // vector<double> * fE=new vector<double>; tree->SetBranchAddress("E",&fE);
		    vector<double> * fdE=new vector<double>; tree->SetBranchAddress("dE",&fdE); // deposited Energy // (MeV)
		    vector<int> * ftrackID=new vector<int>; tree->SetBranchAddress("trackID",&ftrackID);
		    vector<int> * fmotherID=new vector<int>; tree->SetBranchAddress("motherID",&fmotherID);
		    vector<int> * fPDGcode=new vector<int>; tree->SetBranchAddress("PDGcode",&fPDGcode);
		    vector<string> * fprocess=new vector<string>; tree->SetBranchAddress("process",&fprocess);

		    // vector<vector<double> >* fSignalsDigi=new vector<vector<double> >();
		    // tree->SetBranchAddress("SignalsDigi",&fSignalsDigi);
		    // vector<double> * fStampTime=new vector<double>; tree->SetBranchAddress("StampTime", &fStampTime);
		    // vector<int> * fOpChDigi=new vector<int>; tree->SetBranchAddress("OpChDigi",&fOpChDigi);

		    // vector<vector<double> >* fSignalsDeco=new vector<vector<double> >();
		    // tree->SetBranchAddress("SignalsDeco",&fSignalsDeco);
		    // vector<double> * fStampTimeDeco=new vector<double>; tree->SetBranchAddress("StampTimeDeco",&fStampTimeDeco);
		    // vector<int> * fOpChDeco=new vector<int>; tree->SetBranchAddress("OpChDeco",&fOpChDeco);
		    
		    
		    // =======================================================================
		    // Photons arrival time histogram
		    // =======================================================================
		    
		    // TH1D --> 1D histogram with 1 double per channel
		    TH1D* hVUV__total = new TH1D("","",1000,0,10000); // name, title, bins, xdown, xup

		    for(int i=0; i<tree->GetEntries(); i++)
		    {
		      	tree->GetEntry(i);
		      	cout << endl << endl << "Entry: " << i << endl;
		      	cout << "dE size: " << fdE -> size() << endl;
		      
		      	double dE_electron = 0;
		      	double X_start_electron = 0;

		      	//let's see where (in which TPC) the energy is deposited by the electon (if any)
		      	for(int j=0; j<fstepX->size(); j++) 
		      	{
					cout<<"PDGs: "<<fPDGcode->at(j)<<endl;

					if(fPDGcode->at(j)!=-11) continue; // check we have a Michel electron after decay

					dE_electron = fdE->at(j); // (MeV)
					X_start_electron = fstepX->at(j).at(0);	

					cout << "Process: " << fprocess->at(j) <<"	pdg: "
						<< fPDGcode->at(j) << "	motherID: " << fmotherID -> at(j)
						<<"	dE: "<< dE_electron <<"	X_start_electron: "<< X_start_electron 
						<< endl;	
		      	}
		      
		      	if(!(dE_electron >0)) continue;// we need the electron in the active volume
				
				cout<<"Event: "<<event<<endl;
		      
		      	//Selecting in which TPC we will look for the signal
		      	if(X_start_electron <= 0) // takes all the 60 PMT with X<0
					selected_pmtid.assign(realisticPMT_IDs_TPC0, realisticPMT_IDs_TPC0 + 60); 
		      	else // or X>0
					selected_pmtid.assign(realisticPMT_IDs_TPC1, realisticPMT_IDs_TPC1 + 60);

		      	//This loop is for TRUE signal
		      	//Averaging the PMT signals (VUV comnponent) of the event
		      	hVUV__total->Reset();
		      	for(int k=0; k<fSimPhotonsLiteVUV->size(); k++) 
		      	{
						//selecting the PMTs in the array we want to use
						it = find (selected_pmtid.begin(), selected_pmtid.end(), k);
						if (it == selected_pmtid.end()) continue;	

						for(int j=0; j<fSimPhotonsLiteVUV->at(k).size(); j++)
						 	hVUV__total->Fill(fSimPhotonsLiteVUV->at(k).at(j));	
		      	}
		      	//Plotting true averaged signal
		      	TCanvas *can1 = new TCanvas("can1", "Arrival time histogram",200,200,600,500);
		      	can1->cd(1);
		      	hVUV__total->SetTitle(0);
		      	hVUV__total->GetXaxis()->SetTitle("time [ns]");
		      	hVUV__total->GetYaxis()->SetTitle("entries");
		      	hVUV__total->SetLineColor(1);
		      	hVUV__total->Draw();   
		      	can1->Update();
		      	can1->Modified();
		      	can1->WaitPrimitive();

		    }

		    delete tree;
	    	delete fh;
		    break; // to read only the first tree

	    } // over trees in directories

  	} // over directories

  	return;
}

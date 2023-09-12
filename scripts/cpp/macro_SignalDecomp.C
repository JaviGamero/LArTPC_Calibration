// Libraries
#include <string>
#include <iostream>
#include <unistd.h>
#include <filesystem>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TPolyLine3D.h>


// Namespaces
using namespace std::filesystem;
using namespace std;

// MAIN
void macro_SignalDecomp()
{
	char dir[256];
	char skip[256] = "/Users/javigamero/MyMac/DS_Master/TFM/data/sample_particles_v2/.DS_Store";

  	getcwd(dir, 256); 
  	strcat(dir, "/data/sample_particles_v2"); // you must be in main folder of the project

  	for (const auto &folder: directory_iterator(dir)) // iterate above folders
  	{
  		if (folder.path() == skip) continue; // skip .DS_Store/ folder
	  	for (const auto &treePath: directory_iterator(folder)) // above trees in folder
	  	{
	  		// =================================================================
		    // Read tree with root
		    // =================================================================
		    char *path = new char[treePath.path().string().length()+1];
		    strcpy(path, treePath.path().string().c_str());

		    // cout << path << endl;
		    TFile *fh = new TFile(path);		 
			TTree *tree = fh -> Get<TTree>("opanatree/OpAnaPerTrackTree"); // real data

			
			// =================================================================
		    // Read interesting variables 
		    // =================================================================

		    unsigned int run; tree->SetBranchAddress("runID",&run);
		    unsigned int subrun; tree->SetBranchAddress("subrunID",&subrun);
		    unsigned int event; tree->SetBranchAddress("eventID",&event);
		    int trackID; tree->SetBranchAddress("TrackID",&trackID);

		    vector<vector<double> >* fSimPhotonsVIS=new vector<vector<double> >();
		    vector<vector<double> >* fSimPhotonsVUV=new vector<vector<double> >();
		    tree->SetBranchAddress("SimPhotonsVIS",&fSimPhotonsVIS);
		    tree->SetBranchAddress("SimPhotonsVUV",&fSimPhotonsVUV);



			// =================================================================
		    // Signal Decomposition
		    // =================================================================

			TH1D* hVUV_mu = new TH1D("","",1000,0,10000);
		    TH1D* hVUV_e = new TH1D("","",1000,0,10000);
		    TH1D* hVUV_total = new TH1D("","",1000,0,10000);
		    int pre_ID = 1;

		    for (int i=0; i<tree->GetEntries(); i++)
		    {
		    	tree->GetEntry(i);
		    	if(event != pre_ID)
				{
				  	hVUV_total->Add(hVUV_mu);
					  hVUV_total->Add(hVUV_e); 
					  TCanvas *can2 = new TCanvas("can2", "can2",200,200,800,500);
					  can2->cd(1);
					  hVUV_total->SetTitle(0);
					  hVUV_total->GetXaxis()->SetTitle("time [ns]");
					  hVUV_total->GetYaxis()->SetTitle("entries");
					  hVUV_total->SetLineColor(1);
					  hVUV_total->Draw();
					  hVUV_mu->SetLineColor(4);
					  hVUV_mu->Draw("same");
					  hVUV_e->SetLineColor(2);
					  hVUV_e->Draw("same");
					  
					  can2->Update();
					  can2->Modified();
					  can2->WaitPrimitive();
					  
					  hVUV_mu->Reset();	  
					  hVUV_e->Reset();
					  hVUV_total->Reset();
				}
	  
		  		cout<< "it: " << i << ", run:" << run << ", subrun: " << subrun << ", event: " << event << ", trackID:" << trackID <<endl;
		      
			  	for(int k=0; k<fSimPhotonsVUV->size(); k++) 
		  		{
			
					if(trackID == 1)
			  		for(int j=0; j<fSimPhotonsVUV->at(k).size(); j++) 
			  			hVUV_mu->Fill(fSimPhotonsVUV->at(k).at(j));
			  		
					
					else
			  		for(int j=0; j<fSimPhotonsVUV->at(k).size(); j++)
			    		hVUV_e->Fill(fSimPhotonsVUV->at(k).at(j));	
		    	}
	      
	    		pre_ID = event; 
		    	

		    } // loop over entries
		    delete tree;
	    	delete fh;

		} // loop over .roots in a folder
		

  	} // loop over folders

} // macro end


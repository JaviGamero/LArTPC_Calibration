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

using std::filesystem::directory_iterator;

void macro_DigDec(){

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
		    // TFile *fh = new TFile(path);
		    TFile* fh= new TFile("/Users/javigamero/MyMac/DS_Master/TFM/data/sample_particles_v2/60183057_4/opana_tree_wwf_197261bb-f7e5-406d-b34c-4f4dd517c023.root");
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
		    // vector<int> * fmotherID=new vector<int>; tree->SetBranchAddress("motherID",&fmotherID);
		    vector<int> * fPDGcode=new vector<int>; tree->SetBranchAddress("PDGcode",&fPDGcode);
		    vector<string> * fprocess=new vector<string>; tree->SetBranchAddress("process",&fprocess);

		    //Digitised signal variables
		    vector<vector<double> >* fSignalsDigi=new vector<vector<double> >();
		    tree->SetBranchAddress("SignalsDigi",&fSignalsDigi);
		    vector<double> * fStampTime=new vector<double>; tree->SetBranchAddress("StampTime", &fStampTime);
		    vector<int> * fOpChDigi=new vector<int>; tree->SetBranchAddress("OpChDigi",&fOpChDigi);
		    vector<double> x_raw, y_raw;

		    // Deconvolutioned signal variables
		    vector<vector<double> >* fSignalsDeco=new vector<vector<double> >();
		    tree->SetBranchAddress("SignalsDeco",&fSignalsDeco);
		    vector<double> * fStampTimeDeco=new vector<double>; tree->SetBranchAddress("StampTimeDeco",&fStampTimeDeco);
		    vector<int> * fOpChDeco=new vector<int>; tree->SetBranchAddress("OpChDeco",&fOpChDeco);
		    vector<double> x_deco, y_deco;
		    
		    
		    // =========================================================================
	    	// DIGITALISED AND DECONVOLUTIONED SIGNAL
	    	// ========================================================================= 

		    for(int i=0; i<tree->GetEntries(); i++)
		    {
		      	tree->GetEntry(i); //what is an entry exactly? physically
		      
		      	double dE_electron = 0;
      			double X_start_electron = 0;

      			//let's see where (in which TPC) the energy is deposited by the electon (if any)
			    for(int j=0; j<fstepX->size(); j++)  // What is fstepX ???????????????????????????????
			    {
					cout<<"PDGs: "<<fPDGcode->at(j)<<endl;
					if(fPDGcode->at(j)!=-11) continue; // check we have a Michel electron after decay

					dE_electron = fdE->at(j); // (MeV) 
					X_start_electron = fstepX->at(j).at(0);
					
					// to see the steps of the electron
					for (int l = 0; l<fstepX->at(j).size(); l++)
						cout<<fstepX->at(j).at(l)<<endl;

					cout<<fprocess->at(j)<<"  pdg: "<<fPDGcode->at(j)<<"  dE:"<<dE_electron<<" X_start_electron:"<<X_start_electron<<endl;	
      			}
      
		      	if(!(dE_electron >0)) continue;// we need the electron in the active volume
					cout<<"Event: "<<event<<endl;
      
			    //Selecting in which TPC we will look for the signal
			  	if(X_start_electron <= 0) // takes all the 60 PMT with X<0
					selected_pmtid.assign(realisticPMT_IDs_TPC0, realisticPMT_IDs_TPC0 + 60); 
			  	else // or X>0
					selected_pmtid.assign(realisticPMT_IDs_TPC1, realisticPMT_IDs_TPC1 + 60);


	    		// ADC --> Apparent Diffusion Coefficient
	    		// OpChDigi --> associated Photon-Detector ID
	      		//loop over the different PMTs
	      		cout<<"fOpChDigi->size() "<<fOpChDigi->size()<<endl;
		      	for(int k=0; k<fOpChDigi->size(); k++) // for every photon detector
		      	{
					double max_digi=-1;
					double max_deco=-1;
					
					//selecting the PMTs in the array we want to use
					it = find (selected_pmtid.begin(), selected_pmtid.end(), fOpChDigi->at(k));
					if (it == selected_pmtid.end()) continue;
					
					//This loop is for DIGITISED signal
					x_raw.clear(); y_raw.clear();    
					for(int j=0; j<fSignalsDigi->at(k).size(); j++)
					{
			  			double rawADC = fSignalsDigi->at(k).at(j);
			  			rawADC-=fBaseline;

			  			if(max_digi < -rawADC) max_digi = -rawADC; // to limit the plot
			  
			  			double t = fSamplingTime*j+fStampTime->at(k)*1000-shiftStamp;
			  			if(t>-1000 && t<11000)
		  				{
		    				x_raw.push_back(t);
		    				y_raw.push_back(-rawADC);
		  				}	  
					}
			
					///This loop is for DECONVOLVED/RECONSTRUCTED signal
					x_deco.clear(); y_deco.clear();    				
					for(int j=0; j<fSignalsDeco->at(k).size(); j++)
					{
		  				double decoADC = fSignalsDeco->at(k).at(j)/500.; 
		  				double t = fSamplingTime*j+fStampTimeDeco->at(k)*1000-shiftStamp; 
		  
		  				if(max_deco < decoADC) max_deco = decoADC; // to limit the plot
		  
		  				if(t>-1000 && t<11000)
		  				{
		    				x_deco.push_back(t);
		    				y_deco.push_back(decoADC);
		  				}	    
					}	  
		
	        
					const int dim = x_raw.size();
					if(!(dim>0)) continue;

					cout<<"PMT:  "<<fOpChDigi->at(k)<<endl;
					TGraph *g1=new TGraph(dim,&x_raw[0],&y_raw[0]);
		
					const int dim2 = x_deco.size();
					if(!(dim2>0)) continue;
					
					cout<<"PMT:  "<<fOpChDeco->at(k)<<endl;
					
					TGraph *g2=new TGraph(dim2,&x_deco[0],&y_deco[0]);
					TCanvas *can2 = new TCanvas("can2", "can2",200,200,1200,500);
		
					can2->Divide(2,1);
					can2->cd(1);
					g1->SetTitle("Digitised signal");
					g1->GetXaxis()->SetTitle("Time [ns]");
					g1->GetYaxis()->SetTitle("ADC");
					g1->SetMarkerStyle(20);
					g1->GetXaxis()->SetRangeUser(-1000,20000);
					g1->GetYaxis()->SetRangeUser(-100,1.1*max_digi);
					g1->Draw("al");
					can2->cd(2);
					g2->SetTitle("Deconvolutioned signal");
					g2->GetXaxis()->SetTitle("Time [ns]");
					g2->GetYaxis()->SetTitle("ADC");
					g2->SetMarkerStyle(20);
					g2->GetXaxis()->SetRangeUser(-1000,20000);
					g2->GetYaxis()->SetRangeUser(-5,1.1*max_deco);
					g2->Draw("al"); 
					
					can2->Update();
					can2->Modified();
					can2->WaitPrimitive();
	       
					delete g1;
					delete g2;
		      	} 
		    }

		    delete tree;
	    	delete fh;
		    break; // to read only the first tree

	    } // over trees in directories

	  } // over directories

  return;
}

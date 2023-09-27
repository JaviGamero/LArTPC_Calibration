// Libraries
#include <string>
#include <iostream>
#include <fstream>
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

void AddToCSV_T(int v[], int dim, string idx, ofstream &file);

// MAIN
void macro_CreateCSV_decomp()
{
	char dir[256];
	char skip[256] = "/Users/javigamero/MyMac/DS_Master/TFM/data/sample_particles_v2/.DS_Store";

  	getcwd(dir, 256); 
  	strcat(dir, "/data/sample_particles_v2"); // you must be in main folder of the project

	char data_preproc[256];
    getcwd(data_preproc, 256);
    strcat(data_preproc, "/data_preproc/");
	
    ofstream file; // file where all temporary series will be written
    file.open(strcat(data_preproc, "??.csv"));
    string isclosed="False";

	int ntree = 0;

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
			TH1D* hVIS_mu = new TH1D("","",1000,0,10000);
		    TH1D* hVUV_e = new TH1D("","",1000,0,10000);
			TH1D* hVIS_e = new TH1D("","",1000,0,10000);
		    int pre_ID = 1;

		    for (int i=0; i<tree->GetEntries(); i++)
		    {
		    	tree->GetEntry(i);
		    	if(event != pre_ID)
				{
					// Add to csv
					int dim = hVUV_mu->GetNbinsX();
					TAxis *X = hVUV_mu->GetXaxis();

					int timeSerieVUV[dim];
					int timeSerieVIS[dim];
					int idealTimeSerie[dim];
					int x[dim];

					for (int j=0; j<=dim+1; j++) // finishes in Nbins+1, if there are N bins, there will be N+1 points
					{
						timeSerieVUV[j] = hVUV_e -> GetBinContent(j);
						timeSerieVIS[j] = hVIS_e -> GetBinContent(j);
						idealTimeSerie[j] = timeSerieVIS[j] + timeSerieVUV[j];
					}   

					for (int j=0; j<dim; j++)
                    	x[j] = X->GetBinCenter(j+1); // Bin indexes start at 1

					string idx = to_string(ntree) + "_" + to_string(event);
                	AddToCSV_T(timeSerieVIS, dim, idx, file);
					
					hVUV_mu->Reset();	
					hVIS_mu->Reset();	  
					hVUV_e->Reset();
					hVIS_e->Reset();	  
				}
		      
			  	for(int k=0; k<fSimPhotonsVUV->size(); k++) 
		  		{
					if(trackID == 1){
			  			for(int j=0; j<fSimPhotonsVUV->at(k).size(); j++) 
			  				hVUV_mu->Fill(fSimPhotonsVUV->at(k).at(j));

						for(int j=0; j<fSimPhotonsVIS->at(k).size(); j++) 
			  				hVIS_mu->Fill(fSimPhotonsVIS->at(k).at(j));
					}

			  		
					else{
						for(int j=0; j<fSimPhotonsVUV->at(k).size(); j++)
			    			hVUV_e->Fill(fSimPhotonsVUV->at(k).at(j));	

						for(int j=0; j<fSimPhotonsVIS->at(k).size(); j++)
			    			hVIS_e->Fill(fSimPhotonsVIS->at(k).at(j));	
					}
			  		
		    	}
	      
	    		pre_ID = event; 
		    	

		    } // loop over entries
		    delete tree;
	    	delete fh;

			ntree++;
			cout << "Tree: " << ntree << endl;

		} // loop over .roots in a folder
		

  	} // loop over folders

	file.close();

} // macro end


void AddToCSV_T(int v[], int dim, string idx, ofstream &file)
{
    /* This function adds an array to a csv as a line. Hence the csv will have
    "dim" columns. 

    Variables: 
    - v: array to write. ONLY ONE
    - dim: dimension of the array, also the final number of columns in the csv
    - file: ofstream file, it MUST BE OPENED OUTSIDE the function
    */

    file << idx + ";";
    for (int i=0; i<dim-1; i++) 
        file << to_string(v[i]) + ";";

    file << to_string(v[dim-1]) + "\n";

    return;

}
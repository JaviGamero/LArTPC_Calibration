"""
ETL_Process.py
Author: Javier Gamero Muñoz

This file is responsible for extracting data making use of root, for further 
transformation to usual python data types and a final load to a csv. 
This will make the fact of reading data easier. 
"""

# Libraries
import uproot
import matplotlib.pyplot as plt 
import os 
import time 

# MAIN 
ROOT = os.path.join(os.getcwd(), "data/sample_particles_v2/")
skip = os.path.join(ROOT, ".DS_Store"); skip += "/"
t0 = time.time()

for folder in os.listdir(ROOT): 
    PATH = os.path.join(ROOT, folder)
    PATH += "/"
    
    if (PATH == skip): continue
    
    for f in os.listdir(PATH):
        file = os.path.join(PATH, f)
        
        with uproot.open(file) as rootfile: 
            tree = rootfile["opanatree/OpAnaTree"]
            # print(tree.keys())
        
            # branches = tree.arrays() # to active all branches --> inefficient
            branches = tree.arrays(["SimPhotonsLiteVUV", "SimPhotonsLiteVIS"#, 
                                    # "SignalsDeco", "SignalsDigi"#,
                                    # "stepX", "dE",
                                    # "trackID", "motherID", "motherID", "process"
                                    ], library="np")
            
            # print(len(branches["SimPhotonsLiteVUV"]))
            # print(len(branches["SignalsDeco"][0]), len(branches["SignalsDeco"][0]))
            
            # figure, ax = plt.subplots(1,1)
            # ax.hist(list(branches["SimPhotonsLiteVUV"][0]), bins=50)
            # plt.show()
        
print("Execution Time: ", time.time()-t0)


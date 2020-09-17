// Main script
void RecoTreeMaker()
{

  // Set up output tree
  TFile tree_file("wbls3pct_reco.root", "recreate");
  TTree *tree = new TTree("tree", "");

  // Define branch variable
  std::string int_name;
  double mc_x, mc_y, mc_z;
  double reco_x, reco_y, reco_z;
  double mc_time, mc_energy;
  double pmt_hits, inner_hits, veto_hits;
  double reco_pe, reco_time, n_triggers;
  double n9;

  // Define branches
  tree->Branch("int_name", &int_name);  
  tree->Branch("mc_x", &mc_x);  
  tree->Branch("mc_y", &mc_y);  
  tree->Branch("mc_z", &mc_z);  
  tree->Branch("reco_x", &reco_x);  
  tree->Branch("reco_y", &reco_y);  
  tree->Branch("reco_z", &reco_z);  
  tree->Branch("mc_energy", &mc_energy);  
  tree->Branch("mc_time", &mc_time);  
  tree->Branch("pmt_hits", &pmt_hits);  
  tree->Branch("inner_hits", &inner_hits);  
  tree->Branch("veto_hits", &veto_hits);  
  tree->Branch("reco_pe", &reco_pe);  
  tree->Branch("reco_time", &reco_time);  
  tree->Branch("n_triggers", &n_triggers);  
  tree->Branch("n9", &n9);  

  // List of all directories to include
  std::vector<TString> dirs = {
    "IBDPositronHeyshamSig_LIQUID_ibd_p_hs", 
    "IBDPositronHeyshamBkg_LIQUID_ibd_p_hb",
    "208Tl_LIQUID_CHAIN_232Th_NA", 
    "212Bi_LIQUID_CHAIN_232Th_NA", 
    "214Pb_LIQUID_CHAIN_222Rn_NA", 
    "40K_PMT_40K_NA", 
    "208Tl_PMT_CHAIN_232Th_NA", 
    "212Bi_PMT_CHAIN_232Th_NA", 
    "214Pb_PMT_CHAIN_238U_NA", 
    "IBDNeutron_LIQUID_ibd_n", 
    "210Bi_LIQUID_CHAIN_222Rn_NA", 
    "212Pb_LIQUID_CHAIN_232Th_NA", 
    "228Ac_LIQUID_CHAIN_232Th_NA", 
    "IBDPositron_LIQUID_ibd_p", 
    "210Bi_PMT_CHAIN_238U_NA", 
    "212Pb_PMT_CHAIN_232Th_NA", 
    "228Ac_PMT_CHAIN_232Th_NA", 
    "IBD_LIQUID_pn_ibd", 
    "210Tl_LIQUID_CHAIN_222Rn_NA", 
    "214Bi_LIQUID_CHAIN_222Rn_NA", 
    "234Pa_PMT_CHAIN_238U_NA", 
    "li9_LIQUID_A_Z",
    "210Tl_PMT_CHAIN_238U_NA", 
    "214Bi_PMT_CHAIN_238U_NA", 
    "40K_LIQUID_40K_NA",
    "n17_LIQUID_A_Z"};

  // Loop over the input directories
  for(int f = 0; f < dirs.size(); f++){

    int_name = dirs[f];
    // Parent directory
    TString directory = "/p/gpfs1/adg/path_a/wbls_3pct/production_06232020/bonsai_root_files_detectorMedia_wbls_3pct_WM_0420_lightSim/";
    TString dirname = directory+"Watchman_"+dirs[f]+"/";
    // FIXME If any of the files are in a different parent directory...
    if(f <= 1){
      directory = "/usr/WS2/brooks50/temp/bonsai_root_files/";
      dirname = directory+"Watchman+"+dirs[f]+"/";
    }
    
    std::cout<<dirname<<"\n";

    // Get all of the files in the directory
    TString fname; 
    TSystemDirectory dir(dirname, dirname); 
    TList *files = dir.GetListOfFiles(); 
    if (!files) continue;

    // Loop over the files
    TSystemFile *file; 
    TIter next(files); 
    while ((file=(TSystemFile*)next())) { 
      fname = file->GetName(); 
      std::cout<<fname<<"\n";
      // Check it's a root file
      if (file->IsDirectory() || !fname.EndsWith(".root")) continue;
      // FIXME skip any broken files
      if (fname=="run1594184246894525532.root") continue;

      // Get PMT info from the file
      TFile b_file(dirname+fname); //path to your rat-pac files

      // Set up a TTreeReader to access the tree information
      TTreeReader reader("data", &b_file);

      TTreeReaderValue<double> x(reader, "x");
      TTreeReaderValue<double> y(reader, "y");
      TTreeReaderValue<double> z(reader, "z");
      TTreeReaderValue<double> t(reader, "t");
      TTreeReaderValue<double> mcx(reader, "mcx");
      TTreeReaderValue<double> mcy(reader, "mcy");
      TTreeReaderValue<double> mcz(reader, "mcz");
      TTreeReaderValue<double> mcenergy(reader, "mc_energy");
      TTreeReaderValue<double> mct(reader, "mct");
      TTreeReaderValue<double> reco_n9(reader, "n9");
      TTreeReaderValue<double> pe(reader, "pe");
      TTreeReaderValue<int> subid(reader, "subid");
      TTreeReaderValue<int> inner_hit(reader, "inner_hit");
      TTreeReaderValue<int> veto_hit(reader, "veto_hit");

      // Loop through the events
      while(reader.Next()){

        // Only look at the first triggered event (skip retriggers)
        if(*subid != 0) continue;
        // Fill the tree
        mc_x = *mcx; 
        mc_y = *mcy;
        mc_z = *mcz;
        reco_x = *x; 
        reco_y = *y;
        reco_z = *z;
        mc_time = *mct;
        mc_energy = *mcenergy;
        n9 = *reco_n9;
        pmt_hits = *inner_hit + *veto_hit;
        inner_hits = *inner_hit;
        veto_hits = *veto_hit;
        reco_pe = *pe;
        reco_time = *t;

        tree->Fill();

      }

    }

  }

  // Write the tree to file
  tree_file.cd();
  tree->Write();
  tree_file.Close();
}

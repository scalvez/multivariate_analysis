1#!/bin/bash

(
    for isotope_str in 0nu 2nu tl208 bi214
    do
        echo "Isotope : $isotope_str"
        if $isotope_str -eq 2nu; then
            sed -i -e 's@.*Long64.*@for (Long64_t ievt=0; ievt<20000;ievt++) {@g' classification_application.C
        fi


        sed -i -e 's@.*TString fname.*@TString fname = "./root_export_'$isotope_str'_25G.root";@g' classification_application.C
        sed -i -e 's@.*bdt_score.*@TFile *target = new TFile ("bdt_score_'$isotope_str'.root","RECREATE");@g' classification_application.C

        root -l -q classification_application.C
        if $isotope_str -eq 2nu; then
            sed -i -e 's@.*Long64.*@for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {@g' classification_application.C
        fi

    done
)

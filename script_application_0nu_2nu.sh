#!/bin/bash

(
    for isotope_str in 0nu 2nu
    do
        echo "Isotope : $isotope_str"

        isotope=$isotope_str
        # bb='2nu'
        # if [ "$isotope" = "$bb" ]; then
        #     sed -i -e 's@.*Long64.*@for (Long64_t ievt=20000; ievt<40000;ievt++) {@g' classification_application_simple.C
        # fi

        sed -i -e 's@.*Long64.*@for (Long64_t ievt=0; ievt<20000;ievt++) {@g' classification_application_simple.C

        sed -i -e 's@.*TString fname.*@TString fname = "./root_export_'$isotope_str'_25G.root";@g' classification_application_simple.C
        sed -i -e 's@.*bdt_score.*@TFile *target = new TFile ("bdt_score_'$isotope_str'.root","RECREATE");@g' classification_application_simple.C

        root -l -q classification_application_simple.C

        # if [ $isotope = $bb ]; then
        #     sed -i -e 's@.*Long64.*@for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {@g' classification_application_simple.C
        # fi
        # if [ $isotope = $bb ]; then
        #     sed -i -e 's@.*Long64.*@for (Long64_t ievt=0; ievt<10000;ievt++) {@g' classification_application_simple.C
        # fi

    done
    sed -i -e 's@.*Long64.*@for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {@g' classification_application_simple.C

    root -l -b bdt_score_0nu_2nu.C
)

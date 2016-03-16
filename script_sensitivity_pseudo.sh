#!/bin/bash

# halflife_bdt_file=./halflives_bdt_pseudo.txt
# halflife_roi_file=./halflives_roi_pseudo.txt
halflife_file=./halflives_pseudo.txt
(
    count=0
    for pseudo in $(ls -rt ./pseudo_tree/merged/pseudo_1*)
    do
        echo $pseudo
        bdt_score_file="bdt_score_$count.root"
        sed -i -e 's@.*TString fname.*@TString fname = "'$pseudo'";@g' classification_application.C
        sed -i -e 's@.*bdt_score.*@TFile *target = new TFile ("'$bdt_score_file'","RECREATE");@g' classification_application.C

        root -l -q classification_application.C

        spectra_file="spectra_test_$count.root"
        # echo $spectra_file
        root -l -q 'spectra_pseudo.C("'$pseudo'","'$spectra_file'")'

        root -l -q 'sensitivity_pseudo.C("'$bdt_score_file'","'$spectra_file'","'$halflife_file'")'

        # echo $count
        rm $bdt_score_file
        rm $spectra_file
        let count++
    done
     # echo $count
)

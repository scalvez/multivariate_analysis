#!/bin/bash

(
    for isotope_str in 0nu 2nu
    do
        echo "Isotope : $isotope_str"
        isotope=$isotope_str
        bb='2nu'
        signal='0nu'
        n_step=5000
        n_0nu=20000
        n_2nu=100000


        if [ "$isotope" = "$signal" ]; then
            n_app=$(($n_0nu / $n_step))-1
            echo 'Computing '$n_0nu' events'
        fi
        if [ "$isotope" = "$bb" ]; then
            n_app=$(($n_2nu / $n_step))-1
            echo 'Computing '$n_2nu' events'
        fi
        sed -i -e 's@.*TString fname.*@TString fname = "./root_export_'$isotope_str'_25G.root";@g' classification_application_simple.C
        for (( i=0; i<=$n_app; i++ ))
        do
            # echo $i
            i_min=$(($i * $n_step))
            i_max=$((($i + 1) * $n_step))
            # echo $i_min
            # echo $i_max

            sed -i -e 's@.*bdt_score.*@TFile *target = new TFile ("bdt_score_'$isotope_str'_merge_'$i'.root","RECREATE");@g' classification_application_simple.C

            sed -i -e 's@.*Long64.*@for (Long64_t ievt='$i_min'; ievt<'$i_max';ievt++) {@g' classification_application_simple.C
        root -l -q classification_application_simple.C
        done


        hadd -f bdt_score_$isotope_str.root bdt_score_$isotope_str*_merge*.root
        rm bdt_score_*_merge_*.root
    done
    sed -i -e 's@.*Long64.*@for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {@g' classification_application_simple.C

    root -l -q bdt_score_0nu_2nu.C
)

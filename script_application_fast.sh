#!/bin/bash

(
    for isotope_str in 0nu 2nu tl208 bi214 radon
    do
        echo "Isotope : $isotope_str"
        isotope=$isotope_str
        bb='2nu'
        signal='0nu'
        tl='tl208'
        bi='bi214'
        rn='radon'
        n_step=5000
        n_0nu=100000
        n_2nu=1000000
        n_tl=120000
        n_bi=165000
        n_radon=20000

        sed -i -e 's@.*TString fname.*@TString fname = "./data_'$isotope_str'.root";@g' classification_application.C

        if [ "$isotope" = "$signal" ]; then
            n_app=$(($n_0nu / $n_step))-1
            echo 'Computing '$n_0nu' events'
        fi
        if [ "$isotope" = "$bb" ]; then
            n_app=$(($n_2nu / $n_step))-1
            echo 'Computing '$n_2nu' events'
        fi
        if [ "$isotope" = "$tl" ]; then
            n_app=$(($n_tl / $n_step))-1
            echo 'Computing '$n_tl' events'
        fi
        if [ "$isotope" = "$bi" ]; then
            n_app=$(($n_bi / $n_step))-1
            echo 'Computing '$n_bi' events'
        fi
        if [ "$isotope" = "$rn" ]; then
            n_app=$(($n_radon / $n_step))-1
            echo 'Computing '$n_radon' events'
        fi

        for (( i=0; i<=$n_app; i++ ))
        do
            # echo $i
            i_min=$(($i * $n_step))
            i_max=$((($i + 1) * $n_step))
            # echo $i_min
            # echo $i_max

            sed -i -e 's@.*bdt_score.*@TFile *target = new TFile ("bdt_score_'$isotope_str'_merge_'$i'.root","RECREATE");@g' classification_application.C

            sed -i -e 's@.*Long64.*@for (Long64_t ievt='$i_min'; ievt<'$i_max';ievt++) {@g' classification_application.C
        root -l -q classification_application.C
        done

        hadd -f bdt_score_$isotope_str.root bdt_score_$isotope_str*_merge*.root
        rm bdt_score_*_merge_*.root
    done
    sed -i -e 's@.*Long64.*@for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {@g' classification_application.C

    root -l -q bdt_score.C
)

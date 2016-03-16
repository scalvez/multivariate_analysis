#!/bin/bash

(
    shopt -s nullglob #empty arrays if no match
    two_nu_files=(2nu/sub_samples/*)
    for i in {0..99}
    do
        two_nu_1=${two_nu_files[$i]}
        two_nu_2=${two_nu_files[$i+100]}
        output="2nu/root_export_$i.root"
        command='~/.macros/merge_trees.C("snemodata","'$output'","'$two_nu_1'","'$two_nu_2'")'

        root -l -q $command
    done
)

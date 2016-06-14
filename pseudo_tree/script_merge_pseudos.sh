#!/bin/bash

(
    shopt -s nullglob #empty arrays if no match
    tl_files=(tl208/*)
    bi_files=(bi214/*)
    two_nu_files=(2nu/*)
    for i in {0..99}
    do
        tl=${tl_files[$i]}
        bi=${bi_files[$i]}
        two_nu=${two_nu_files[$i]}
        output="merged/pseudo_$i.root"
        command='~/.macros/merge_trees.C("snemodata","'$output'","'$two_nu'","'$bi'","'$tl'")'
        # command='~/.macros/merge_trees.C("snemodata","'$output'","'$two_nu'","'$bi'","'$tl'")'
        # command='~/.macros/merge_trees.C("snemodata","'$output'","'$bi'","'$tl'")'
        # echo $command
        root -l -q $command
    done
)

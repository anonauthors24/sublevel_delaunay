#! /bin/sh

BOTTLENECK_DIST="bottleneck_dist"
if [[ -f "$HOME/old_git/hera/bottleneck/bottleneck_dist" ]]; then
  BOTTLENECK_DIST="$HOME/old_git/hera/bottleneck/bottleneck_dist"
fi;

echo "*****Call flip_it"
./main -random ${1} ${2}
echo "*****Compute persistence of flip-complex"
./persistence_up_to_density output.scc ${3}
echo "*****Compute persistence of delaunay complex"
python3 ./alpha_from_gudhi.py last_instance.txt ${3}
echo "*****Bottleneck distance"
$BOTTLENECK_DIST  Gudhi_pers_pairs.txt Flipper_pers_pairs.txt

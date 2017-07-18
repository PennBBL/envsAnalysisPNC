#!/bin/bash

# directory of model
dir=$1
# comma separated threshold for overlay
thresh=$2
# list pattern for contrasts to run. a wild card goes after this name
contrasts=$3
image=cluster.sig.p05fdr.mgh # sig.mgh # s5000logp3.sig.masked.mgh # sig.mgh is the non cluster corrected map
threshname=$(echo $thresh | tr "," "-")

export SUBJECTS_DIR=/data/joy/BBL/studies/pnc/processedData/structural/freesurfer53
mkdir $dir/pngs 2>/dev/null

h=rh
himg=$(ls $SUBJECTS_DIR/fsaverage5/surf/${h}.inflated)
modelname=$(basename $dir)

for con in $(ls ${dir}/${contrasts}* -d); do
	imagename=$(echo $image | sed "s+.mgh++g")
	sig=$(ls $con/$image)
	conname=$( basename $con)
	# doing these both at once didn't work
	if [ "$h" == "lh" ]; then
		freeview -viewport 3d -zoom 1.8 -ss $dir/pngs/$modelname.lateral_${conname}_logp${threshname}_${imagename}.png -f $himg:overlay=$sig:overlay_threshold="$thresh":edgethickness=0
		freeview -viewport 3d -zoom 1.8 -cam Elevation 200 -ss $dir/pngs/${modelname}.medial_${conname}_logp${threshname}_${imagename}.png -f $himg:overlay=$sig:overlay_threshold="$thresh":edgethickness=0
	else
		freeview -viewport 3d -zoom 1.8 -cam Elevation -20 -ss $dir/pngs/$modelname.lateral_${conname}_logp${threshname}_${imagename}.png -f $himg:overlay=$sig:overlay_threshold="$thresh":edgethickness=0
		freeview -viewport 3d -zoom 1.8 -cam Elevation 180 -ss $dir/pngs/${modelname}.medial_${conname}_logp${threshname}_${imagename}.png -f $himg:overlay=$sig:overlay_threshold="$thresh":edgethickness=0
	fi
done

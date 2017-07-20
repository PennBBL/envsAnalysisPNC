#! /bin/bash

# Set parameters
meas=area
fwhm=20
template=fsaverage5

# Create model: for lm you have to write all terms of interaction, e.g. age + sex + age:sex
lm='"~ age + ageSqrd + sex + race + averageManualRating + pedu"' # used as argument for R script. since it's a string it must be in quotes (twice).

# Set paths
data=/data/joy/BBL/projects/envsMeduAnalysis/fs_analysis/freesurfer53/n1375_envs_demos.csv
scriptsdir=/data/joy/BBL/projects/envsMeduAnalysis/fs_analysis/freesurfer53
logdir=$scriptsdir/logs
rm -f $logdir/*
export SUBJECTS_DIR=/data/joy/BBL/studies/pnc/processedData/structural/freesurfer53

# Place where concatenated files will be stored
homedir=/data/joy/BBL/projects/envsMeduAnalysis/fs_analysis/freesurfer53
# Place where output is stored
outdir=$homedir

# Number and order of subjects
n=$(cat $data | wc -l)
n=$(echo "$n - 1" | bc)

name=$(echo $lm | sed "s/~//g" | sed "s/\"//g" | sed "s/ //g" | sed "s/+/_/g" | sed "s/:/BY/g" | sed "s/1/mean/g" )
echo $name
model=n${n}.$meas.$name.fwhm$fwhm.$template

# Where design files are stored
wd=$outdir/designs/$model

# Where model results are written
lhoutdir=$outdir/lh.$model
mkdir $lhoutdir 2>/dev/null

# Create design matrix and contrast matrices
Rvar=$(which R)
$Rvar --slave --file=$scriptsdir/design_matrix.R --args "$lm" "$data" "$wd"

# Create concatenated surface image for sample
## create subject call for mris_preproc command
subjs=""
for sub in $(cat $wd/subjlist.txt); do
	subjs=$(echo $subjs --s $sub)
done

## create image
if [ ! -e "$homedir/lh.$meas.fwhm$fwhm.$template.n$n.mgh" ]; then
	qsub -V -b y -q all.q -S /bin/bash -o $logdir -e $logdir mris_preproc $subjs --hemi lh --meas $meas --target $template --out $homedir/lh.$meas.fwhm$fwhm.$template.n$n.mgh
fi


# Store this file as a log
cp $scriptsdir/group_fs_analysis.sh $wd
cp $data $wd

cons=""
for i in $(ls $wd/contrast*.mat); do
	cons=$(echo "$cons --C $i ")
done

# Run the model
if [ -e "$homedir/lh.$meas.fwhm$fwhm.$template.n$n.mgh" ]; then
	mri_glmfit --glmdir $lhoutdir --y $homedir/lh.$meas.fwhm$fwhm.$template.n$n.mgh  --X $wd/X.mat $cons --surf fsaverage5 lh --fwhm $fwhm --save-yhat
	for img in $(ls $lhoutdir/contrast*/sig.mgh); do
		output=$(dirname $img)
		mri_surfcluster --in $img --subject fsaverage5 --fdr 0.05 --hemi lh --ocn $output/cluster.id.p05fdr.mgh --o $output/cluster.sig.p05fdr.mgh --sum $output/cluster.sum.p05fdr.txt
		mri_surfcluster --in $img --subject fsaverage5 --fdr 0.01 --hemi lh --ocn $output/cluster.id.p01fdr.mgh --o $output/cluster.sig.p01fdr.mgh --sum $output/cluster.sum.p01fdr.txt
		mris_convert -c $output/cluster.id.p05fdr.mgh $SUBJECTS_DIR/fsaverage5/surf/lh.sphere $output/cluster.id.p05fdr.asc
		mris_convert -c $output/cluster.id.p01fdr.mgh $SUBJECTS_DIR/fsaverage5/surf/lh.sphere $output/cluster.id.p01fdr.asc
	done
	mris_convert -c $lhoutdir/yhat.mgh $SUBJECTS_DIR/fsaverage5/surf/lh.sphere $lhoutdir/yhat.asc
	mris_convert -c $lhoutdir/rvar.mgh $SUBJECTS_DIR/fsaverage5/surf/lh.sphere $lhoutdir/rvar.asc
	cp $homedir/lh.$meas.fwhm$fwhm.$template.n$n.mgh  $lhoutdir
fi

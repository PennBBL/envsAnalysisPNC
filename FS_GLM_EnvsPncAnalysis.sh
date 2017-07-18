# LB 20170717
# Use this script to generate GLM for specific parameters

#! /bin/bash
meas=area
hemi=rh
template=fsaverage5
fwhm=20
analysis=rhGLMedu # other example: rhGLMses
lm='"~ age + ageSqrd + sex + race + averageManualRating + pedu"' # This is used as argument for R script. Must be in quotes twice because it is a string.
data=/data/joy/BBL/projects/envsMeduAnalysis/fs_analysis/rhGLMedu_ses/n1375_demos_20170717.csv
scriptsdir=/data/joy/BBL/projects/envsMeduAnalysis/fs_analysis/rhGLMedu_ses/scripts
logdir=$scriptsdir/logs

rm -f $logdir/*
export SUBJECTS_DIR=/data/joy/BBL/studies/pnc/processedData/structural/freesurfer53

# Place where concatenated files will be stored
homedir=/data/joy/BBL/projects/envsMeduAnalysis/fs_analysis/$analysis

# Place where output is stored
outdir=$homedir

# Number and order of subjects
n=$(cat $data | wc -l)
n=$(echo "$n - 1" | bc)

name=$(echo $lm | sed "s/~//g" | sed "s/\"//g" | sed "s/ //g" | sed "s/+/_/g" | sed "s/:/BY/g" | sed "s/1/mean/g" )
echo $name
model=n${n}.$meas.$name.$template

# Place where design files are stored
mkdir $outdir/designs/$model
wd=$outdir/designs/$model

# Place where model results are written
${hemi}outdir=$outdir/$hemi.$model
mkdir ${hemi}outdir 2>/dev/null

# Create a design matrix and contrast matrices
Rvar=$(which R)
$Rvar --slave --file=$scriptsdir/FS_GLM_designMatrix.R --args "$lm" "$data" "$wd"

# Create an image file to be analyzed
subjs=""
for sub in $(cat $wd/subjlist.txt); do
	subjs=$(echo $subjs --s $sub)
done

## Submit subject list to create the final surface file
if [ ! -e "$homedir/$hemi.$meas.fwhm$fwhm.$template.n$n.mgh" ]; then
	qsub -V -b y -q all.q -S /bin/bash -o $logdir -e $logdir mris_preproc $subjs --hemi $hemi --meas $meas --target $template --out $hemi.$meas.fwhm$fwhm.$template.n$n.mgh
fi

# Store this file as a log
cp $scriptsdir/FS_GLM_EnvsPncAnalysis.sh $wd
cp $data $wd

# Set the contrasts (created by design file)
cons=""
for i in $(ls $wd/contrast*.mat); do
	cons=$(echo "$cons --C $i ")
done

# Run the model
if [ -e "$homedir/$hemi.$meas.fwhm$fwhm.$template.n$n.mgh" ]; then
	mri_glmfit --glmdir $rhoutdir --y $homedir/$hemi.$meas.fwhm$fwhm.$template.n$n.mgh --X $wd/X.mat $cons --surf fsaverage5 rh --save-yhat

	# Cluster threshold at p<0.01 and p<0.05 500 mm^2
	for img in $(ls $rhoutdir/contrast*/sig.mgh); do
		output=$(dirname $img)
		mri_surfcluster --in $img --subject fsaverage5 --fdr 0.05 --hemi rh --ocn $output/cluster.id.p05fdr.mgh --o $output/cluster.sig.p05fdr.mgh --sum $output/cluster.sum.p05fdr.txt
		mri_surfcluster --in $img --subject fsaverage5 --fdr 0.01 --hemi rh --ocn $output/cluster.id.p01fdr.mgh --o $output/cluster.sig.p01fdr.mgh --sum $output/cluster.sum.p01fdr.txt
		mris_convert -c $output/cluster.id.p05fdr.mgh $SUBJECTS_DIR/fsaverage5/surf/rh.sphere $output/cluster.id.p05fdr.asc
		mris_convert -c $output/cluster.id.p01fdr.mgh $SUBJECTS_DIR/fsaverage5/surf/rh.sphere $output/cluster.id.p01fdr.asc
	done

	# Write cluster maps to csvs for plotting in R
	mris_convert -c $rhoutdir/yhat.mgh $SUBJECTS_DIR/fsaverage5/surf/rh.sphere $rhoutdir/yhat.asc
	mris_convert -c $rhoutdir/rvar.mgh $SUBJECTS_DIR/fsaverage5/surf/rh.sphere $rhoutdir/rvar.asc
	cp $homedir/$hemi.$meas.fwhm$fwhm.$template.n$n.mgh  $rhoutdir
fi

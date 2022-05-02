#!/bin/bash -e

# Processing script for neonatal dMRI data.
# 
# Usage:  ./preproc.sh sub-CC00069XX12 ses-26300
#
# 

# Global settings: force single-threaded processing
export MRTRIX_NTHREADS=2

# Paths
INPATHSUB=$INPATH/$1/$2/Dy-Di
OUTPATHSUB=$OUTPATH/$1/$2

# Find input data
INPUT=""
if [ -d $INPATHSUB ]; then
    INPUT=$(find $INPATHSUB -type f -name "*dhcp[5678]mbdti2sense_Re.nii")
fi
if [ -z "$INPUT" ]; then
    echo "Warning: input data not found for $1 $2 !"
    exit 1
fi
PHASE=${INPUT/Re.nii/RePh.nii}
CHI2E=${INPUT/Re.nii/Ch.nii}

# # Check input file dimensions; exit if not according to standard.
# if [ "$(mrinfo -size $INPUT -quiet)" != "100 100 64 300" ]; then
#     echo "Warning: incomplete dMRI data for $1 $2 ; skip"
#     exit 1
# fi

# Setup file directories and Lock for processing
if [ -e $OUTPATHSUB/lock ] || [ -e $OUTPATHSUB/premc.mif.md5 ]; then
    exit 1
else
    mkdir -p $OUTPATHSUB
    touch $OUTPATHSUB/lock
fi

echo "Run preprocessing in $1 $2"


# Log all commands and output
exec &>> $OUTPATHSUB/log.txt
set -x
echo "-- LAUNCHING PREPROCESSING --"
echo "host: $(hostname)"


# Check if dataset is complete
NVOL=$(mrinfo -size $INPUT -quiet | awk '{print $4}')
if [ $NVOL -eq 300 ]; then
    GRAD="dhcp300.txt"
    INDX="indx.txt"
else
    echo "Warning: incomplete dataset."
    GRAD="dhcp$NVOL.txt"
    head -n $NVOL dhcp300.txt > $GRAD
    INDX=$(mktemp -t indx-XXX.txt)
    cut -d' ' -f1-$NVOL indx.txt > $INDX 
fi


# Gibbs ringing removal
echo "---  GIBBS RINGING REMOVAL  ---"
mrdegibbs $INPUT - | mrconvert - - -stride -1,+2,+3,+4 | mrconvert - -grad $GRAD -import_pe_eddy acqp.txt $INDX $OUTPATHSUB/premc.mif


# Calculate voxel weights
echo "---  CALCULATING VOXEL WEIGHTS  ---"
# Detrend Log-Chi2 error; run MRF-GMM segentation; then smooth output
LOGRES=$(mrcalc $CHI2E -log - | mrfilter - smooth -fwhm 15,15,0 - | mrcalc $CHI2E -log - -sub - | mrconvert -coord 3 0:$((NVOL-1)) - -)
TMPVOXW=$(mktemp -t voxweights-XXXXX.mif)
./voxorgmm.py $LOGRES $TMPVOXW
mrfilter $TMPVOXW smooth -fwhm 1.5,1.5,0 $OUTPATHSUB/voxelweights.mif.gz
rm $TMPVOXW $LOGRES


echo "--- ESTIMATING THE FIELD MAP ---"
# Select cleanest b=0 volumes
PHASEMIF=$(mrconvert $PHASE -grad $GRAD -import_pe_eddy acqp.txt $INDX -)
TMPMASK=$(dwiextract -bzero $OUTPATHSUB/premc.mif - | mrmath -axis 3 - mean - | mrthreshold - -)
VOLS=$(./selectb0volumes.py $PHASEMIF $TMPMASK)
rm $PHASEMIF $TMPMASK

# Extract selected images for Topup
refdir=$(mktemp -dt topupref-XXXXX)
mrconvert $OUTPATHSUB/premc.mif -coord 3 $VOLS - | mrconvert - $refdir/ref.nii -export_pe_table $refdir/b0acqp.txt -stride -1,+2,+3,+4 

# Run FSL Topup
time fsl5.0-topup --imain=$refdir/ref.nii --datain=$refdir/b0acqp.txt --config=b02b0 --out=$refdir/topup_out --fout=$refdir/fieldmap --iout=$refdir/b0ref

# save field map and selected volumes
mrconvert $refdir/fieldmap.nii.gz -set_property fieldidx $VOLS $OUTPATHSUB/fieldmap.mif.gz


# Brain extraction
echo "---  BRAIN MASKING  ---"
mrmath -axis 3 $refdir/b0ref.nii.gz mean $refdir/b0.nii
fsl5.0-bet2 $refdir/b0.nii $refdir/b0m -m
mrconvert $refdir/b0m_mask.nii.gz $OUTPATHSUB/reconmask.mif.gz

# clean up temporary files
rm -r $refdir


# Checksum and compression
echo "---  COMPRESSION  ---"
md5sum $OUTPATHSUB/premc.mif > $OUTPATHSUB/premc.mif.md5
gzip $OUTPATHSUB/premc.mif || :

# Done
echo "--- PREPROCESSING DONE  ---"
rm $OUTPATHSUB/lock


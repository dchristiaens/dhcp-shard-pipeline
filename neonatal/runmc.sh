#!/bin/bash -e

# Run motion correction for neonatal dMRI data.
# 
# Usage:  ./runmc.sh sub-CC00069XX12 ses-26300
#
# 

# Paths
INPATHSUB=$OUTPATH/$1/$2
OUTPATHSUB=$OUTPATH/$1/$2

# Check input
if [ ! -e $INPATHSUB/premc.mif.md5 ] || [ ! -e $INPATHSUB/fieldmap.mif.gz ]; then
    echo "ERROR: input files not found."
    exit 1
fi

# Lock for processing
if [ -e $OUTPATHSUB/lock ] || [ -e $OUTPATHSUB/postmc-mssh.mif.md5 ]; then
    exit 1
else
    mkdir -p $OUTPATHSUB
    touch $OUTPATHSUB/lock
fi

echo "Motion correction in $1 $2"


# Log all commands and output
exec &>> $OUTPATHSUB/log.txt
set -x
echo "-- LAUNCHING MOTION CORRECTION --"
echo "host: $(hostname)"


echo "---  MOTION CORRECTION  ---"

# Crop fixed inputs if incomplete data
NVOL=$(mrinfo -size $INPATHSUB/premc.mif.gz -quiet | awk '{print $4}')
if [ $NVOL -eq 300 ]; then
    # Select prior weights of correct size
    PRW="priorweights.txt"
    # Rescale for spin history
    tmpfile=$(mrcalc $INPATHSUB/premc.mif.gz scaling-recon07.mif -multiply -)
else
    # Select prior weights of correct size
    PRW=$(mktemp -t priorweights-XXX.txt)
    cut -d' ' -f1-$NVOL priorweights.txt > $PRW
    # Rescale for spin history
    tmpfile=$(mrconvert -coord 3 0:$((NVOL-1)) scaling-recon07.mif - | mrcalc $INPATHSUB/premc.mif.gz - -multiply -)
fi

# Retrieve field map alignment index
FIDX=$(mrinfo $INPATHSUB/fieldmap.mif.gz -property fieldidx | cut -d ',' -f1)

# SHARD reconstruction
time dwimotioncorrect $tmpfile -mask $INPATHSUB/reconmask.mif.gz -lmax 0,4,6,8 $OUTPATHSUB/postmc-mssh.mif -setup neonatal.config -sspfile ssp.txt -fieldmap $INPATHSUB/fieldmap.mif.gz -fieldidx $FIDX -rlmax 4,2,0 -mb 4 -sorder 3,2 -priorweights $PRW -voxelweights $INPATHSUB/voxelweights.mif.gz -export_motion $OUTPATHSUB/motion.txt -export_weights $OUTPATHSUB/sliceweights.txt -scratch $TMPDIR

# Evaluate DWIs 
mrinfo $INPATHSUB/premc.mif.gz -export_grad_mrtrix grad.txt -force
mssh2amp $OUTPATHSUB/postmc-mssh.mif grad.txt $OUTPATHSUB/postmc-dwi.mif -nonnegative


echo "---  COMPRESSION  ---"
md5sum $OUTPATHSUB/postmc-mssh.mif > $OUTPATHSUB/postmc-mssh.mif.md5
md5sum $OUTPATHSUB/postmc-dwi.mif > $OUTPATHSUB/postmc-dwi.mif.md5
gzip $OUTPATHSUB/postmc-mssh.mif || :
gzip $OUTPATHSUB/postmc-dwi.mif || :


echo "---  STAGE 2 DONE  ---"
rm $OUTPATHSUB/lock



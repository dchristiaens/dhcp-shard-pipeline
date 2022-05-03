#!/bin/bash -e

# Run slice intensity modulation correction (dStripe) for neonatal dMRI data.
# 
# Requires dStripe Docker image, see https://github.com/maxpietsch/dStripe
# 
# Usage:  ./rundstripe.sh sub-CC00069XX12 ses-26300
# 
# Note: script follows the pipeline run for release 3 but is untested
#


# set up environment
source local-config.sh

INPATHSUB=$OUTPATH/$1/$2
OUTPATHSUB=$DSOUTPATH/$1/$2

# skip if already done
[ -e "${OUTPATHSUB}"/postmc_dstriped-dwi.nii.gz ] && exit 0

# Check input
if [ ! -e $INPATHSUB/postmc-mssh.mif.gz ] || [ ! -e $INPATHSUB/reconmask.mif.gz ]; then
    echo "ERROR: input files not found."
    exit 1
fi

# dStripe Docker image setup
image=maxpietsch/dstripe:1.0
if [[ ! "$(docker images -q "${image}" 2> /dev/null)" ]]; then
  echo dStripe Docker image is missing
  exit 1
fi

device=0 # comma separated for multi-GPU or "cpu"
batch_size=30 # 30 is fine for use on Tesla V100 (32GB), reduce if using the CPU or if GPU is out of memory
grad=$(dirname "$0")/dhcp300.txt

# Log all commands and output
exec &>> $OUTPATHSUB/log.txt
set -x
echo "-- LAUNCHING ${image} for $1 $2 --"
echo "host: $(hostname)"

# project to full 300 volumes
mssh2amp "$OUTPATHSUB"/postmc-mssh.mif.gz "${grad}" "${OUTPATHSUB}"/amp300.mif

# estimate destripe field
docker run --rm --volume "${OUTPATHSUB}/":/data "${image}" \
  dwidestripe /data/amp300.mif /data/reconmask.mif.gz /data/dstripe_field.mif.gz \
  -device "${device}" -batch_size "${batch_size}" -scratch "${TMPDIR}"

# apply field
mrcalc "${OUTPATHSUB}"/amp300.mif "${OUTPATHSUB}"/dstripe_field.mif.gz -mult - | mrconvert - -stride +1,+2,+3,+4 "${OUTPATHSUB}"/postmc_dstriped-dwi300.mif

# truncate dwi to original length
vols=( $(mrinfo $OUTPATHSUB/postmc-dwi.mif.gz -size) )
vols=${vols[3]}
if [[ $vols == 300 ]]; then
	mrconvert "${OUTPATHSUB}"/postmc_dstriped-dwi300.mif "${OUTPATHSUB}"/postmc_dstriped-dwi.nii.gz -export_grad_fsl "${OUTPATHSUB}"/postmc_dstriped-dwi.bvec "${OUTPATHSUB}"/postmc_dstriped-dwi.bval -export_grad_mrtrix "${OUTPATHSUB}"/postmc_dstriped-dwi.grad_mrtrix
else
	mrconvert "${OUTPATHSUB}"/postmc_dstriped-dwi300.mif -coord 3 0:$(( vols - 1 )) "${OUTPATHSUB}"/postmc_dstriped-dwi.nii.gz -export_grad_fsl "${OUTPATHSUB}"/postmc_dstriped-dwi.bvec "${OUTPATHSUB}"/postmc_dstriped-dwi.bval -export_grad_mrtrix "${OUTPATHSUB}"/postmc_dstriped-dwi.grad_mrtrix

# cleanup
rm "${OUTPATHSUB}"/amp300.mif "${OUTPATHSUB}"/dstripe_field.mif.gz "${OUTPATHSUB}"/postmc_dstriped-dwi300.mif

#!/usr/bin/env bash
#. /appcg/data2/SPECULOOSPipeline/.bashrc

set -o nounset
set -o errexit
set -o pipefail

abspath() {
    python -c "import os; print(os.path.realpath('${1}'))"
}

# Which tasks to run. Set to "1" if the task should be run otherwise "0".
# Initialize all tasks to 1 by default
T1="1" # create input lists, default: 1
T2="1" # create masterbias, default: 1
T3="1" # create masterdark, default: 1
#T4="1" # copy temporary shutter map, default: 1
T5="1" # create masterflat, default: 1
T6="1" # reduce science images, default: 1
T7="1" # check if each fits image has astrom in it's header, default:1
T8="1" # create stack image, default: 1
T9="1" # perform photometry, default: 1
T10="1" # find SPXXXXXX ID of target star and plot diff LC of target
#readonly T11="0" # water vapour correction
T12="1" # generate PDF report for the night

# Parse and remove task control arguments, leaving positional args intact
filtered_args=()
only_task=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --no_T1) T1="0"; shift ;;
        --no_T2) T2="0"; shift ;;
        --no_T3) T3="0"; shift ;;
        --no_T5) T5="0"; shift ;;
        --no_T6) T6="0"; shift ;;
        --no_T7) T7="0"; shift ;;
        --no_T8) T8="0"; shift ;;
        --no_T9) T9="0"; shift ;;
        --no_T10) T10="0"; shift ;;
        --no_T12) T12="0"; shift ;;
        --only_T1) only_task="1"; shift ;;
        --only_T2) only_task="2"; shift ;;
        --only_T3) only_task="3"; shift ;;
        --only_T5) only_task="5"; shift ;;
        --only_T6) only_task="6"; shift ;;
        --only_T7) only_task="7"; shift ;;
        --only_T8) only_task="8"; shift ;;
        --only_T9) only_task="9"; shift ;;
        --only_T10) only_task="10"; shift ;;
        --only_T12) only_task="12"; shift ;;
        *) filtered_args+=("$1"); shift ;; # Keep non-task arguments
    esac
done

# Handle --only_TX logic
if [[ -n "$only_task" ]]; then
    T1="0"; T2="0"; T3="0"; T5="0"; T6="0"; T7="0"; T8="0"; T9="0"; T10="0"; T12="0"
    case $only_task in
        1) T1="1" ;;
        2) T2="1" ;;
        3) T3="1" ;;
        5) T5="1" ;;
        6) T6="1" ;;
        7) T7="1" ;;
        8) T8="1" ;;
        9) T9="1" ;;
        10) T10="1" ;;
        12) T12="1" ;;
    esac
fi

# Reassign positional parameters with filtered arguments
set -- "${filtered_args[@]}"

# Make variables readonly after parsing
readonly T1 T2 T3 T5 T6 T7 T8 T9 T10 T12

#if [[ $# -ne 8 ]] && [[ $# -ne 8 ]]; then
#    cat >&2 <<-EOF
#  Usage: $0 <runname> <root-directory> <input-catalogue> <initial-wcs-solution> <confidence-map> <shuttermap> <wcsfit-reference-frame> [master-flat]
#  Argument descriptions:
#  * runname
#  The name of the run to use. This is so multiple runs can be performed on different
#  data, and the outputs can be unique.
#  * root-directory
#  This is the location of the input files. The directory structure must be as follows:
#  root-directory
#      Telescope
#          images
#              <one directory per date>
#                    IMAGE*.fits.bz2.fts
#  * input-catalogue
#  The list of coordinates to place apertures at
#  * initial_wcs_solution
#  The initial wcs solution computed by Toms MCMC code to compute distortion
#  parameters
#  * confidence-map
#  * shuttermap
#  * date in the format YYYYMMDD
#  * wcsfit-reference-frame
#  * master-flat
#  Custom master flat to use, overriding the flat computed from the supplied data
#EOF
#    exit 1
#fi

# command line arguments
readonly SCRIPTDIR=${ORCHARD_PATH:-/Users/matthewhooton/PycharmProjects/orchard/src}
readonly RUNNAME=${1}
#readonly GIVEN_INPUTCATALOGUE=$3
readonly WCSSOLUTION=''
readonly DATES=$3
readonly CTHRESH=$4
readonly STHRESH=$5
readonly TEL=$6
readonly INPUT_TARGETS=${7:-}

#readonly WCSFIT_REFERENCE_FRAME=$5
readonly WCSFIT_REFERENCE_FRAME=''
#readonly TARGET=''

#MASTER_FLAT=${9:-}
#if [[ ! -z "${MASTER_FLAT}" ]]; then
#    MASTER_FLAT="$(abspath "${MASTER_FLAT}")"
#fi

#GIVEN_INPUTCATALOGUE=${9:-}
#if [[ ! -z "${GIVEN_INPUTCATALOGUE}" ]]; then
#    GIVEN_INPUTCATALOGUE="$(abspath "${GIVEN_INPUTCATALOGUE}")"
#fi


# there are ** cores on the server
echo "$(python -c "import multiprocessing; print(multiprocessing.cpu_count())")"
#readonly CORES=$(python -c "import multiprocessing; print multiprocessing.cpu_count()"-1)
readonly CORES=1
#readonly CORES=$(($(python -c "import multiprocessing; print multiprocessing.cpu_count()") / 2))
#if
#readonly CORES=32
#$(python -c "import multiprocessing; print multiprocessing.cpu_count()")
readonly APSIZE=4
readonly NUMSTACK=50
readonly IPIX=6
readonly XMATCH='vizgaia3'
readonly OUTEXT="fits"
readonly VERSION='v2'
#readonly TLIST="/appct/data/SPECULOOSPipeline/tests/target_list_ids_201905.txt"
readonly TLIST=$(abspath ${2})/ml_40pc.txt
#readonly EXT='fts'

if [ "$TEL" = "ARTEMIS" ]; then
  readonly GAIN=1.1
elif [ "$TEL" = "SAINT-EX" ]; then
  readonly GAIN=3.48
else
  readonly GAIN=1.0029
fi

echo "Using ${CORES} cores"
#echo "$(ldd /appct/data/SPECULOOSPipeline/miniconda2/bin/solve-field)"

#readonly T8="0" # run image subtraction, default: 0
#readonly T9="0" # detrend, default: 1
#readonly T10="0" # detrend with lightcurves, default: 1c
#readonly T11="0" # Make qa plots, default: 1
#readonly T13="0" # create lightcurve plotting data for each image, default: 1
#readonly T15="0" # correct the shape of the master flat, default: 0
#readonly T16="0" # QC on bias and dark images, default: 0
#readonly T18="0" #  plot differential LC of target


# Zero Level Pipeline
# Here all the commands are listed.
# this script can be run  from command line. to do the whole pipeline.


# for every image on a night classify as science, bias, dark or flat.
create_input_lists() {
    echo "START T1"
    echo "Create lists with Images"
    ensure_directory "${OUTPUTDIR}/${DATE}"
    CMD="python ${SCRIPTDIR}/main/createlists.py \"$IMGDIRS\" ${DATDIR} \"${OUTPUTDIR}\" \"${REPORTDIR}\" IMAGE ${EXT} $RUNNAME $DATE $TEL"
    echo $CMD
    ${CMD}
    ret=$?
    echo $ret
#        if [ $ret -ne 0 ]; then
#            echo "reduction failed"
#            continue
#             #Handle failure
#             #exit if required
#        fi
    echo "END T1"
}

# get the names of all targets observed this night
find_target_names(){
    IMAGELISTS=${OUTPUTDIR}/${DATE}/reduction/${RUNNAME}_image_*.list
    TARGET=()  # Initialize as an empty array

    for IMAGELIST in ${IMAGELISTS}
    do
        IMAGELIST=${IMAGELIST#${OUTPUTDIR}/${DATE}/reduction}
        SUBLIST=${IMAGELIST#/${RUNNAME}_image_}
        TARGET_NAME=${SUBLIST%.list}
        TARGET+=("${TARGET_NAME}")  # Add to array with proper quoting
    done

    echo "${TARGET[@]}"  # Output array elements with proper quoting
}

# Replace the space_to_underscore function with a space_to_double_hyphen function
space_to_double_hyphen() {
    echo "$1" | sed 's/ /--/g'
}

#correct_target_names(){
#
#    # save the targets in an array
#    CORRECT_TARGET=""
#
#    for j in ${TARGET[*]}
#    do
#        IMAGELISTS=${OUTPUTDIR}/${DATE}/${j}/${RUNNAME}/*_processed.dat
#        if $(any_filelists ${IMAGELISTS}); then
#            for IMAGELIST in ${IMAGELISTS}
#            do
#                CORRECT_TARGET+=("$(python ${SCRIPTDIR}/correct_target_names.py ${IMAGELIST} ${j})")
#            done
#        fi
#    done
#
#    echo ${CORRECT_TARGET[*]}
#}


# find the extension of images in directory - to ensure we have fits/fts files
find_ext(){

    IMAGELISTS=(${IMGDIRS}/*.*)

    for IMAGELIST in ${IMAGELISTS[@]}
    do
        filename=$(basename -- "$IMAGELIST")
        if [ "${filename##*.}" != "log" ] && [ "${filename##*.}" != "Z" ] && [ "${filename##*.}" != "txt" ] && [ "${filename##*.}" != "gz" ]
        then
            extension="${filename##*.}"
        fi
        filename="${filename%.*}"
    done
    echo ${extension}
}

# find the extension of images in directory - to ensure we have fits/fts files
find_ext2(){

    IMAGELISTS=(${IMGDIRS}/*/*.*)

    for IMAGELIST in ${IMAGELISTS}
    do
          filename=$(basename -- "$IMAGELIST")
          extension="${filename##*.}"
          filename="${filename%.*}"
    done

    echo ${extension}
}

create_master_bias() {
    echo "START T2"
    # Create MasterBias
    echo "Create MasterBias"
    CMD="python ${SCRIPTDIR}/calibration/pipebias.py $BIASLIST ${RUNNAME}_MasterBias.fits ${OUTPUTDIR}/${DATE}/reduction ${REPORTDIR} ${GAIN} ${RUNNAME} ${i}"
    echo ${CMD}
    ${CMD}

    CMD="python ${SCRIPTDIR}/reporting/report.py ${report} ${DATE} ${i} ${TEL} ${RUNNAME} 2"
    ${CMD}
    echo "END T2"
}

create_master_dark() {
    echo "START T3"
    #Create MasterDark
    echo "Create MasterDark"
    CMD="python ${SCRIPTDIR}/calibration/pipedark.py $DARKLIST ${RUNNAME}_MasterBias.fits ${RUNNAME}_MasterDark.fits ${OUTPUTDIR}/${DATE}/reduction ${REPORTDIR} ${GAIN} ${RUNNAME} ${i}" #/Reduction/output/${RUNNAME}
    echo ${CMD}
    ${CMD}

    CMD="python ${SCRIPTDIR}/reporting/report.py ${report} ${DATE} ${i} ${TEL} ${RUNNAME} 3"
    ${CMD}
    echo "END T3"
}
#
#qc_bias_dark(){
#
#    echo "Perform Quality Checks on Bias and Dark Images"
#    python ${SCRIPTDIR}/reduction/QCBiasDark.py \
#            $OUTPUTDIR \
#            -s 20180221 \
#            -e ${DATE}
#}


create_master_flat() {
    echo "START T5"
    #Create MasterFlat
    echo "Create MasterFlat"
    CMD="python ${SCRIPTDIR}/calibration/pipeflat.py $FLATLIST ${T2} ${T3} ${RUNNAME}_MasterBias.fits ${RUNNAME}_MasterDark.fits ${RUNNAME}_MasterFlat.fits ${OUTPUTDIR}/${DATE}/reduction ${REPORTDIR}"
    echo ${CMD}
    ${CMD}

    CMD="python ${SCRIPTDIR}/reporting/report.py ${report} ${DATE} ${i} ${TEL} ${RUNNAME} 5"
    ${CMD}
    echo "END T5"
}

#copy_custom_master_flat() {
#    local MFLAT_DEST=${WORKINGDIR}/Reduction/output/${RUNNAME}/${RUNNAME}_MasterFlat.fits
#    echo "Copying custom master flat ${MASTER_FLAT} => ${MFLAT_DEST}"
#    cp "${MASTER_FLAT}" "${MFLAT_DEST}"
#}

reduce_images() {
    # Helper function to reduce a list of lists of images
    # Function submits jobs asynchronously and returns the list of job names
    #   used to run the analysis so the wait step can halt other processing.
    counter="0"

    MASTERFLATS=${OUTPUTDIR}/${DATE}/reduction/${RUNNAME}_MasterFlat_*.fits
    MASTERFLAT=""
    for M in ${MASTERFLATS}
    do
        MFLAT=${M#${OUTPUTDIR}/${DATE}/reduction/}
        echo ${MFLAT}
        if [ -z "$MASTERFLAT" ]; then
            MASTERFLAT="$MFLAT"
        else
            MASTERFLAT="$MASTERFLAT $MFLAT"
        fi
    done

    # Use target name directly since files are created with underscores
    IMAGELIST=${OUTPUTDIR}/${DATE}/reduction/${RUNNAME}_image_${i}.list

    printf "For target = ${i}\n"
    ensure_directory "${OUTPUTDIR}/${DATE}/${i}/${RUNNAME}" #${IMAGELIST%.*}
    CMD="python ${SCRIPTDIR}/calibration/pipered.py ${IMAGELIST} --biasname ${RUNNAME}_MasterBias.fits --darkname ${RUNNAME}_MasterDark.fits --flatname $MASTERFLAT --caldir ${OUTPUTDIR}/${DATE}/reduction --outdir ${OUTPUTDIR}/${DATE}/${i}/${RUNNAME} --gain ${GAIN} --version ${VERSION} --usebias 1 --usedark 1"
    echo ${CMD}
    ${CMD}
}


any_filelists() {
    local readonly IMAGELISTS=$1
    ls ${IMAGELISTS} 2>/dev/null >/dev/null
}

# calibrate the science images using master bias/dark/flat images
reduce_science_images() {
    echo "START T6"
    echo "Reduce Science Images"
    IMAGELISTS=${OUTPUTDIR}/${DATE}/reduction/${RUNNAME}_image_*.list
    if $(any_filelists ${IMAGELISTS}); then
        reduce_images "${IMAGELISTS}" #"${1}"
        CMD="python ${SCRIPTDIR}/reporting/report.py ${report} ${DATE} ${i} ${TEL} ${RUNNAME} 6"
        ${CMD}
        echo "END T6"
    else
        echo "Reduction failed for ${i}"
        continue
    fi
#    if ~ls ${OUTPUTDIR}/${DATE}/${i}/${RUNNAME}/*_processed.dat &>/dev/null; then
#        echo "Reduction failed for ${i}"
#        continue
#    else
#        CMD="python ${SCRIPTDIR}/reporting/report.py ${report} ${DATE} ${i} ${TEL} ${RUNNAME} 6"
#        ${CMD}
#        echo "END T6"
#    fi
}

# perform astrometric fit of each image using ASTROMETRY.NET code
check_astrometry(){
    echo "START T7"
    printf "\n**Check if astrometry in header of images**\n"
#    for i in ${TARGET[*]}
#    do
    printf "For Target = ${i}\n"
    IMAGELISTS=${OUTPUTDIR}/${DATE}/${i}/${RUNNAME}/*_processed.dat

    if $(any_filelists ${IMAGELISTS}); then
        for IMAGELIST in ${IMAGELISTS}
        do
            if [ -e "$IMAGELIST" ]; then
                echo ${IMAGELIST}
                CMD="python ${SCRIPTDIR}/astrom/astrometry.py
                    $IMAGELIST
                    --nproc=${CORES}\
                    --ext=${OUTEXT}"
                 echo ${CMD}
                ${CMD}
            else
                echo "${IMAGELIST} does not exist for ${i}"
            fi
        done

        CMD="python ${SCRIPTDIR}/reporting/report.py ${report} ${DATE} ${i} ${TEL} ${RUNNAME} 7"
        ${CMD}
        echo "END T7"

    else
        echo "${IMAGELISTS} does not exist for ${i}"
        continue
    fi
}

# create stack image to generate catalogue of field stars
create_stack_image() {
    echo "START T8"
    printf "\n**Create Stack Image**\n"
    ensure_directory "${OUTPUTDIR}/StackImages"

    if [ -d "${DATDIR}/catcache" ]; then
        rm -rf ${DATDIR}/catcache
    fi
    IMAGELISTS=${OUTPUTDIR}/${DATE}/${i}/${RUNNAME}/*_processed.dat


    if $(any_filelists ${IMAGELISTS}); then

        for IMAGELIST in ${IMAGELISTS}
        do
            echo ${IMAGELIST}
            SUBLIST=${IMAGELIST#${OUTPUTDIR}/${DATE}/${i}/${RUNNAME}/}
            FILTER=${SUBLIST%_processed.dat}

#            correct_i="$(python ${SCRIPTDIR}/correct_target_names.py ${IMAGELIST} ${i})"
#            echo "${correct_i}"

            #check if stack catalogue already exists
            if [ -e ${OUTPUTDIR}/StackImages/${correct_i}_stack_catalogue_${FILTER}.${OUTEXT} ]; then
                echo "The stack catalogue for target " ${correct_i} " with filter " ${FILTER} " already exists"
            else
#                python ${SCRIPTDIR}/reporting/email_eon.py \
#                ${DATE} \
#                -r ${RUNNAME} \
#                -rep ${report} \
#                -t ${i} \
#                -tel ${TEL} \
#                -e 0

                CMD="python ${SCRIPTDIR}/photometry/ZLP_create_cat.py
                    $IMAGELIST
                    ${correct_i}_outstack_${FILTER}.${OUTEXT}
                    ${correct_i}_stack_catalogue_${FILTER}.${OUTEXT}
                    ${OUTPUTDIR}/StackImages
                    ${REPORTDIR}
                    $FILTER
                    $DATE
                    False
                    $CTHRESH
                    $STHRESH
                    False
                    $NUMSTACK
                    $IPIX
                    $XMATCH
                    $APSIZE
                    $CORES
                    $OUTEXT"
                echo ${CMD}
                ${CMD}

            fi
        done

    else
        continue
    fi
#    done
#    fi
    CMD="python ${SCRIPTDIR}/reporting/report.py ${report} ${DATE} ${i} ${TEL} ${RUNNAME} 8"
    ${CMD}
    echo "END T8"
}

wait_for_jobs() {
    if [ -z "$1" ]
    then
        echo "Error in invokation; wait_for_jobs <jobids>" >&2
        exit 1
    fi

    JOBIDS="${1}"

    echo "Wait until jobs '${JOBIDS}' finish"
    qsub -hold_jid "${JOBIDS}" -N WAIT  -sync y -cwd ${WORKINGDIR}/wait.sh
}

iterate_and_act_on_lists() {
    local readonly lists="$1"
    local readonly action="$2"
    #local readonly target="$3"

    if $(any_filelists ${lists}); then
        for fname in ${lists}; do
            eval "${action} ${fname}"
        done
    fi
}

# perform aperture photometry on a single target
single_perform_aperture_photometry() {
#    local filelist=$1
    local readonly filelists=${OUTPUTDIR}/${DATE}/${1}/${RUNNAME}/*_processed.dat
    local readonly output_directory=${OUTPUTDIR}/${DATE}/${1}

    if $(any_filelists ${filelists}); then
        for filelist in ${filelists}
        do
            if [ -d "${DATDIR}/catcache" ]; then
                rm -rf ${DATDIR}/catcache
            fi
            SUBLIST=${filelist#${OUTPUTDIR}/${DATE}/${1}/${RUNNAME}/}
            FILTER=${SUBLIST%_processed.dat}
            echo "For filter = ${FILTER}"
            echo "Using Stack Catalogue ${OUTPUTDIR}/StackImages/${2}_stack_catalogue_${FILTER}.${OUTEXT}"

            CMD="python ${SCRIPTDIR}/photometry/ZLP_app_photom.py \
                --catfile ${OUTPUTDIR}/StackImages/${2}_stack_catalogue_${FILTER}.${OUTEXT}\
                --filter ${FILTER} \
                --nproc ${CORES} \
                --filelist ${filelist} \
                --outdir ${output_directory} \
                --date ${DATE} \
                --apsize ${APSIZE} \
                --ext ${OUTEXT}\
                --ipix ${IPIX} \
                --s_thresh ${STHRESH} \
                --catsrc ${XMATCH}"
            echo ${CMD}
            ${CMD}
                #--wcsref ${WCSFIT_REFERENCE_FRAME}

            #PIPELINESHA=$(extract_pipeline_sha $(dirname $0))
            condense_photometry ${filelist}_phot ${2} ${FILTER}
        done
    fi
}

condense_photometry() {
    #Condense the photometry
   # for every target create an output.fts file for that specific night in the output/date/target folder
   # also check if the overall output.fts file exists for that target and if it does then append this night's date to it
   # if not create a new output.fts file for the target

    local filelist=$1
    local outname="output.${OUTEXT}"
    local filter=$3

    printf "\n**Condensing images into one ${outname}**\n"

    CMD="python ${SCRIPTDIR}/condense/zlp_condense.py \
        --outputdir ${OUTPUTDIR}
        --obsdir ${OBSDIR}
        --target ${2}
        --oldtarget ${i}
        --outname ${outname}
        $( cat ${filelist} | sed 's/.${OUTEXT}/.${OUTEXT}.phot/' )
        --date ${DATE}
        --filter ${FILTER}
        --ext ${OUTEXT}
        --version ${VERSION}
        --tlist ${TLIST}"
    echo "${CMD}"

    # add to global output fits file:
    if [ -e ${OUTPUTDIR}/${outname} ]; then
        echo "UPDATE"
        CMD="${CMD} --update"
    fi

    ${CMD}
}

perform_aperture_photometry() {
    echo "START T9"
    printf "\n**Running aperture photometry**\n"

    printf "For Target = ${i}\n"
    single_perform_aperture_photometry ${i} ${correct_i}

    CMD="python ${SCRIPTDIR}/reporting/report.py ${report} ${DATE} ${i} ${TEL} ${RUNNAME} 9"
    ${CMD}
    echo "END T9"
}

#run_detrending() {
#    SYSREM=sysrem
#    if hash ${SYSREM} 2>/dev/null; then
#        echo "Detrending with SYSREM"
#
#        local readonly photomfile=$(find ${OUTPUTDIR}/output -name '${targ}_output.fts')
#        if [ ! -z "${photomfile}" ]; then
#            local readonly output_directory=$(dirname $photomfile)
#            local readonly outfile=${output_directory}/tamout.fits
#            echo "Running sysrem to create ${outfile}"
#
#            if [ -z ${NOSYSREM:-} ]; then
#                ${SYSREM} ${photomfile} ${outfile}
#                python ${BASEDIR}/scripts/combine_with_sysrem.py -v -p ${photomfile} -t ${outfile}
#            else
#                echo '*** sysrem has been disabled with envar: NOSYSREM. unset to run sysrem'
#            fi
#        else
#            echo "Cannot find photometry output files" >&2
#        fi
#    else
#        echo "Cannot find sysrem binary ${SYSREM}" >&2
#    fi
#}

perform_differential_photometry(){
    echo "Plot the differential LC of the target star"
    local readonly output=${OUTPUTDIR}/${DATE}/${i}/*_${DATE}_output.${OUTEXT}
    ensure_directory "${OUTPUTDIR}/lightcurves"

    for out in ${output}
    do
      SUBLIST=${out#${OUTPUTDIR}/${DATE}/${i}/}
      THING1=${SUBLIST%_${DATE}_output.${OUTEXT}}
      filt=${THING1#*_}
      gaia_id=${THING1%_${filt}}
      echo ${filt}
      echo ${gaia_id}

      if [[ "${gaia_id}" == "*" ]]; then
        echo "No output fits files found, skip this target!"

      else
#        j=5
#        CMD="python ${SCRIPTDIR}/diff_photometry/SPlightcurve.py \
#          ${out} \
#          -o ${OUTPUTDIR} \
#          -l ${OUTPUTDIR}/${DATE}/${i}/lightcurves \
#          -f ${filt} \
#          -b 5 \
#          -d ${DATE} \
#          -a ${j} \
#          -v ${VERSION} \
#          -targ ${TLIST} \
#          --basedir ${BASEDIR}"  # Add this line
#        echo ${CMD}
#        ${CMD}

        for ((j=3;j<=8;j++));
        do
            CMD="python ${SCRIPTDIR}/diff_photometry/SPlightcurve.py \
                ${out} \
                -o ${OUTPUTDIR} \
                -l ${OUTPUTDIR}/${DATE}/${i}/lightcurves \
                -f ${filt} \
                -b 5 \
                -d ${DATE} \
                -a ${j} \
                -v ${VERSION} \
                -targ ${TLIST}"
            echo ${CMD}
            ${CMD}
        done
      fi
    done
#    done

}

#water_vapour_correction(){
#
#    echo "Performing water vapour correction"
#
#    for ((j=3;j<=8;j++));
#    do
#        python ${SCRIPTDIR}/water_vapour/water_vapour_nightly.py \
#            -t ${TEL} \
#            -d ${DATE} \
#            -a ${j} \
#            /appct/data/SPECULOOSPipeline/
#
#    done
#}

#PDF reports
pdf_report(){
    echo "Create PDF report for the night"

    # Convert TARGET array items from double-hyphens to spaces for display
    DISPLAY_TARGETS=""
    for t in ${TARGET[*]}; do
        if [ -z "$DISPLAY_TARGETS" ]; then
            # Convert back from double-hyphens to spaces
            DISPLAY_TARGETS="$(echo "$t" | sed 's/--/ /g')"
        else
            DISPLAY_TARGETS="$DISPLAY_TARGETS $(echo "$t" | sed 's/--/ /g')"
        fi
    done

    CMD="python ${SCRIPTDIR}/reporting/pdf_report_catriona.py \
        --datdir \"${DATDIR}\" \
        --obsdir \"${OBSDIR}\" \
        --date \"${DATE}\" \
        --target \"${DISPLAY_TARGETS}\" \
        --ap \"5\" \
        --telescope \"${TEL}\" \
        --version \"${VERSION}\""
    echo "${CMD}"
    eval "${CMD}"
}
# Some helper functions
ensure_directory() {
    DIR=${1}
    test -d ${DIR} || mkdir -p ${DIR}
}

setup_environment() {
#    if [ -z ${DISABLE_ANACONDA:-} ]; then
#        # Allow the user to override the anaconda path variable
#        if [ -z ${ANACONDA_PATH:-} ]; then
#        # If anaconda is available, use it
#        case `hostname` in
#            ngtshead*)
#            ANACONDA_PATH=/home/sw/anaconda
#            ;;
#            *)
#            ANACONDA_PATH=${HOME}/anaconda
#            ;;
#        esac
#        fi
#
#        PARANAL_ANACONDA_PATH=/usr/local/anaconda
#
#        if [[ -d ${ANACONDA_PATH} ]]; then
#            export PATH=${PARANAL_ANACONDA_PATH}/bin:${ANACONDA_PATH}/bin:${PATH}
#        fi
#    fi
#    source ~/.bashrc

    set +o nounset

    echo "Using python: $(which python)"
    echo "Using WCSFIT: $(which wcsfit)"
    echo "Using IMSTACK: $(which imstack)"
    echo "Using IMCORE: $(which imcore)"
#    echo "Setting # of user processes:"
#    ulimit -u 2061801

#    case "$(hostname -s)" in
#        ngts*)
#            export IERS_DATA=/usr/local/pipeline/data
#            export JPLEPH_DATA=${IERS_DATA}/linux_p1550p2650.430t
#            ;;
#        mbp*)
#            export IERS_DATA=${HOME}/.local/data
#            export JPLEPH_DATA=${IERS_DATA}/linux_p1550p2650.430t
#            ;;
#    esac

#    if [ ! -z ${IERS_DATA} ]; then
#        echo "IERS data path: ${IERS_DATA}"
#    fi
#
#    if [ ! -z ${JPLEPH_DATA} ]; then
#        echo "JPLEPH data path: ${JPLEPH_DATA}"
#    fi

    # LD_LIBRARY_PATH for sysrem
#    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/appcg/data2/SPECULOOSPipeline/python27/lib"
#    LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}:/opt/intel/composer_xe_2013_sp1.0.080/compiler/lib/intel64

    echo "Environment set up"
    set -o nounset
}

setup_directory_structure() {
    for subdir in ${OUTPUTDIR} ${OUTPUTDIR}/${DATE} ${OUTPUTDIR}/${DATE}/reduction; do #AperturePhot Reduction Reduction/output; do
        local dirpath=${subdir}
        ensure_directory ${subdir}
    done
}

#run_lightcurves_detrending() {
#    # Aperture index for casu lightcurves:
#    #  0 = auto (bad)
#    #  1 = 1 * rcore
#    #  2 = sqrt(2) * rcore
#    #  3 = 2 * rcore
#    local ap_index=1
#
#    if hash lightcurves-casu 2>/dev/null; then
#        local readonly ref=$GIVEN_INPUTCATALOGUE
#        local readonly photomfile=$(find ${WORKINGDIR}/AperturePhot/output -name 'output.fits' -print)
#        if [ ! -z "${photomfile}" ]; then
#            local readonly output_directory=$(dirname $photomfile)
#            local readonly outfile=${output_directory}/casu-lightcurves-out.fits
#            local readonly number_of_coefficients=2
#            local readonly source_files_dir=${WORKINGDIR}/Reduction/output/${RUNNAME}
#            echo "Running casu lightcurves file to create ${outfile}"
#
#            lightcurves-casu -f ${number_of_coefficients} -a ${ap_index} -o ${outfile} -p ${ref} $(find ${source_files_dir} -name 'proc*.phot')
#            python ${BASEDIR}/scripts/combine_with_casu_detrended.py -v -p ${photomfile} -d ${outfile}
#        fi
#    else
#        echo "Cannot find CASU lightcurves binary" >&2
#    fi
#}
#
#generate_qa_plots() {
#    bash ${BASEDIR}/scripts/zlp-qa/run.sh \
#        ${WORKINGDIR} \
#        ${WORKINGDIR}/QualityAssessment
#}

extract_pipeline_sha() {
    local readonly dirname="$1"
    (cd $dirname && git rev-parse HEAD)
}

#function contains() {
#    local n=$#
#    local value=${!n}
#    for ((i=1;i < $#;i++)) {
#        if [ "${!i}" == "${value}" ]; then
#            echo "n"
##            return 0
#        fi
#    }
#    echo "y"
##    return 1
#}

containsElement () {
  local e match="$1"
  shift
  for e; do [[ "$e" == "$match" ]] && return 0; done
  return 1
}

# Do photometry on subtracted Images

main() {
    setup_environment
    setup_directory_structure

    cd ${DATDIR}
    date
    echo "Using ${CORES} cores"

    EXT="$(find_ext)"

    echo "Extension: ${EXT}"

    if [ "${EXT}" != "fts" ] && [ "${EXT}" != "fits" ] && [ "${EXT}" != "fz" ]
    then
      EXT="$(find_ext2)"
    fi

    echo "${EXT}"
    QC="${REPORTDIR}/QC"
    report="${REPORTDIR}/Pipeline_Report"
    count=0

    if [ "${EXT}" == "fts" ] || [ "${EXT}" == "fits" ] || [ "${EXT}" == "fz" ]
    then

        [ "$T1" = "1" ] && create_input_lists

        # Use input targets if provided, otherwise fall back to find_target_names
        if [ -n "$INPUT_TARGETS" ]; then
            TARGET=($INPUT_TARGETS)
            echo "Using specified targets: ${TARGET[@]}"
        else
            TARGET="$(find_target_names)"
            echo "Using discovered targets: ${TARGET[@]}"
        fi

        # if there are no science images then attempt to create master calibration images and produce report
        if [[ " ${TARGET[@]} " =~ "logs" ]]; then
            echo "There are no science images for ${DATE}!"
            TARGET="N/A"
            echo $TARGET
            [ "$T2" = "1" ] && create_master_bias
            [ "$T3" = "1" ] && create_master_dark
            [ "$T5" = "1" ] && create_master_flat
            [ "$T12" = "1" ] && pdf_report
        else
            echo $TARGET
            for i in ${TARGET[*]}
            do
                printf "\nTARGET: ${i}\n"
                CMD="python ${SCRIPTDIR}/reporting/report.py ${report} ${DATE} ${i} ${TEL} ${RUNNAME} 1"
                ${CMD}

                count=$((count+1))

                #cd ${WORKINGDIR}/Reduction
                [ "$T2" = "1" ] && [ $count = 1 ] && create_master_bias
                [ "$T3" = "1" ] && [ $count = 1 ] && create_master_dark
                [ "$T5" = "1" ] && [ $count = 1 ] && create_master_flat
#                if [[ ! -z "${MASTER_FLAT}" ]]; then
#                    copy_custom_master_flat
#                fi
                [ "$T6" = "1" ] && reduce_science_images

                IMAGELISTS=${OUTPUTDIR}/${DATE}/${i}/${RUNNAME}/*_processed.dat
                if $(any_filelists ${IMAGELISTS}); then
                    for IMAGELIST in ${IMAGELISTS}
                    do
                        correct_i="$(python ${SCRIPTDIR}/utils/correct_target_names.py ${IMAGELIST} ${i} ${OUTPUTDIR})"
                        correct_i=$(echo "${correct_i}" | tr ' ' '--')
                        gaia_id="$(python ${SCRIPTDIR}/utils/gaia_id_from_schedule.py ${OBSDIR} ${DATE} ${i} --silent)"
                        correct_i="$(python ${SCRIPTDIR}/utils/already_run_targs.py ${TLIST} ${OUTPUTDIR} ${i} ${gaia_id})"
                    done
                fi

                [ "$T7" = "1" ] && check_astrometry
                [ "$T8" = "1" ] && create_stack_image
                [ "$T9" = "1" ] && perform_aperture_photometry
                [ "$T10" = "1" ] && perform_differential_photometry

                CMD="python ${SCRIPTDIR}/reporting/email_eon.py \
                ${DATE} \
                -r ${RUNNAME} \
                -rep ${report} \
                -t ${i} \
                -tel ${TEL} \
                -e 1"
                echo ${CMD}
                ${CMD}
            done

            [ "$T12" = "1" ] && pdf_report
        fi
    else
        if [[ "${EXT}" == *"output"* ]]; then
          echo "There are no images for ${DATE}!"
          TARGET="N/A"
          [ "$T12" = "1" ] && pdf_report
        else
            echo "Extension was found to be ${EXT}. Please remove the corrupted file from ${IMGDIRS}"
        fi
    fi

#    [ "$T18" = "1" ] && plot_diff_lc

#    [ "$T9" = "1" ] && run_detrending
#
#    [ "$T10" = "1" ] && run_lightcurves_detrending
#
#    [ "$T11" = "1" ] && generate_qa_plots

    date
    echo "PIPELINE COMPLETE"
}

echo "PIPELINE VERSION ${VERSION}"
echo "For Telescope: ${TEL}"
ALL_DATES=(${DATES})
echo "All Dates: ${DATES[*]}"
for DATE in ${DATES[*]}
do
  echo "D: ${DATE}"
done

for DATE in ${DATES[*]}
do
    echo "For Date: ${DATE}"
    BASEDIR=$(abspath ${2})
    OBSDIR=${BASEDIR}/Observations/${TEL}
    DATDIR=${BASEDIR}/PipelineOutput/${VERSION}/${TEL}
    IMGDIRS=${OBSDIR}/images/${DATE} #** #*/*
    OUTPUTDIR=${DATDIR}/output
    REPORTDIR=${DATDIR}/reports
    ensure_directory "${REPORTDIR}"
    ensure_directory "${OUTPUTDIR}"

    BIASLIST=${OUTPUTDIR}/${DATE}/reduction/${RUNNAME}_bias.list
    DARKLIST=${OUTPUTDIR}/${DATE}/reduction/${RUNNAME}_dark.list
    FLATLIST=${OUTPUTDIR}/${DATE}/reduction/${RUNNAME}_flat_*.list

    main 2>&1 | tee ${DATDIR}/logs/${DATE}_${RUNNAME}_${VERSION}${INPUT_TARGETS:+_$(echo $INPUT_TARGETS | tr ' ' '_')}.log
done
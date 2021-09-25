#!/bin/bash

GNUPLOTBIN=${GNUPLOT-gnuplot}
HAVE_GS=0
HAVE_IM=0
HAVE_JPEG_TERM=0
sopt=0
NAME=""

preflight() {

	# gnuplot is mandatory
	test -x "$(command -v $GNUPLOTBIN)"
	if [ 0 -ne $? ] ; then
		echo "ERROR: gnuplot is required"
		exit 1
	fi

	echo "set term" | $GNUPLOTBIN 2>&1 | grep -wq jpeg
	if [ 0 -eq $? ] ; then
		HAVE_JPEG_TERM=1
	fi

	test -x "$(command -v gs)"
	if [ 0 -eq $? ] ; then
		HAVE_GS=1
	fi

	test -x "$(command -v convert)"
	if [ 0 -eq $? ] ; then
		HAVE_IM=1
	fi

	if [ $HAVE_JPEG_TERM -eq 0 -a $HAVE_GS -eq 0 -a $HAVE_IM -eq 0 ] ; then
		echo "ERROR: Ghostscript or ImageMagik is required"
		exit 1
	fi

}

usage() {
	local me=$(basename $0)
	echo "Usage: $me [OPTIONS] FILE"
	echo "Options:"
	echo "    -h          Show help"
	echo "    -N          Specify output file base name"
	echo "    -s          Sequential processing"
	echo "    -t 'TITLE'  Space-separated variable names"
}

preflight

while getopts "c:hN:st:T:" o; do
    case "${o}" in
        c)
			declare -a cols=( ${OPTARG} )
            ;;
        h)
            usage
			exit 0
            ;;
		N)
			NAME="${OPTARG}"
			;;
        s)
            sopt=1
            ;;
        t)
			declare -a title=( ${OPTARG} )
            ;;
		T)
			TITLE="${OPTARG}"
			;;
        *)
            usage
			exit 1
            ;;
    esac
done

shift $((OPTIND-1))

FNAME="$1"

test -f "$FNAME" || {
	usage
	exit 1
}

NCOLS=$(egrep -v '^#' $FNAME | head -1 | awk '{ print NF }')

ID=$(echo $(basename $FNAME) | sed 's/\..*$//')_$(date +%Y%m%d%H%M%S)_$(printf "%05d" $RANDOM)

[ -z "$NAME" ] || ID=$NAME

xplot() {
	local fname=$1
	local col=$2
	local title=$3
	local j=$(( $col - 1 ))
	local f=${ID}_x$(printf "%02d" $j).jpg
	local ptitle=$(echo $ID | sed 's:_: :g')

	[ -z "$TITLE" ] || ptitle="$TITLE"

	if [ $HAVE_JPEG_TERM -ne 0 ] ; then
		echo "set term jpeg size 1280,960 ; set output '$f' ; set title '$ptitle' ; plot '$fname' u 1:$col w l ls 1 title 'x(${j}) $title'" | $GNUPLOTBIN > /dev/null 2>&1 || exit 1
	else
		local tmp=$(mktemp --tmpdir=/tmp tmp.XXXXX.ps)
		echo "set term postscript ; set output '$tmp' ; set title '$ptitle' ; plot '$fname' u 1:$col w l title 'x(${j}) $title'" | $GNUPLOTBIN > /dev/null 2>&1 || exit 1
		if [ $HAVE_GS -ne 0 ] ; then
			gs -q -dSAFER -dBATCH -dNOPAUSE -sDEVICE=jpeg -r300 \
			   -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -dMaxStripSize=8192 \
			   -sOutputFile=$f -dAutoRotatePages=/None \
			   -c '<</Orientation 3>> setpagedevice' -f $tmp || exit 1
		else
			convert -density 300 $tmp -flatten -background white -resize 1024 -rotate 90 $f || exit 1
		fi
		rm -f $tmp
	fi
}

if [ 0 -eq $sopt ] ; then
	maxjobs=$(nproc)
	(( maxjobs-- ))
else
	maxjobs=1
fi

for i in $(seq 2 $NCOLS) ; do
	j=$(( $i - 2))
	c=$(( $i - 1))

	if [ 0 -eq ${#cols[@]} ] ; then
		plot=1
	else
		plot=0
		for k in ${cols[@]} ; do
			if [ $c -eq $k ] ; then
				plot=1
				break
			fi
		done
	fi

	if [ 1 -eq $plot ] ; then
	    xplot $FNAME $i ${title[$j]} &
	    while [ $(jobs -r | wc -l) -gt $maxjobs ] ; do
		sleep 10
	    done
	fi

done

wait

exit 0



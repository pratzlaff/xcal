basedir=/data/legs/rpete/flight/xcal
datadir=/data/legs/rpete/flight/xcal/data

sources() {
    echo mkn421 pks2155
}

obsid2source() {
    local obsid="$1"
    for src in `sources`
    do
	source2obsids $src | grep -q $obsid && { echo $src; return; }
    done
    echo "unrecognized obsid=$obsid" 1>&2
    return 1
}

band_range() {
    local band="$1"
    case "$band" in
	0) echo 0.5 8.0 ; return ;;
	1) echo 0.33 0.54 ; return ;;
	2) echo 0.54 0.85 ; return ;;
	3) echo 0.85 1.5 ; return ;;
	4) echo 1.5 4.0 ; return ;;
	5) echo 4.0 10.0 ; return ;;
	6) echo 0.26 0.514 ; return ;;
	*) echo "band '$band' is unrecognized" 1>&2 ; return 1 ;;
    esac
}

source2obsids()
{
    local src="$1"
    grep '^[0-9]' $basedir/obsids_${src} | cut -f 1 | grep -v 16474
}

instruments()
{
    local obsid="$1"
    local src=`obsid2source $obsid` ||  { echo "unrecognized obsid=$obsid" 1>&2; return 1; }
    local f=$(ls $datadir/"$src"/"$obsid"/tg_reprocess/*_evt2.fits 2>/dev/null)

    if [ -z "$f" ]
    then
	grep "$obsid" $basedir/obsids_${src} | perl -anle 'print $F[-1]'
	return
    fi

    echo $(detnam "$f")/$(grating "$f")
}
    
asol_stack()
{
    local obsid="$1"
    local source=$(obsid2source $obsid)

    ls $datadir/$source/$obsid/primary/pcadf*asol1.fits* | perl -le 'chomp(@f=<>); @t=map { /pcadf(\d{9})/ } @f; @t{@t}=(); @o=(); print join(",", map { $t=$_; (grep/\Q$t\E/,@f)[-1]} sort keys %t)'
}

pbk_select()
{
    local dir="$1"
    ls $dir/*pbk0.fits* | tail -1
}

msk_select()
{
    local dir="$1"
    ls $dir/acisf*msk1.fits* | tail -1
}

detnam()
{
    local evt2="$1"
    punlearn dmkeypar
    dmkeypar "$evt2" detnam echo+
}

grating()
{
    local evt2="$1"
    punlearn dmkeypar
    dmkeypar "$evt2" grating echo+
}

evt2_bin()
{
    local evt2="$1"
    local bin="$2"

    local dir=`dirname "$evt2"`
    local olist=$(tg_order_list "$(detnam "$evt2")" "$(grating "$evt2")")
    local i=0

    punlearn dmcopy

    : ${PERL:='perl'}
    $PERL $basedir/bintimes.pl "$evt2" "$bin" |
    while read start stop elapsed
    do
	local outevt2="$dir/evt2_"$(printf "%02d" $i)'.fits'
	local outpha2="$dir/pha2_"$(printf "%02d" $i)'.fits'
	dmcopy "$evt2""[time=$start:$stop]" "$outevt2" clobber=yes opt=all
	tgextract "$outevt2" "$outpha2" \
	    outfile_type=pha_typeII \
	    tg_srcid_list=all \
	    tg_part_list=header_value \
	    inregion_file=CALDB \
	    tg_order_list="$olist" \
	    ancrfile=none \
	    respfile=none \
	    clobber=yes
	((i++))
	[[ $i > 99 ]] && { echo "stopping with i=$i" 1>&2; return 1; }
    done
}

tg_order_list()
{
    local det="$1"
    local grating="$2"

    case $det in
	ACIS*) echo -1,-2,-3,1,2,3 ; return ;;
	HRC*) echo -1,1,-2,2,-3,3,-4,4,-5,5,-6,6,-7,7,-8,8,-9,9,-10,10 ; return ;;
	*) echo "tg_order_list() - unknown detector '$det'" 1>&2 ; return 1 ;;
    esac
}

xspec_fit()
{
    local newgain=''
    [ $NEWGAIN -eq 1 ] && newgain=_newgain

    local src="$1"
    local obsid="$2"
    local elow=0.5
    local ehigh=8.0

    local tgdir=$datadir/$src/$obsid/tg_reprocess$newgain
    local outdir=$datadir/$src/$obsid/fit$newgain
    mkdir -p $outdir

    local pha2=`ls "$tgdir"/hrcf*_pha2.fits`
    local base=`echo "$pha2" | sed -e 's/_pha2.fits//'`
    local base=`basename "$base"`
    local postfix=_LEG

    # ugly
    export HEADAS=/usr/local/heasoft-6.13/x86_64-unknown-linux-gnu-libc2.17-0
    . $HEADAS/headas-init.sh
    rm -f $tgdir/${base}_{npha,bkg}2.fits
    /data/legs/rpete/flight/acis_letg_cedge/grating2xspec \
                "$pha2" both

    for ostr in neg pos
    do
	case "$ostr" in
	    neg) title="$obsid: TG\\\\_M = -1" ; row=1 ; x='-' ;;
	    pos) title="$obsid: TG\\\\_M = +1" ; row=2 ; x='' ;;
	    *) echo "unrecognized ostr = '$ostr'" 1>&2 ; return 1 ;;
	esac
    
	local outxcm=$outdir/xsfit${ostr}${postfix}.xcm
	local outps=$outdir/xsfit${ostr}${postfix}.ps/cps
	local outlog=$outdir/xsfit${ostr}${postfix}.log

	cat - <<EOP > $outxcm
set dir $tgdir
set base $base
set maxorder 10 
set elow $elow
set ehigh $ehigh

data \${dir}/\${base}_npha2.fits{${row}}
back \${dir}/\${base}_bkg2.fits{${row}}

for { set i 1 } { \$i <= \$maxorder } { incr i } {
  response \$i: \${dir}/LEG_${x}\$i.rmf
  arf \$i: \${dir}/\${base}LEG_${x}\${i}_garf.fits
}

ignore 1:0.-\${elow},\${ehigh}-**
statistic lstat
#abund wilm
model 1:m1 tbabs*logpar
.0161    -0.01
2.0       0.01
0.5       0.01
1.0      -0.01
0.3       0.01

for { set i 2 } { \$i <= \$maxorder } { incr i } {
  model  \$i:m\$i tbabs*logpar
  = m1:1
  = m1:2
  = m1:3
  = m1:4
  = m1:5
}
fit
log $outlog
error 1.0 m1:2,3,5
flux \${elow},\${ehigh}
show all
log none
cpd ${outps}
for { set i 1 } { \$i <= \$maxorder } { incr i } {
  set j [expr \$i+2]
  setplot command color [expr \$i+1] on \$j
}
setplot energy
setplot rebin 20 100
plot data ratio

exit
EOP

	xspec - $outxcm
    
    done

    # ugly!!!
    unset LD_LIBRARY_PATH
}

sherpa_fit()
{
    local newgain=''
    [ $NEWGAIN -eq 1 ] && newgain=_newgain

    local obsid="$1"
    local src=`obsid2source "$obsid"` || { echo "unrecognized obsid=$obsid" 1>&2; return 1; }

    local tgdir=$datadir/$src/$obsid/tg_reprocess$newgain

    local band
#    for band in 0 1 2 3 4 5
    for band in 6
    do
	local fitdir=$datadir/$src/$obsid/fit$newgain/B"$band"

	rm -rf "$fitdir"
	mkdir -p "$fitdir"

	local inst=`instruments "$obsid"`
	case "$inst" in
	    *HRC*)

		sherpa_fit_generic \
		    $obsid \
		    $datadir/$src/$obsid \
		    LEG \
		    1 \
		    neg \
		    $band

		[ $band -ne 6 ] && \
		sherpa_fit_generic \
		    $obsid \
		    $datadir/$src/$obsid \
		    LEG \
		    2 \
		    pos \
		    $band
		;;

	    *LETG*)

		sherpa_fit_generic \
		    $obsid \
		    $datadir/$src/$obsid \
		    LEG \
		    1 \
		    neg \
		    $band \
		    $tgdir/acisf*LEG_-1_garf.fits \
		    $tgdir/LEG_-1.rmf

		[ $band -ne 6 ] && \
		sherpa_fit_generic \
		    $obsid \
		    $datadir/$src/$obsid \
		    LEG \
		    4 \
		    pos \
		    $band \
		    $tgdir/acisf*LEG_1_garf.fits \
		    $tgdir/LEG_1.rmf
		;;
	    # row component order
	    #   1   HEG      -1
	    #   4   HEG      +1
	    #   7   MEG      -1
	    #  10   MEG      +1

	    *HETG*)
		sherpa_fit_generic \
		    $obsid \
		    $datadir/$src/$obsid \
		    HEG \
		    1 \
		    neg \
		    $band \
		    $tgdir/acisf*HEG_-1_garf.fits \
		    $tgdir/HEG_-1.rmf
		sherpa_fit_generic \
		    $obsid \
		    $datadir/$src/$obsid \
		    HEG \
		    4 \
		    pos \
		    $band \
		    $tgdir/acisf*HEG_1_garf.fits \
		    $tgdir/HEG_1.rmf
		sherpa_fit_generic \
		    $obsid \
		    $datadir/$src/$obsid \
		    MEG \
		    7 \
		    neg \
		    $band \
		    $tgdir/acisf*MEG_-1_garf.fits \
		    $tgdir/MEG_-1.rmf
		sherpa_fit_generic \
		    $obsid \
		    $datadir/$src/$obsid \
		    MEG \
		    10 \
		    pos \
		    $band \
		    $tgdir/acisf*MEG_1_garf.fits \
		    $tgdir/MEG_1.rmf
		;;
	    *) echo "unrecognized instruments='$inst', returning" 1>&2 ; return 1 ;;
	esac

    done # for band
}

sherpa_fit_generic() {

    local newgain=''
    [ $NEWGAIN -eq 1 ] && newgain=_newgain

    local obsid="$1"
    local dir="$2"
    local arm="$3"
    local row="$4"
    local ostr="$5"
    local band="$6"
    local garf="$7"
    local rmf="$8"

    read elow ehigh <<<$( band_range "$band" )

    case "$ostr" in
	neg) local title="$obsid: TG\\_M = -1"; local x='-x';;
	pos) local title="$obsid: TG\\_M = +1"; local x='x';;
	*) echo "unrecognized ostr = '$ostr'" 1>&2 ; return 1 ;;
    esac
    
    local inst=`instruments "$obsid"` || return 1
    local resptext
    case "$inst" in
	*ACIS*)
	    resptext=$(cat <<EOF
    load_arf("$garf")
    load_arf("$garf", bkg_id=1)
    load_arf("$garf", bkg_id=2)

    load_rmf("$rmf")
EOF
		)
	    ;;
	*HRC*)
	    local pha2=`ls "$dir"/tg_reprocess"$newgain"/hrcf*_pha2.fits`
	    local base=`echo "$pha2" | sed -e 's/_pha2.fits//'`
	    resptext=$(cat <<EOF
    orders=range(1, 11)
    garfs=[ "${base}LEG_%d_garf.fits" % ($x) for x in orders ]
    rmfs=[ "${dir}/tg_reprocess${newgain}/LEG_%d.rmf" % ($x) for x in orders ]

    load_multi_arfs(garfs, orders)
    load_multi_rmfs(rmfs, orders)
EOF
		)
	    ;;
	*) echo "shouldn't be here" 1>&2; return 1 ;;
    esac

    local fitdir=$dir/fit$newgain/B$band
    local outpy=$fitdir/${arm}_${ostr}.py

    cat - <<EOP > $outpy
import myfit

import glob
import os
import re
import sherpa.utils.err
from sherpa.astro.ui import *
from pychips import *

base = '$dir'

for d in ['$fitdir', '$fitdir/binned']:
    try: os.mkdir(d)
    except: pass

pha2 = glob.glob(base+'/tg_reprocess${newgain}/*_pha2.fits')
pha2.extend( glob.glob(base+'/tg_reprocess${newgain}/pha2_[0-9][0-9].fits'))

for pha2 in pha2:

    out='$fitdir'
    match = re.match( r'.*/pha2_(\\d+).fits', pha2 )
    if match:
        out = out + '/binned/${arm}_${ostr}_' + match.groups()[0]
    else:
        out = out + '/${arm}_${ostr}'

    set_default_id($row)
    load_pha(pha2, use_errors=True)

$resptext

    snr=5
    group_snr(snr, bkg_id=1)
    group_snr(snr, bkg_id=2)
    group_snr(snr)

    subtract()

    elow=$elow
    ehigh=$ehigh
    ignore(0., 10000.)
    notice(elow, ehigh)

    myfit.source_def("$src", $band)

    fit()
    show_fit( outfile=out+'.params', clobber=True )

    flux = calc_energy_flux(elow, ehigh)
    f=open( out+'.flux', 'w' )
    f.write("energy flux = %s\n" % flux)
    f.close()

    try:
	conf()
	show_conf( outfile=out+'.conf', clobber=True )
    except sherpa.utils.err.EstErr:
	pass

    show_all( outfile=out+'.show_all', clobber=True )

    myfit.plot_range(elow, ehigh, "$title")

    print_window( out+'.pdf',  [ "clobber", True, 'fittopage', True, 'orientation', 'landscape'] )
#    print_window( out+'.ps',  [ "clobber", True, 'fittopage', True, 'orientation', 'landscape'] )
#    print_window( out+'.png', [ "clobber", True ] )
EOP

    sherpa -b $outpy

    # rotate PDF files
    for f in $fitdir/${arm}_${ostr}.pdf $fitdir/binned/${arm}_${ostr}_??.pdf 
    do
	pdf270 $f -o $f
    done
}

hrcs_qefile() {
    local pha2="$1"
    local qe="$2"

    local qedir=$CALDB/data/chandra/hrc/qe
    [ $qe -gt 14 ] && qedir=/data/legs/rpete/flight/hrcs_qe/post_N0014

    if [ $qe -lt 12 ]
    then
	printf "$qedir/hrcsD1999-07-22qeN%04d.fits" $qe
	return
    fi

    local obsid=$(dmkeypar $pha2 obs_id echo+)
    local dateobs=$(dmkeypar $pha2 date-obs echo+)
    local qefile

    case $obsid in
	14238)
	    qefile=$(printf "$qedir/hrcsD2012-03-29qeN%04d.fits" $qe)
	    ;;
	14234|14396|14397)
	    qefile=$(printf "$qedir/hrcsD1999-07-22qeN%04d.fits" $qe)
	    ;;
	*)
	    qefile=$(printf "$qedir/hrcsD1999-07-22qeN%04d.fits" $qe)
	    [ $dateobs \> 2012-03-29 ] && qefile=$(printf "$qedir/hrcsD2012-03-29qeN%04d.fits" $qe)
	    ;;
    esac

    echo $qefile
}

hrcs_garfs() {
    local obsid="$1"
    local qe="$2"

    local source=$(obsid2source $obsid)

    local qeufile=

    case $obsid in
	14238)
	    qeufile=$CALDB/data/chandra/hrc/qeu/hrcsD2012-03-29qeuN0007.fits ;;
	14234|14396|14397)
	    qeufile=$CALDB/data/chandra/hrc/qeu/hrcsD2012-01-01qeuN0007.fits ;;
	*) ;;
    esac

    local o=$obsid
    local s=$source

    local outdir=$datadir/$s/$o/tg_reprocess/hrcs_garfs_qe$qe
    mkdir -p $outdir

    local pha2=$(ls $datadir/$s/$o/tg_reprocess/hrcf*_pha2.fits)
    local evt2=$(ls $datadir/$s/$o/tg_reprocess/hrcf*_evt2.fits)
    local dtf1=$(ls $datadir/$s/$o/primary/hrcf*_dtf1.fits* | tail -1)
    local bpix1=$(ls $datadir/$s/$o/secondary/hrcf*_bpix1.fits* | tail -1)
    local asol=$(asol_stack $obsid)

    local qefile=$(hrcs_qefile $pha2 $qe)

    local order=1
    while [ $order -le 10 ]
    do
	for junk in neg pos
	do
	    if [ $junk = neg ]
	    then
		orderstr=-
		rmf=$outdir/LEG_-$order.rmf
		row=$(( $order*2-1 ))
	    else
		orderstr=
		rmf=$outdir/LEG_$order.rmf
		row=$(( $order*2 ))
	    fi

	    punlearn ardlib

	    pset ardlib AXAF_HRC-S1_QE_FILE=$qefile'[AXAF_QE1]'
	    pset ardlib AXAF_HRC-S2_QE_FILE=$qefile'[AXAF_QE2]'
	    pset ardlib AXAF_HRC-S3_QE_FILE=$qefile'[AXAF_QE3]'

	    if [ -n "$qeufile" ]
	    then
		pset ardlib AXAF_HRC-S1_QEU_FILE=$qeufile'[AXAF_QEU1]'
		pset ardlib AXAF_HRC-S2_QEU_FILE=$qeufile'[AXAF_QEU2]'
		pset ardlib AXAF_HRC-S3_QEU_FILE=$qeufile'[AXAF_QEU3]'
	    fi

	    if [ $qe -lt 14 ]
	    then
		pset ardlib AXAF_LETG_1111_LSF_FILE=$CALDB/data/chandra/hrc/lsfparm/hrcsleg${orderstr}1D1999-07-22lsfparmN0003.fits
	    fi

	    punlearn mkgrmf

	    mkgrmf \
		grating_arm=LEG \
		order=$orderstr$order \
		outfile=$rmf \
		srcid=1 \
		detsubsys=HRC-S2 \
		threshold=1e-06 \
		obsfile=$pha2'[SPECTRUM]' \
		regionfile=$pha2 \
		wvgrid_arf=compute \
		wvgrid_chan=compute \
		verbose=0 \
		clobber=yes

	    punlearn mkgarf
	    punlearn fullgarf
	    fullgarf \
		$pha2 \
		$row \
		$evt2 \
		$asol \
		"grid($rmf[cols ENERG_LO,ENERG_HI])" \
		$dtf1 \
		$bpix1 \
		$outdir/ \
		maskfile=NONE \
		clobber=no

	done # for pos and neg

	(( order++ ))
    done # for each order

    punlearn ardlib

    rm $outdir/_*

}

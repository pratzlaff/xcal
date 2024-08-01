#! /usr/bin/perl

use warnings;
use strict;
use Astro::FITS::CFITSIO;
use PDL;

@ARGV==2 or die;
my ($f, $bin) = @ARGV;

my $h = Astro::FITS::CFITSIO::fits_read_header( $f. '[1]' );

my $gti =
  $h->{DETNAM} =~ /ACIS-\d*7\d*/ ? 'gti7' :
  $h->{DETNAM} =~ /HRC/ ? 'gti' :
  die $h->{DETNAM};

my ($tstart, $tstop) = get_tstart_tstop( $f, $gti );

my ($bstart, $bstop) = bin_times( $tstart, $tstop, $bin );
wcols $bstart, $bstop, $bstop-$bstart;

sub bin_times {
  my ($t1, $t2, $tbin) = @_;

  dies $tbin unless $tbin > 0;
  die "$t1 $t2" if $t1 >= $t2;

  my (@start, @stop);
  my $t = $t1;
  do {
    push @start, $t;
    $t += $tbin;
    $t = $t2 if $t > $t2;
    push @stop, $t;

  } while $t < $t2;

  if (@start>1 and ($stop[-1]-$start[-1]) < 0.5 * $tbin) {
    $stop[-2] = $stop[-1];
    pop(@start) ; pop(@stop);
  }

  return pdl(\@start), pdl(\@stop);
}

sub get_tstart_tstop {
  my ($f, $gti) = @_;

  my $s = 0;
  my $fptr = Astro::FITS::CFITSIO::open_file($f . "[$gti]", Astro::FITS::CFITSIO::READONLY(), $s);

  my $nrows = 0;
  $fptr->get_num_rows($nrows, $s);

  my $tstart = [];
  my $tstop = [];

  $fptr->read_col_dbl( 1, 1, 1, $nrows, 0, $tstart, undef,  $s);
  $fptr->read_col_dbl( 2, 1, 1, $nrows, 0, $tstop, undef,  $s);

  $fptr->close_file($s);

  $s == 0 or die $s;

  return pdl($tstart)->min, pdl($tstop)->max;

}


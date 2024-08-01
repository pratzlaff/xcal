#! /usr/bin/perl

use strict;
use warnings;

use PDL;
use PDL::Graphics::PGPLOT;
use Astro::FITS::CFITSIO;
use Data::Dumper;

my $version = '0.1';

use FindBin;
use Config;
use Carp;

use Getopt::Long;
my %default_opts = (
		    dev => '/xs',
		    binned => 0,
		    errb => 1,
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'binned!', 'errb!',
	   'dev=s', 'rdb=s', 'newgain!',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

@ARGV == 2 or die "Usage: $0 [options] mkn421|pks2155 0|1|2|3|4|5\n";

my ($src, $band) = @ARGV;

$band =~ /^[0-6]$/ or die "invalid band='$band', should be 0|1|2|3|4|5\n";

my ($nx, $ny) = $band == 0 ? (2, 2) : (1, 3);
dev $opts{dev}, $nx, $ny, { hardch => 1.8, hardlw => 2 };

my @hrc = (
	   'hrc',     # Sherpa HRC fits
#	   'xshrc',   # XSPEC HRC fits
	   );
my @instruments = ( qw/ leg meg heg /, @hrc );

my $title;
for ($src) {
  $_ eq 'pks2155' and $title = 'PKS 2155', last;
  $_ eq 'mkn421' and $title = 'Mkn 421', last;
  die $src;
}

my %ytitles = (
	       alpha => '\ga',
	       beta => '\gb',
	       norm => 'norm',
	       eflux => sprintf("%s-%s keV flux (erg cm\\u-2\\d s\\u-1\\d)", @{band_range($band)}),
	       );

# splits observations into related groups of obsids
my ($obsid_groups, $date_strings) = obsid_groups( $src );

# for writing flux data to an RDB file
my %eflux;
for my $det (qw/ hrc leg meg heg /) {
  $eflux{$det} = {
		  time => pdl([]),
		  obsid => pdl([]),
		  pos => {
			  eflux => pdl([]),
			  eflux_low => pdl([]),
			  eflux_high => pdl([]),
			  },
		  neg => {
			  eflux => pdl([]),
			  eflux_low => pdl([]),
			  eflux_high => pdl([]),
			  },
		 };
}

# for each group, make a plot
#for my $i (0..$#{$obsid_groups}) {
for my $i (0) {

  my @obsids = @{$obsid_groups->[$i]};

  my $title = $title . ': ' . substr( $date_strings->[$i], 1, 10 );

  my ($params, $d, $all) = read_fit_params($src, \@obsids, $band);

  my @instruments = keys %$d;

  my %d = %{$d};
  my %all = %{$all};

  use Data::Dumper;
  #print Dumper \%d;

  # convert to piddles
  for my $key (qw/ obsid time /, @$params) {
    my $data = $all{$key};
    $all{$key} = {};
    $all{$key}{good} = pdl( map { defined $_ ? 1 : 0 } @$data );
    $all{$key}{data} = pdl( $data );
  }

  die if $all{time}{good}->sum != $all{time}{good}->nelem;
  die if $all{obsid}{good}->sum != $all{obsid}{good}->nelem;

=begin comment

%d = (
      leg => {
	      obsid => [],
	      time =>  {
			pos => [],
			neg => [],
		       },
	      alpha => {
			pos => {
				good => [],
				data => [],
			       },
			neg => {
				good => [],
				data => [],
				},
		       },
      )

=cut

  for my $inst (@instruments) {
    $d{$inst}{$_} = pdl($d{$inst}{$_}) for qw/ obsid time /;
    for my $param (map { $_, $_.'_low', $_.'_high' } @$params) {
      for my $order (qw/ pos neg /) {
	my $data = $d{$inst}{$param}{$order};
	$d{$inst}{$param}{$order} = {};
	$d{$inst}{$param}{$order}{good} = pdl( [ map { defined $_ ? 1 : 0 } @$data ] );
	$d{$inst}{$param}{$order}{data} = pdl( $data );
      }
    }
  }

  # make time manageable
  my ($time_min, $time_max) = $all{time}{data}->minmax;

  $d{$_}{time_days} = $d{$_}{time} / 86400 for @instruments;
  $d{$_}{time} -= $time_min for @instruments;

  $all{time}{data} -= $time_min;

  # now slightly offset neg/pos order times
  my $tfudge = 0.003 * ($time_max - $time_min);
  for my $inst (@instruments) {
    my $time = $d{$inst}{time};
    $d{$inst}{time} = {};
    $d{$inst}{time}{neg} = $time - $tfudge;
    $d{$inst}{time}{pos} = $time + $tfudge;
  }

  my @params = @$params;

  my @yranges = (
		 [-0.2, 2.4],
		 [0.2, 2.5],
		 [1e-10, 5e-10],
		 );

  for my $pi (0..$#params) {
    my $param = $params[$pi];
    my $yrange = $yranges[$pi];
    my $possym = 6;
    my $negsym = 7;

    my %colors = (
		  leg => 1,
		  meg => 2,
		  heg => 6,
		  hrc => 4,
		  xshrc => 5,
		 );

    for my $inst (@instruments) {

      my $good;

#      print $inst, ' ', $param, "\n";
#      print $d{$inst}{time}{pos}, $d{$inst}{$param}{pos}, $d{$inst}{obsid}, "\n";

      my $ytitle = exists $ytitles{$param} ? $ytitles{$param} : $param;

      $good = $d{$inst}{$param}{pos}{good};
      for my $color (1, $colors{$inst}) {
      points $d{$inst}{time}{pos}->where($good), $d{$inst}{$param}{pos}{data}->where($good),
	{
	 ( $pi==0 or ($nx==2 and $pi==1)) ? (title => $title) : (),
	 ( $pi==2 or $pi==3 ) ? (xtitle => 'Time offset (\fis\fr)') : (),
	 ytitle => $ytitle,
	 xrange => [$all{time}{data}->minmax],
	 yrange => $yrange,
#	 yrange => [$all{$param}{data}->where($all{$param}{good})->minmax],
	 border => 1,
	 color => $color,#$colors{$inst},
	 symbol => $possym,
	};
      hold;
    }

      $good = $d{$inst}{$param}{neg}{good};
      points $d{$inst}{time}{neg}->where($good), $d{$inst}{$param}{neg}{data}->where($good),
	{
	 color => $colors{$inst},
	 symbol => $negsym,
	};

      if ($opts{errb}) {
      for my $order (qw/ pos neg /) {

	my $linestyle = $order eq 'neg' ? 4 : 1;

	my $good;

	my $errb_opts = {
			 color => $colors{$inst},
			 linestyle => $linestyle,
			 term => 0,
			};

	$good = $d{$inst}{$param.'_low'}{$order}{good};
       	errb(
	     $d{$inst}{time}{$order}->where($good),
	     $d{$inst}{$param}{$order}{data}->where($good),
       	     undef, undef,
       	     $d{$inst}{$param}{$order}{data}->where($good)
	     - $d{$inst}{$param.'_low'}{$order}{data}->where($good),
	     undef,
	     $errb_opts,
       	    );

	$good = $d{$inst}{$param.'_high'}{$order}{good};
       	errb(
	     $d{$inst}{time}{$order}->where($good),
	     $d{$inst}{$param}{$order}{data}->where($good),
       	     undef, undef,
	     undef,
       	     $d{$inst}{$param.'_high'}{$order}{data}->where($good)
	     - $d{$inst}{$param}{$order}{data}->where($good),
	     $errb_opts,
       	    );
      }
      }

    } # for @instruments

    release;

  } # for @params


  for my $inst (qw/ leg meg heg hrc /) {
    next unless exists $d{$inst};

    $eflux{$inst}{time} = $eflux{$inst}{time}->append($d{$inst}{time_days});
    $eflux{$inst}{obsid} = $eflux{$inst}{obsid}->append($d{$inst}{obsid});
    for my $order (qw/ pos neg /){

      for my $param (qw/ eflux eflux_low eflux_high/) {
	my $data = ones($d{$inst}{$param}{$order}{data}) * -1;
	my $good = $d{$inst}{$param}{$order}{good};
	(my $tmp = $data->where($good)) .= $d{$inst}{$param}{$order}{data}->where($good);
	$eflux{$inst}{$order}{$param} = $eflux{$inst}{$order}{$param}->append($data);
      }
    }
  }

} # for @$obsid_groups

if (exists $opts{rdb}) {

  my @cols = qw/ obsid days pos pos_low pos_high neg neg_low neg_high /;

  my $header = <<EOP;
# days     - number of days since 1998.0
#
# pos      - +1 order energy flux (ergs/sec/cm^2)
# pos_low  - one-sigma confidence interval, lower bound
# pos_high -    "         "          "      upper   "
#
# neg, neg_low, neg_high are same as above but for -1 order
#
EOP
  $header .= join( "\t", @cols ) . "\n";
  $header .= join( "\t", ('N') x @cols ) . "\n";

  my $fmt = "%d\t%.3f" . "\t%.4g"x6;


  for my $det (qw/ hrc leg meg heg /) {
      wcols $fmt,
	$eflux{$det}{obsid},
	$eflux{$det}{time},
	$eflux{$det}{pos}{eflux},
	$eflux{$det}{pos}{eflux_low},
	$eflux{$det}{pos}{eflux_high},
	$eflux{$det}{neg}{eflux},
	$eflux{$det}{neg}{eflux_low},
	$eflux{$det}{neg}{eflux_high},
	$opts{rdb} . $det . '.rdb',
	{ header => $header };
    }

}


exit;

sub read_fit_params {

  my ($src, $obsids, $band) = @_;

  my (@params, %d, %all); # return data

  for my $obsid (@$obsids) {

    # my ($pha2) = glob("data/$src/$obsid/tg_reprocess/*_pha2.fits") or die $obsid;
    # my $h = Astro::FITS::CFITSIO::fits_read_header($pha2.'[spectrum]');
    # my $tstart = $h->{TSTART};
    # my $tstop = $h->{TSTOP};

    my @globs = $opts{binned} ? ("data/$src/$obsid/tg_reprocess/evt2_??.fits")
      : ("data/$src/$obsid/tg_reprocess/*_evt2.fits");
    $opts{newgain} and unshift(@globs,
			       $opts{binned} ?
			       "data/$src/$obsid/tg_reprocess_newgain/evt2_??.fits"
			       : "data/$src/$obsid/tg_reprocess_newgain/*_evt2.fits");

    my @files;
    for my $glob (@globs) {
      @files = glob($glob) and last;
    }
    @files or die join("\n", @globs);

    for my $evt2 (@files) {

      my ($tstart, $tstop) = first_last_event_times($evt2);
      my $time = 0.5 * ($tstart + $tstop);

      push @{$all{time}}, $time;
      push @{$all{obsid}}, $obsid;

      my @k;

      my $h = Astro::FITS::CFITSIO::fits_read_header($evt2.'[events]');
      if ($h->{DETNAM} =~ /hrc/i) { @k = @hrc; }
      elsif ($h->{GRATING} =~ /hetg/i) {
	@k = qw/ meg heg /;
      } elsif ($h->{GRATING} =~ /letg/i) {
	@k = qw/ leg /;
      } else {
	die $h->{DETNAM} . ' / ' . $h->{GRATING} ;
      }

      my $date = substr($h->{'DATE-OBS'}, 0, 7);

      for my $k (@k) {

	my $psub = $k =~ /xs/ ? \&xspec_params : \&sherpa_params;

	push @{$d{$k}{time}}, $time;
	push @{$d{$k}{obsid}}, $obsid;
	push @{$d{$k}{date}}, $date;

	for my $o (qw/ neg pos /) {

	  my $dir = "data/$src/$obsid/fit/B$band";

	  $dir = "data/$src/$obsid/fit_newgain/B$band" if $k eq 'leg' and $opts{newgain};

	  $dir .= '/binned' if $opts{binned};

	  my $f =
	    $k =~ /xs/ ? "$dir/xsLEG_${o}.log" :
	      $k eq 'hrc' ? "$dir/LEG_${o}" :
		$dir . '/' .uc($k). "_${o}";

	  if ($opts{binned}) {
	    my ($n) = $evt2 =~ /(_\d{2})\.fits/ or die $evt2;
	    $f =~ s/_$o/_$o$n/;
	  }

	  my ($eflux, $names, $best, $low, $high) = $psub->($f) or next;
	  if (! @params) {
	    @params = @$names;
	  }

	  for my $i (0..$#{$names}) {
	    my $name = $names->[$i];
	    my $best = $best->[$i];
	    my $low = $low->[$i];
	    my $high = $high->[$i];

	    $low = $high = $best unless $opts{errb};

	    warn "$obsid, $name" unless defined $best;

	    push @{$d{$k}{$name}{$o}}, $best;
	    push @{$d{$k}{$name.'_low'}{$o}}, $low;
	    push @{$d{$k}{$name.'_high'}{$o}}, $high;

	    push @{$all{$name}}, $best, $low, $high;

	    if ($name eq 'norm') {
	      my ($eflux_low, $eflux_high);
	      if (defined $best) {
		$eflux_low = $eflux*$low/$best if defined $low;
		$eflux_high = $eflux*$high/$best if defined $high;
	      }

	      push @{$d{$k}{eflux}{$o}}, $eflux;
	      push @{$d{$k}{eflux_low}{$o}}, $eflux_low;
	      push @{$d{$k}{eflux_high}{$o}}, $eflux_high;
	      push @{$all{eflux}}, $eflux, $eflux_low, $eflux_high;
	      warn "$obsid, eflux" unless defined $eflux;

	      push @params, 'eflux' if @params = @$names;

	    }

	  } # for names

	} # for order

      } # for key

    } # for @files

  } # for @obsids

  return \(@params, %d, %all);
}

sub obsid_groups {

  my $src = shift;

  my $obsid_file = "obsids_$src";
  open(my $fh, '<', $obsid_file) or die "could not open '$obsid_file': $!";

  my (@obsids, @times, @date_obs);

  while (<$fh>) {
    /^\d+/ or next;
    my ($obsid) = split;

    # FIXME
#    $obsid eq '7291' and next;
#    $obsid eq '8380' and next;
#    $obsid eq '9705' and next;
#    $obsid eq '9710' and next;
#    $obsid eq '11974' and next;

    my ($pha2) = glob("data/$src/$obsid/tg_reprocess/*_pha2.fits") or die;
    my $h = Astro::FITS::CFITSIO::fits_read_header($pha2.'[spectrum]');
    my $tstart = $h->{TSTART};
    my $date_obs = $h->{'DATE-OBS'};

    push @obsids, $obsid;
    push @times, $tstart;
    push @date_obs, $date_obs;
  }

  my $times = double \@times;

  my $sorti = $times->qsorti;
  $_ = $_->index($sorti)->sever for $times;

  @obsids = @obsids[$sorti->list];
  @times = @times[$sorti->list];
  @date_obs = @date_obs[$sorti->list];

  my $tdiff = $times->slice('1:-1') - $times->slice('0:-2');
  my $spliti = which( $tdiff > 7*86400 );

  my @obsid_groups;
  my @date_strings;

  my $lasti = 0;
  for my $i ($spliti->list, $#obsids) {
    push @obsid_groups, [ @obsids[$lasti..$i] ];
    push @date_strings, $date_obs[$lasti];
    $lasti = $i+1;
  }

  return \( @obsid_groups, @date_strings );

}

sub parse_xspec_conf {
  my $line = shift;
  my @tmp = split ' ', $line;
  my ($low, $high) = @tmp[2,3];
  my ($blow, $bhigh) = $tmp[4] =~ /\((.*),(.*)\)/;
  return $low-$blow, $low, $high;
}

sub xspec_params {
  die "this needs to be completely re-written, see sherpa_params() below";
  my $f = shift;

  my ($a, $alo, $ahi, $b, $blo, $bhi, $n, $nlo, $nhi);
  my ($pflux, $eflux);

  open(my $fh, '<', $f) or warn("cannot open '$f' for reading"), return;
  while (defined (my $line = <$fh>)) {
    $line =~ /^#\s+parameter\s+confidence range/i or next;
    ($a, $alo, $ahi) = parse_xspec_conf(scalar <$fh>);
    ($b, $blo, $bhi) = parse_xspec_conf(scalar <$fh>);
    ($n, $nlo, $nhi) = parse_xspec_conf(scalar <$fh>);
    last;
  }

  while (defined (my $line = <$fh>)) {
    $line =~ /^#\s+model flux/i or next;
    ($pflux, $eflux) = $line =~ /^#\s+model flux\s+(\S+)\s+photons\s+\((\S+)\s+/i;
    last;
  }

  defined $a and defined $pflux or die;

#  return $a, $b, $n;
  return $a, $alo, $ahi, $b, $blo, $bhi, $n, $nlo, $nhi, $eflux;
}

# uses the base name and reads '.{conf,flux}' files
sub sherpa_params {
  my $base = shift;
  my $fh;

  my $fluxf = $base . '.flux';
  open($fh, '<', $fluxf) or warn("cannot open '$fluxf' for reading"), return;
  my $eflux = ( split( ' ', scalar(<$fh>) ) )[3];
  close $fh;

  my $conff = $base . '.conf';

  open($fh, '<', $conff) or warn("cannot open '$conff' for reading"), return;

  my (@names, @best, @low, @high);

  while (defined( my $line = <$fh>)) {
    next if $line =~ /^warning/i;

    if ($line =~ /\s+Param/) {
      <$fh>; # skip a line
      while (defined(my $line = <$fh>)) {

	next unless $line =~ /\w/; # there are two empty lines at the end of the file

	my ($name, $best, $low, $high) = split ' ', $line;

	# undetermined values show up as '----'
	for my $val ($best, $low, $high) { $val = undef unless $val =~ /\d/ }

	$name =~ s/^\w+\.//; # e.g., turn lp.alpha into alpha

	if (defined $best) {
	  $low = $best + $low if defined $low;
	  $high = $best + $high if defined $high;
	}

	push @names, $name;
	push @best, $best;
	push @low, $low;
	push @high, $high;
      }
    }
  }

  return ($eflux, \@names, \@best, \@low, \@high);
}

sub first_last_event_times {
  my $file = shift;

  my $fptr = Astro::FITS::CFITSIO::open_file($file, Astro::FITS::CFITSIO::READONLY(), my $status = 0);

  $fptr->movnam_hdu(Astro::FITS::CFITSIO::BINARY_TBL(), 'events', 0, $status);

  my $time_col;
  $fptr->get_colnum(Astro::FITS::CFITSIO::CASEINSEN(), 'time', $time_col, $status);

  my $nrows;
  $fptr->get_num_rows($nrows, $status);

  my ($first, $last);
  $fptr->read_col_dbl($time_col, 1, 1, 1, 0, $first, undef, $status);
  $fptr->read_col_dbl($time_col, $nrows, 1, 1, 0, $last, undef, $status);

  $fptr->close_file($status);

  $status and die $file;

  return $first->[0], $last->[0];
}

sub band_range {
  my $band = shift;
  my %bands = (
	       0 => [0.5, 8.0],
	       1 => [0.33, 0.54],
	       2 => [0.54, 0.85],
	       3 => [0.85, 1.5],
	       4 => [1.5, 4],
	       5 => [4, 10],
	       6 => [0.26, 0.514],
	       );
  exists $bands{$band} or die $band;
  return $bands{$band};
}

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
}

=head1 NAME

template - A template for Perl programs.

=head1 SYNOPSIS

cp template newprog

=head1 DESCRIPTION

blah blah blah

=head1 OPTIONS

=over 4

=item --help

Show help and exit.

=item --version

Show version and exit.

=back

=head1 AUTHOR

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> Feb 2014

=head1 SEE ALSO

perl(1).

=cut


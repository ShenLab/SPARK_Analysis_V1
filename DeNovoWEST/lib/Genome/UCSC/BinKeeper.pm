package Genome::UCSC::BinKeeper;

use strict;
use warnings;
use Carp;
use Storable qw(dclone);
use Iterator::Simple qw|iterator|;
use IO::File;
use IO::Detect qw|is_filehandle|;
use Scalar::Util qw|reftype|;
use List::Util qw(max min sum);
use Data::Dumper;
use Iterator::DBI;
use Genome::Ranges qw|validate_rng validate_elem range_overlap|;
use Utils::Hash qw|merge_opts|;
use Utils::File qw|open_file|;
use Exporter;

use base qw(Exporter);
our @EXPORT_OK = qw( binLevelsExtended
				  binLevels
				  binFirstShift
				  binNextShift
				  binOffsetOldToExtended
				  BINRANGE_MAXEND_512M
				  binOffsetExtended
				  binOffset
				  binFromRangeStandard
				  binFromRangeExtended
				  binFromRange
				  binFromRangeBinKeeperExtended
			   );
our %EXPORT_TAGS = ( 'all' => \@EXPORT_OK );


=head1 NAME

Genome::UCSC::BinKeeper - UCSC's Binning Scheme for Rapid Range Query.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 NOTE

UCSC's genome intervals are 0-based half-open.
To conform with other inhouse modules, we use 1-based inclusive intervals.

=head1 DESCRIPTION

UCSC's BinKeeper approach allow us to fast pin-point our interest
into the given genomic ranges from a set of possibly overlapping
genomic intervals. Compare with the perl's Array::IntSpan techniques,
it abandons the assumption that ranges must not overlap with each
other.

The underlying core of this module is UCSC's notion of bin, which is
a small integer index calculated from the start and end. Each interval
will be stored into a bin slot.

The object is internally implemented as hash ref of bins. Each range
is stored as array reference with: chromStart, chromEnd, and optionally
value fields. We also retain the interface of Genome::Ranges class.


=head1 SYNOPSIS

    use Genome::UCSC::BinKeeper;

    my $bk = Genome::UCSC::BinKeeper->new();
    while(<STDIN>) {
        ($chr, $st, $ed, $score) = split;
        $bk->add($chr, $st, $ed, $score); 
    }
    my @overlaps = $bk->find_range("chr1", 12345, 34567);
    # ...

=head1 EXPORT

Thoes exported functions are only for internal use.

=over 5 

=item * binLevelsExtended

=item * binLevels

=item * binFirstShift

=item * binNextShift

=item * binOffsetOldToExtended

=item * BINRANGE_MAXEND_512M

=item * binOffsetExtended

=item * binOffset

=item * binFromRangeStandard

=item * binFromRangeExtended

=item * binFromRange

=item * binFromRangeBinKeeperExtended

=back

=cut

{
  # Modified from UCSC's binRange.h and binRange.c
  ################### ************** IMPORTANT ********************#######################
  # Note that UCSC's intervals are zero-based half-open, while the coordinations systems
  # in my inhouse perl modules are one-based inclusive.
  # We will make use UCSC's convention to store intervals internally within this module
  # but all the interfaces will use 1-based intervals with both ends closed
  ################### ************** IMPORTANT ********************#######################

  # binRange Stuff to handle binning - which helps us restrict
  # our attention to the parts of database that contain info
  # about a particular window on a chromosome. This scheme
  # will work without modification for chromosome sizes up
  # to half a gigaBase.  The finest sized bin is 128k (1<<17).
  # The next coarsest is 8x as big (1<<13).  There's a hierarchy
  # of bins with the chromosome itself being the final bin.
  # Features are put in the finest bin they'll fit in.

  my @binOffsetsExtended = (4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0);
  my @binOffsets = (512+64+8+1, 64+8+1, 8+1, 1, 0);
  my $_binFirstShift = 17;
  my $_binNextShift = 3;
  my $_binOffsetOldToExtended = 4681;
  my $BINRANGE_MAXEND_512M = (512*1024*1024);

  sub binLevelsExtended { return scalar(@binOffsetsExtended); }
  sub binLevels     { return scalar(@binOffsets); }
  sub binFirstShift { return $_binFirstShift; }
  sub binNextShift  { return $_binNextShift;  }
  sub binOffsetOldToExtended { return $_binOffsetOldToExtended; }
  sub BINRANGE_MAXEND_512M { return $BINRANGE_MAXEND_512M; }

  sub binOffsetExtended {
	my ($level) = @_;
	return $binOffsetsExtended[$level];
  }

  sub binOffset {
	my ($level) = @_;
	return $binOffsets[$level];
  }

  sub binFromRangeStandard {
	my ($start, $end) = @_;
	# Given start,end in chromosome coordinates assign it
	# a bin.   There's a bin for each 128k segment, for each
	# 1M segment, for each 8M segment, for each 64M segment,
	# and for each chromosome (which is assumed to be less than
	# 512M.)  A range goes into the smallest bin it will fit in.
	#
	my ($startBin, $endBin) = ($start, $end-1);
	$startBin >>= $_binFirstShift;
	$endBin >>= $_binFirstShift;
	for ( my $i = 0; $i< @binOffsets; ++$i) {
	  if ($startBin == $endBin) {
		return $binOffsets[$i] + $startBin;
	  }
	  $startBin >>= $_binNextShift;
	  $endBin >>= $_binNextShift;
	}
	croak "start $start, end $end out of range in findBin (max is 512M)";
	return 0;
  }

  sub binFromRangeExtended {
	my ($start, $end) = @_;
	# Given start,end in chromosome coordinates assign it
	# a bin.   There's a bin for each 128k segment, for each
	# 1M segment, for each 8M segment, for each 64M segment,
	# for each 512M segment, and one top level bin for 4Gb.
	#      Note, since start and end are int's, the practical limit
	#      is up to 2Gb-1, and thus, only four result bins on the second
	#      level.
	# A range goes into the smallest bin it will fit in.
	my ($startBin, $endBin) = ($start, $end-1);
	$startBin >>= $_binFirstShift;
	$endBin >>= $_binFirstShift;
	for ( my $i = 0; $i< @binOffsetsExtended; ++$i) {
	  if ($startBin == $endBin) {
		return $_binOffsetOldToExtended + $binOffsetsExtended[$i] + $startBin;
	  }
	  $startBin >>= $_binNextShift;
	  $endBin >>= $_binNextShift;
	}
	croak "start $start, end $end out of range in findBin (max is 2Gb)";
	return 0;
  }

  sub binFromRange {
	my ($start, $end) = @_;
	if ($end <= $BINRANGE_MAXEND_512M) {
	  return binFromRangeStandard($start, $end);
	}
	else {
	  return binFromRangeExtended($start, $end);
	}
  }

  sub binFromRangeBinKeeperExtended {
	my ($start, $end) = @_;
	# This is just like binFromRangeExtended() above, but it doesn't limit
	# the answers to the range from _binOffsetOldToExtended and up.
	# It simply uses the whole new bin scheme as if it was the only one.
	my ($startBin, $endBin) = ($start, $end-1);
	$startBin >>= $_binFirstShift;
	$endBin >>= $_binFirstShift;
	for ( my $i = 0; $i< @binOffsetsExtended; ++$i) {
	  if ($startBin == $endBin) {
		return $binOffsetsExtended[$i] + $startBin;
	  }
	  $startBin >>= $_binNextShift;
	  $endBin >>= $_binNextShift;
	}
	croak "start $start, end $end out of range in findBin (max is 2Gb)";
	return 0;
  }

}


=head1 SUBROUTINES/METHODS

=head2 $class->new [RANGE]

Create a new empty object or new object from ranges.

Also allow construction from file, assuming first three columns will be
chrom, start, and end; and the remaining will be values. Users can optionally
provide a call back to construct value fields. The input for the callback
will be the content from column 4 and later.

=cut

sub new {
	my ($class, $ranges, $argref) = @_;
	my $arg = merge_opts($argref, callback => undef, bed => undef);
	if (defined $ranges) {
		if (reftype $ranges && reftype $ranges eq 'HASH') {
			validate_rng($ranges);
			my %bk;
  			while (my ($chrom, $intvs) = each %$ranges) {
				foreach my $intv (@$intvs) {
					add(\%bk, $chrom, @$intv);
				}
  			}
  			return bless \%bk, $class;
		}
		else {
			my $fin;
			my %bk;
			if (is_filehandle($ranges)) {
				$fin = $ranges;
			}
			elsif (-f $ranges) {
				$fin = open_file($ranges);
			}
			else {
				croak "Cannot determine the type of $ranges";
			}
			while(<$fin>) {
				my @row = split;
				$row[1] += 1 if $arg->{bed};
				validate_elem(@row[0,1,2]);
				my $chrom = shift @row;
				if (defined $arg->{callback}) {
					my $start = shift @row;
					my $end = shift @row;
					add(\%bk, $chrom, $start, $end, $arg->{callback}->(@row));
				}
				else {
					add(\%bk, $chrom, @row);
				}
			}
			return bless \%bk, $class;
		}
	}
	else {
		return bless {}, $class;
	}
}

=head2 $obj->add CHR, START, END [,VAL]

Add an item into the binKeeper, by appending it to the appropriate bin slot.

=cut

sub add {
	my ($self, $chrom, $start, $end, @vals) = @_;
	$end = $start unless defined $end;
	$start = 1 if $start < 1;
	if (validate_elem($chrom, $start, $end)) {
 		# 1-based both-closed interval to 0-based half-open conversion
 		# e.g. [11,14] should become [10,14)
 		my $bin = binFromRangeBinKeeperExtended($start-1, $end);
 		push @{$self->{$chrom}{$bin}}, [$start, $end, @vals];
 		return $self;
 	}
}

=head2 $self->write FILE [, CALLBACK]

This method has the same interface as L<Genome::Ranges>

=cut

sub write {
	my $self = shift @_;
	$self->Genome::Ranges::write(@_);
}



=head2 $self->iter

Return an iterator. The value returned by calling iterator is
[chrom, start, end, val]. 

=cut

sub iter {
  my ($self) = @_;

  my @chroms = sort keys %$self;
  return undef unless @chroms;

  my $chr = 0;
  my $chrom = $chroms[$chr];
  # bin will be sorted numerically
  my @bins  = sort {$a<=>$b} keys %{$self->{$chrom}};
  my $bi  = 0;
  my $bin   = $bins[$bi];
  my $index = 0;
  my $ii = 0;

  return iterator {
		# look for the next defined interval
		if ( !defined $self->{$chrom}{$bin}[$index] ) {
		  while ( ++$bi < @bins ) {
			$bin = $bins[$bi];
			$index = 0;
			if ( defined $self->{$chrom}{$bin}[$index] ) {
				if (wantarray) {
					return ($chrom, @{$self->{$chrom}{$bin}[$index++]});
				}
				else {
					return [$chrom, @{$self->{$chrom}{$bin}[$index++]}];
				}
			}
		  }
		  while ( ++$chr < @chroms ) {
			$chrom = $chroms[$chr];
			$bi = 0;
			@bins  = sort {$a<=>$b} keys %{$self->{$chrom}};
			$bin = $bins[$bi];
			$index = 0;
			$ii = 0;
			if ( defined $self->{$chrom}{$bin}[$index] ) {
				if (wantarray) {
					return ($chrom, @{$self->{$chrom}{$bin}[$index++]});
				}
				else {
					return [$chrom, @{$self->{$chrom}{$bin}[$index++]}];
				} 
			}
		  }
		  return;
		}
		else {
			if (wantarray) {
				return ($chrom, @{$self->{$chrom}{$bin}[$index++]});
			}
		  	else {
		  		return [$chrom, @{$self->{$chrom}{$bin}[$index++]}];
		  	}
		}
	 };
}

=head2 $self->to_ranges

 Transform to Genome::Ranges object.

=cut

sub to_ranges {
  my $self = shift;
  my %ranges;
  my $iter = $self->iter();
  while (my $dat = $iter->()) {
	my ($chrom, @intv) =@$dat;
	push @{$ranges{$chrom}}, [@intv];
  }
  return bless \%ranges, "Genome::Ranges";
}

=head2 $self->count [CHROM]

Return the number of intervals.

=cut

sub count {
	my ($self, $chrom) = @_;
	my @chrom;
	if (defined $chrom) {
		push @chrom, $chrom;
	}
	else {
		@chrom = keys %$self;
	}
	my $count = 0;
	foreach my $k (@chrom) {
		foreach my $b (keys %{$self->{$k}}) {
			next unless defined $self->{$k}{$b};
			$count += @{$self->{$k}{$b}};
		}
	}
	return $count;
}

=head2 $self->size [CHROM]

Return the size of intervals.

=cut

sub size {
	my ($self, $chrom) = @_;
	my @chrom;
	if (defined $chrom) {
		push @chrom, $chrom;
	}
	else {
		@chrom = keys %$self;
	}
	my $size = 0;
	foreach my $k (@chrom) {
		foreach my $b (keys %{$self->{$k}}) {
			next unless defined $self->{$k}{$b};
			$size += sum(map { $_->[1]-$_->[0]+1 } @{$self->{$k}{$b}});
		}
	}
	return $size;
}


=head2 $self->find_range CHROM, START, END [, OPTIONS]

Return an array of all intervals that overlap with a given query interval. 

=head3 Options

=over 5

=item qCover / tCover

Coverage threshold for query and target.

=item calc

Also also calculate the overlapping proportion of.
It will be appended to the value fields as a hashref with keys:
pT, pQ, and iLen.

=back

=cut

sub find_range {
	my ($self, $chrom, $qStart, $qEnd, $argref) = @_;

	$qEnd = $qStart unless defined $qEnd;
 	my $arg = merge_opts($argref, tCover => 1e-10, qCover => 1e-10, calc => undef);

	validate_elem($chrom, $qStart, $qEnd);

	return unless defined $self->{$chrom};
	
	my @overlaps;
	$self->_pinpoint_range($chrom, $qStart, $qEnd, sub {
		my ($interval, $aux) = @_;
		my $tLen = $interval->[1]-$interval->[0]+1;
		my $qLen = $qEnd-$qStart+1;
		my $iLen = range_overlap($interval->[0],$interval->[1],$qStart,$qEnd);
		if ( $iLen >= $arg->{tCover} * $tLen &&
			 $iLen >= $arg->{qCover} * $qLen
			) {
			push @$aux, $interval;
		}
		}, \@overlaps);

  	# Append Q and T overlap as hash reference to each overlapping interval
  	if ($arg->{calc}) {
  		my @overlapcalc;
  		foreach my $intv (@overlaps) {
  			my ($tStart, $tEnd) = ($intv->[0], $intv->[1]);
  			my $iLen = range_overlap($tStart, $tEnd, $qStart, $qEnd);
  			my %stat;
  			$stat{T} = $iLen/($tEnd-$tStart+1);
  			$stat{Q} = $iLen/($qEnd-$qStart+1);
  			$stat{TxQ} =  $stat{T}*$stat{Q};
  			$stat{iLen} = $iLen;
  			my $intvcln = dclone($intv);
  			push @{$intvcln} => \%stat;
  			push @overlapcalc, $intvcln;
  		}
  		return @overlapcalc;
  	}
  	else {
  		return @overlaps;
  	}
}

# Helper function to remove code duplication
# $action is a code ref, whose first arg is the interval element
# the second arg is user provided.
sub _pinpoint_range {
	my ($self, $chrom, $start, $end, $action, $args) = @_;
  	# Note the transformation of coordination system
  	$start = 1 if $start < 1;
  	$end = 1 if $end < 1;
  	my $startBin = ($start-1) >> binFirstShift;
  	my $endBin   = ($end-1) >> binFirstShift;
  	for ( my $ii = 0; $ii < binLevelsExtended; $ii ++ ) {
  		my $offset = binOffsetExtended($ii);
  		for ( my $jj = $startBin + $offset; $jj <= $endBin + $offset; $jj ++ ) {
  			next unless defined $self->{$chrom} && defined $self->{$chrom}{$jj};
  			foreach ( @{$self->{$chrom}{$jj}} ) {
  				if ( range_overlap($_->[0],$_->[1],$start,$end) > 0 ) {
  					$action->($_, $args);
  				}
  			}
  		}
  		$startBin >>= binNextShift;
  		$endBin   >>= binNextShift;
  	}
  	return 0;
}


=head2 $self->any_overlap CHROM, START, END

Return TRUE if start/end overlaps with any item in binKeeper.

=cut 

sub any_overlap {
	my ($self, $chrom, $start, $end) = @_;
	$end = $start unless defined $end;
	validate_elem($chrom, $start, $end);

	return 0 unless defined $self->{$chrom};
	my $startBin = ($start-1) >> binFirstShift;
	my $endBin   = ($end-1) >> binFirstShift;
	for ( my $i = 0; $i < binLevelsExtended; ++$i ) {
		my $offset = binOffsetExtended($i);
		for ( my $j = $startBin + $offset; $j <= $endBin + $offset; ++$j ) {
			next unless defined $self->{$chrom}{$j};
			foreach ( @{$self->{$chrom}{$j}} ) {
				if ( range_overlap($_->[0],$_->[1],$start,$end) > 0 ) {
					return 1;
				}
			}
		}
		$startBin >>= binNextShift;
		$endBin   >>= binNextShift;
	}
	return 0;
}


=head1 AUTHOR

Xueya Zhou, C<< <xueyazhou at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-genome at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Genome>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Genome::UCSC::BinKeeper


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Genome>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Genome>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Genome>

=item * Search CPAN

L<http://search.cpan.org/dist/Genome/>

=back



=head1 LICENSE AND COPYRIGHT

Copyright 2018 Xueya Zhou.

This program is free software; you can redistribute it and/or modify it
under the terms of the the Artistic License (1.0). You may obtain a
copy of the full license at:

L<http://www.perlfoundation.org/artistic_license_1_0>

Aggregation of this Package with a commercial distribution is always
permitted provided that the use of this Package is embedded; that is,
when no overt attempt is made to make this Package's interfaces visible
to the end user of the commercial distribution. Such use shall not be
construed as a distribution of this Package.

The name of the Copyright Holder may not be used to endorse or promote
products derived from this software without specific prior written
permission.

THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.


=cut

1; # End of Genome::UCSC::BinKeeper

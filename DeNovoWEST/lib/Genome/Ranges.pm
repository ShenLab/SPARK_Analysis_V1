package Genome::Ranges;

use strict;
use warnings;
use Carp;
use List::Util qw|min max sum|;
use Scalar::Util qw|reftype|;
use Storable qw(dclone);
use IO::Detect;
use Sort::Versions;
use Iterator::Simple qw(iterator);
use Utils::File qw|open_file|;
use Utils::Hash qw|merge_opts|;
use Exporter;

our @EXPORT_OK = qw(validate_elem validate_rng
					range_to_spec spec_to_range range_overlap);
use base qw(Exporter);


=head1 NAME

Genome::Ranges - General class for genomic intervals.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 DESCRIPTION

This is a simple object-oriented interface to genomic intervals.

B<NOTE>: This inhouse class, genomic intervals are 1-based, inclusive.
It differs from UCSC, where intervals are half-open, 0-based.

The object is internally represented as hash ref, using chrom as key. Ranges on each
chromosome are stored as array reference. Each interval is an array ref with >=2 
fields: [chromStart, chromEnd, (val, ...)]. The value part is optional, and can flexibly 
stores additional data associated with the interval.


=head1 SYNOPSIS

	use Genome::Ranges;

	$rng = Genome::Ranges->new();
	while(<FIN>) {
		my $spec = (split)[1];
		my ($chr, $start, $end) = split;
		$rng->add($chr, $start, $end, {ind => $count ++})
	}

	$rng->sort;
	$it = $rng->iter();
	while($dat = <$it>) {
		print join("\t", @{$dat}[0,1,2]), "\n";
	}

=head1 EXPORT

=over 5

=item * validate_elem

=item * validate_rng

=item * range_to_spec 

=item * spec_to_range 

=item * range_overlap

=back

=head1 SUBROUTINES/METHODS

=head2 $class->new [RANGES]

Create an empty new object or from a hashref.

Also allow construction from file, assuming first three columns will be
chrom, start, and end; and the remaining will be values. Users can optionally
provide a call back to construct value fields. The input for the callback
will be the content from column 4 and later.

Options when input ranges are read from file:
	* callback - callback subroutine to collect data starting from the 4th column
	* bed - indicating the input file follows BED format, i.e. intervals are half-open, zero-based.

=cut


sub new {
	my ($class, $ranges, $argref) = @_;
	my $arg = merge_opts($argref, callback => undef, bed => undef);
	if (defined $ranges) {
		if (reftype $ranges && reftype $ranges eq 'HASH') {
			return bless validate_rng($ranges), $class;
		}
		else {
			my $fin;
			if (is_filehandle($ranges)) {
				$fin = $ranges;
			}
			elsif (-f $ranges) {
				$fin = open_file($ranges);
			}
			else {
				croak "Cannot determine the type of $ranges";
			}
			my %rngs;
			while(<$fin>) {
				my @row = split;
				$row[1] += 1 if $arg->{bed};
				validate_elem(@row[0,1,2]);
				my $chrom = shift @row;
				if (defined $arg->{callback}) {
					my $start = shift @row;
					my $end = shift @row;
					push @{$rngs{$chrom}}, [$start, $end, $arg->{callback}->(@row)];
				}
				else {
					push @{$rngs{$chrom}}, \@row;
				}
			}
			return bless \%rngs, $class;
		}
	}
	else {
		return bless {}, $class;
	}
}


=head2 $self->add CHROM, START, END [,VAL]

Add element into ranges. Value part is optional, can be an array.

=cut

sub add {
	my ($self, $chrom, $start, $end, @vals) = @_;
	$end = $start unless defined $end;
	validate_elem($chrom, $start, $end);
	push @{$self->{$chrom}}, [$start, $end, @vals];
	return $self;
}

=head2 $self->write FILE, ARGREF

Write genomic intervals to file.

Options:
	* callback -- Users can provide callback to format output string of values.
	* bed -- write intervals in BED format. Default is to output value fields as is. 

=cut

sub write {
	my ($self, $file, $argref) = @_;
	my $arg = merge_opts($argref, callback => undef, bed => undef);
	my $fout;
	if (is_filehandle ($file)) {
		$fout = $file;
	}
	else {
		$fout = open_file($file, {write => 1});
	}
	my $it = $self->iter();
	while(my $dat = $it->()) {
		$dat->[1] -= 1 if defined $arg->{bed};
		if (defined $arg->{callback}) {
			my $chrom = shift @$dat;
			my $start = shift @$dat;
			my $end = shift @$dat;
			print $fout join("\t", $chrom, $start, $end, $arg->{callback}->(@$dat)), "\n";
		}
		else {
			print $fout join("\t", @$dat), "\n";
		}
	}
	return $self;
}


=head2 sort [cb_score]

Sort intervals in each chromosome. It is in-situ sorting.

Default is to sort based on the ascernding order of chromStart. Users can 
provide optional call back scoring function. The input for the callback is
interval element [st, ed, val]

=cut

sub sort {
	my ($self, $score) = @_;
	$score = sub { $_[0][0] } unless defined $score; 
	foreach my $chrom (sort { versioncmp($a, $b) } keys %$self) {
		my @sorted_intervals =
			map { $_->[1] }
				sort { $a->[0] <=> $b->[0] }
					map { [ $score->($_), $_ ] }
						@{$self->{$chrom}};
		$self->{$chrom} = [@sorted_intervals];
	}
	return $self;
}

=head2 $self->clone

Make a deep copy of itself.

=cut

sub clone {  return dclone $_[0]; }

=head2 $self->iterator

Return an iterator. The value returned by calling iterator is unpacked
[chrom, start, end, val].

=cut

sub iter {
	my ($self) = @_;

	my @chroms = sort { versioncmp($a, $b) } keys %$self;
	my $chr = 0;
	my $index = 0;
	my $chrom = $chroms[$chr];

  	return iterator {
  		if ( !defined $self->{$chrom}[$index] ) {
			while ( ++$chr < @chroms ) {
				$chrom = $chroms[$chr];
				$index = 0;
				if ( defined $self->{$chrom}[$index] ) {
					if (wantarray) {
						return ($chrom, @{$self->{$chrom}[$index++]});
					}
					else {
						return [$chrom, @{$self->{$chrom}[$index++]}];
					}
				}
		  	}
		  	return;
		}
		else {
			if (wantarray) {
				return ($chrom, @{$self->{$chrom}[$index++]});
			}
			else {
				return [$chrom, @{$self->{$chrom}[$index++]}];
			}
		}
  	};
}

=head2 $self->count [CHROM]

Get the total number of intervals.

=cut

sub count {
	my ($self, $chrom) = @_;
	if (defined $chrom) {
		return defined $self->{$chrom} ? scalar(@{$self->{$chrom}}) : 0;
	}
	else {
		my $count = 0;
		foreach my $k (keys %$self) {
			next unless defined $self->{$k};
			$count += @{$self->{$k}};
		}
		return $count;
	}
}

=head2 $self->size [CHROM]

Get the total size of intervals.

=cut 

sub size {
	my ($self, $chrom) = @_;
	if (defined $chrom) {
		return defined $self->{$chrom} ? sum(map { $_->[1]-$_->[0]+1 } @{$self->{$chrom}}) : 0;
	}
	else {
		my $count = 0;
		foreach my $k (keys %$self) {
			next unless defined $self->{$k};
			$count += sum(map { $_->[1]-$_->[0]+1 } @{$self->{$k}});
		}
		return $count;
	}
}

=head2 Non-method subroutines


=head3 validate_elem CHROM, START [,END]

Validate an interval.

=cut

sub validate_elem {
	my ($chrom, $start, $end) = @_;
	my $flag = 1;
	unless ( $chrom =~ /^\w\S*$/ && $start =~ /^\d+$/ && $end =~ /^\d+$/ ) {
		$flag = 0;
	}
	else {
		if (defined $end) {
			$flag = 0 unless $start <= $end;
		}
	}
	unless ($flag) {
		croak "Incorrect genome range: $chrom:$start-$end";
	}
	else {
		return $flag;
	}
}

=head3 validate_rng RANGE

Validate a range data structure.

=cut

sub validate_rng {
	my ($self) = @_;
	croak "Incorrect reference type" unless reftype $self eq 'HASH';
	while(my ($chrom, $ranges) = each %$self) {
		croak "Incorrect chromosome name: $chrom" unless $chrom =~ /^\w\S*$/;
		croak "Ranges should be found in arrayref" unless reftype $ranges eq 'ARRAY';
		foreach my $intv (@$ranges) {
			unless ($intv->[0] =~ /^\d+$/ && $intv->[1] =~ /^\d+$/ &&
				$intv->[1] >= $intv->[0] ) {
	  			croak "Invalid interval -- $chrom:$intv->[0]-$intv->[1]";
			}
		}
	}
	return $self;
}


=head3 range_to_spec/spec_to_range

Convertion between range and position specification

=cut

sub range_to_spec {
	my ($chrom, $start, $end) = @_;
	#croak "Incorrect genome range: $chrom:$start-$end"
	#	unless $chrom =~ /^\w+$/ && $start =~ /^[\d]+$/
	#		&& (!defined $end || $end =~ /^[\d]+$/);
	validate_elem($chrom, $start, $end);
	if (!defined $end) {
		return sprintf("%s:%d",$chrom,$start);
	}
	else {
		croak "$chrom: $start > $end" if $start > $end;
		return sprintf("%s:%d-%d", $chrom,$start,$end);	
	}
}

sub spec_to_range {
	my ($spec) = @_;
	if ($spec =~ /^(\w\S*):([\d,]+)\-([\d,]+)$/) {
		my ($chrom, $chromStart, $chromEnd) = ($1,$2,$3);
		$chromStart =~ s/,//g; $chromEnd =~ s/,//g;
		croak "$chrom : $chromStart > $chromEnd" if $chromStart > $chromEnd;
		return ($chrom,$chromStart,$chromEnd);
	}
	elsif ($spec =~ /^(\w\S*):([\d,]+)$/) {
		my ($chrom,$position) = ($1,$2,$3);
		$position =~ s/,//g;
		return ($chrom,$position);
	}
	else {
		croak "Incorrect spec format: $spec";
	}
}

=head3 range_overlap START1, END1, START2, END2

Return amount of bases two ranges overlap with each over, 0 if no overlap.

=cut

sub range_overlap {
  my ($start1, $end1, $start2, $end2) = @_;
  my $s = max($start1, $start2);
  my $e = min($end1, $end2);
  my $o = $e-$s+1; # note that our intervals are inclusive
  return $o < 0 ? 0 : $o;
}



=head1 AUTHOR

Xueya Zhou, C<< <xueyazhou at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-genome at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Genome>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Genome::Ranges


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

1; # End of Genome

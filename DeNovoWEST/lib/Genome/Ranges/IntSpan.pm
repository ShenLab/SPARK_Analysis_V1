package Genome::Ranges::IntSpan;
use base qw(Genome::Ranges);

use strict;
use warnings;
use Carp;
use Array::IntSpan;
use List::BinarySearch qw|binsearch_pos|;
use Genome::Ranges qw|validate_elem|;

=head1 NAME

Genome::Ranges::IntSpan - *Non-overlapping* genomic intervals.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 DESCRIPTION

This is a subclass of L<Genome::Ranges>, used for *non-overlapping* genomic
intervals with region-specific information.

The same data structure is used as L<Genome::Ranges>, except that genomic 
intervals in each chromosome will be blessed into L<Array::IntSpan> object.
All the caveats of C<Array::IntSpan> object will also applies.

=head1 SYNOPSIS

    use Genome::Ranges::IntSpan;

    my $span = Genome::Ranges::IntSpan->new($ranges);
    my @overlap = $span->get_range($chrom, $start, $end);

=head1 SUBROUTINES/METHODS

=head2 $class->new RANGES

Constructor. Ranges must should be provided,
We use the validation procedure implemented in base class to validate
the input. The non-overlapness will be checked by L<Array::IntSpan>
module internally.

=cut

sub new {
	my $class = shift @_;
	my $self = $class->SUPER::new(@_);
	$self->sort;
	foreach my $chrom (keys %$self) {
		$self->{$chrom} = Array::IntSpan->new(map { [$_->[0], $_->[1], $_->[2]] } 
													@{$self->{$chrom}}); # only one value field is allowed,
    }
	return $self;
}

=head2 $self->add CHR, START, END [, VAL]

Reload add. To keep validity of the data structure, it will first test the 
interval does not overlap any existing invervals, then user binary search
to determine the index to splice in this interval. Given the performance 
penalty, it is recommended to only use this method to add a few individual
elements. Value part is an scalar. This is the convention to Array::IntSpan.

=cut

sub add {
	my ($self, $chrom, $start, $end, $val) = @_;
	$end = $start unless defined $end;
	validate_elem($chrom, $start, $end);
	my $test = [];
	if (defined $self->{$chrom}) {
		$test = $self->{$chrom}->get_range($start, $end);
	}
	unless (@$test) {
		# Determine the splice in position
		my @stpos = map { $_->[0] } @{$self->{$chrom}};
		my $index = binsearch_pos {$a <=> $b} $start, @stpos;
		splice @{$self->{$chrom}}, $index, 0, [$start, $end, $val];
		return $self;
	}
	else {
		carp "Interval $chrom:$start-$end cannot be added";
		return $self;
	}
}

#no need for sorting, just return itself.
sub sort {
	my $self = shift @_;
	return $self;
}

# all other methods clone, count, size will be reloaded.

=head2 $self->lookup CHR, POS

Wrapper of L<Array::IntSpan>'s lookup method,
Return the value at CHR:POS position.

=cut 

sub lookup {
	my ($self, $chrom, $pos) = @_;
	return undef unless defined $self->{$chrom};
  	$self->{$chrom}->lookup($pos);
}


=head2 $self->lookup CHR, POS

This is an extension of lookup method. For genomic positions that can not be found in the stored intervals.
We will instead return the value from the nearest interval. If chromosome is not stored, then the interval
will be skipped.

=cut

sub lookup_nearest {
	my ($self, $chrom, $pos) = @_;
	return undef unless defined $self->{$chrom};
	my $val = $self->{$chrom}->lookup($pos);
	if (defined $val) {
		return $val;
	}
	else {
		my @starts = map { $_->[0] } @{$self->{$chrom}};
		my @ends = map { $_->[1] } @{$self->{$chrom}};
		my $index = binsearch_pos { $a <=> $b } $pos, @starts;
		if ($index >= @starts) {
			$index --;
		}
		elsif ($index > 0) {
			my $dist_left = $pos-$ends[$index-1];
			my $dist_right = $starts[$index]-$pos;
			if ($dist_left < 0 || $dist_right < 0) {
				die "Incorrect index found at position $chrom:$pos";
			}
			if ($dist_left < $dist_right) {
				$index --;
			}
		}
		$val = $self->{$chrom}->lookup($starts[$index]);
		return $val;
	}
}


=head2 $self->get_range CHR, START, END, [...]

Wrapper of L<Array::IntSpan>'s get_range method.
When overlaps found, it will return an array for all overlapping intervals.
The interval boundaries will always be split at START and END.

Note: when value part is a reference to data structure. All split intervals will
carry the reference to the same data structure. Users can provide necessary
copy subroutines to override this behavior.

=cut

sub get_range {
	my ($self, $chrom, $start, $end, @opts) = @_;
	$end = $start unless defined $end;
	croak "Invalid query range" unless $end >= $start;
	
  	return unless defined $self->{$chrom};	
  	my $overlap = $self->{$chrom}->get_range($start, $end, @opts);
  	return @$overlap;
}

=head2 $self->set_range CHR, START, END, VAL, [...]

Wrapper of L<Array::IntSpan>'s set_range method.
By setting VAL to undef, it removes the given range.
Additional callback can be used to deal with split intervals.

The method can also be used to insert new intervals on existing
L<Array::IntSpan> object.

=cut

sub set_range {
	my ($self, $chrom, $start, $end, $val, $cb) = @_;
	croak "Invalid query range" unless $end >= $start;

	if (defined $self->{$chrom}) {
		$self->{$chrom}->set_range($start, $end, $val, $cb);
	}
	else {
		$self->{$chrom} = Array::IntSpan->new([$start, $end, $val]);
	}
	return $self;
}

=head1 AUTHOR

Xueya Zhou, C<< <xueyazhou at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-genome at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Genome>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Genome::Ranges::IntSpan


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

1; # End of Genome::Ranges::IntSpan

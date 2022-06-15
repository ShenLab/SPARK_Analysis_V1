package Genome::Ranges::IntSet;

use strict;
use warnings;
use Carp;
use Genome::Ranges;
use List::MoreUtils qw|all all any none uniq|;
use Storable qw(dclone);
use Set::IntSpan;
use Sort::Versions;
use Data::Dumper;
use Scalar::Util qw|reftype|;
use IO::Detect;
use Iterator::Simple qw|iterator|;
use Genome::Ranges qw|validate_rng validate_elem|;
use Utils::File qw|open_file|;
use Utils::Hash qw|merge_opts|;

use overload
  '+'    => 'union'     ,
  '-'    => 'diff'      ,
  '*'    => 'intersect' ,
  '^'    => 'xor'       ,
  '/'    => 'xor'       ,
  '~'    => 'complement',
  '+='   => 'U'	  ,
  '-='   => 'D'	  ,
  '*='   => 'I'	  ,
  '^='   => 'X'	  ,
  '/='   => 'X'   ,
  'bool' => sub { not shift->empty() };


=head1 NAME

Genome::Ranges::IntSet - Integer set operations on genome ranges.


=head1 DESCRIPTION

The object is implemented a hashref keyed by chromosome names. Value part is 
L<Set::IntSpan> object. The set can be constructed from a Genome::Ranges object.
The value part of each entry must be ignored.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

	use Genome::Ranges::IntSet;

	$target = Genome::Ranges::IntSet->new(\%trng);
	$cds = Genome::Ranges::IntSet->new(\%crng);

	$intersect = $target * $cds;
	$rngs = $intersect->to_ranges();


=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 $class->new [RANGES]

Create an empty new object or from L<Genome::Ranges> object.

Also allow construction from file, assuming first three columns will be
chrom, start and end. 

Add an option to read intervals from a file with BED format.

=cut

sub new {
	my ($class, $ranges, $argref) = @_;
	my $arg = merge_opts($argref, bed => undef);
	if (defined $ranges) {
		if (reftype $ranges && reftype $ranges eq 'HASH') {
			validate_rng($ranges);
			my %intset;
  			while (my ($chrom, $intvs) = each %$ranges) {
				next unless @$intvs > 0;
				my @ranges = map { [ $_->[0], $_->[1] ] } @$intvs; # remove value part
				$intset{$chrom} = Set::IntSpan->new( \@ranges );
  			}
  			return bless \%intset, $class;
		}
		elsif (reftype $ranges && reftype $ranges eq 'ARRAY') {
			my %rngs;
			foreach my $chrintv (@$ranges) {
				my ($chrom, $start, $end) = @{$chrintv}[0,1,2];
				validate_elem($chrom, $start, $end);
				push @{$rngs{$chrom}}, [$start, $end];
			}
			my %sets;
			while (my ($chrom, $intvs) = each %rngs) {
				$sets{$chrom} = Set::IntSpan->new($intvs);
  			} 
			return bless \%sets, $class;
		}
		else {
			my $fin;
			my %rngs;
			if (is_filehandle($ranges)) {
				$fin = $ranges;
			}
			elsif (-f $ranges) {
				$fin = open_file($ranges);
				#$fin = IO::File->new($ranges) or croak "Cannot open file $ranges";
			}
			else {
				croak "Cannot determine the type of $ranges";
			}
			while(<$fin>) {
				my ($chrom, $start, $end) = (split)[0,1,2];
				$start += 1 if $arg->{bed};
				validate_elem($chrom, $start, $end);
				push @{$rngs{$chrom}}, [$start, $end];
			}
			my %sets;
			while (my ($chrom, $intvs) = each %rngs) {
				$sets{$chrom} = Set::IntSpan->new($intvs);
  			} 
			return bless \%sets, $class;
		}
	}
	else {
		return bless {}, $class;
	}
}

=head2 $self->add CHROM, START, END

Add a genomic interval.

=cut

sub add {
	my ($self, $chrom, $start, $end) = @_;
	$end = $start unless defined $end;
	validate_elem($chrom, $start, $end);
	if (defined $self->{$chrom}) {
		if ($start == $end) {
			$self->{$chrom} += "$start";
		}
		else {
			$self->{$chrom} += "$start-$end";
		}
	}
	else {
		$self->{$chrom} = Set::IntSpan->new("$start-$end");
	}
	return $self;
}

=head2 $self->to_ranges 

Convert to L<Genone::Ranbges> object, with no value part in each entry.

=cut

sub to_ranges {
	my ($self) = @_;
	my %ranges;
	while ( my ($chrom, $set) = each %$self ) {
		next unless $set;
		$ranges{$chrom} = [$set->spans];
	}
	return Genome::Ranges->new(\%ranges);
}

=head2 $self->write FILE

Write genomic intervals to file.

Add option to write intervals to file in BED format.

=cut

sub write {
	my ($self, $file, $argref) = @_;
	my $arg = merge_opts($argref, bed => undef);
	my $fout;
	if (is_filehandle $file) {
		$fout = $file;
	}
	else {
		$fout = open_file($file, {write => 1});
	}
	my $it = $self->to_ranges->iter();
	while(my $dat = $it->()) {
		$dat->[1] -= 1 if $arg->{bed};
		print $fout join("\t", @$dat), "\n"; 
	}
	return $self;
}


=head2 $self->iter

Return an iterator. The value returned by calling iterator is
[chrom, pos].
To call iter, all intervals should be closed (i.e., no +/-inf in star/end).

=cut

sub iter {
	my ($self) = @_;

	my @chroms = sort { versioncmp($a, $b) } keys %$self;
	my $chr = 0;
	my $chrom = $chroms[$chr];

  	return iterator {
  		my $pos = $self->{$chrom}->next;
  		if (!defined $pos) {
  			while ( ++$chr < @chroms ) {
				$chrom = $chroms[$chr];
				$pos = $self->{$chrom}->next;
				if (defined $pos) {
					if (wantarray) {
						return ($chrom, $pos);
					}
					else {
						return [$chrom, $pos];
					}			
				}
			}
			return;
  		}
  		else {
  			if (wantarray) {
  				return ($chrom, $pos);
  			}
  			else {
  				return [$chrom, $pos];
  			}
  		}
  	};
}


=head2 $self->count [CHROM]

Number of non-overlapping intervals in the set.

=cut

sub count {
	my ($self, $chrom) = @_;
	if (defined $chrom) {
		return defined $self->{$chrom} ? scalar($self->{$chrom}->spans) : 0;
	}
	else {
		my $total = 0;
		while ( my ($chrom, $set) = each %$self ) {
			$total += scalar($set->spans);
		}
		return $total;
	}
}


=head2 $self->size [CHROM]

Number of base pairs covered by the set. Return -1 for infinite sets.

=cut

sub size {
	my ($self, $chrom) = @_;
	if (defined $chrom) {
		return defined $self->{$chrom} ? $self->{$chrom}->size : 0;
	}
	else {
		my $total = 0;
		while ( my ($chrom, $set) = each %$self ) {
			if ( $set->size() == -1 ) {
				return -1;
			} 
			else {
				$total += $set->size();
			}
		}
		return $total;
	}
}


=head2 $self->clone

  Keep a deep copy of the Genome::Range::IntSet object.

=cut

sub clone { return dclone($_[0]); }



=head2 Set operations

  Wrapper of the Set::IntSpan's subroutines:
    intersect, union, xor, diff, complement
    I        , U    , X  , D   , C

=cut

sub I {
	my ($self, $other) = @_;
	my @remove;
	foreach my $chrom (keys %$self) {
		if (!defined $other->{$chrom}) {
			push @remove, $chrom;
		}
		else {
			$self->{$chrom}->I( $other->{$chrom} );
			push @remove, $chrom unless $self->{$chrom};
		}
	}
	delete $self->{$_} foreach (@remove);
	return $self;
}

sub intersect {
	my ($self, $other) = @_;
	my $self_copy = $self->clone();
	return $self_copy->I($other);
}

sub U {
	my ($self, $other) = @_;
	foreach my $chrom (keys %$self) {
		next unless defined $other->{$chrom};
		$self->{$chrom}->U( $other->{$chrom} );
	}
	foreach my $chrom (keys %$other) {
		if ( !defined $self->{$chrom} ) {
			$self->{$chrom} = dclone($other->{$chrom});
		}
	}
	return $self;
}

sub union {
	my ($self, $other) = @_;
	my $self_copy = $self->clone();
	return $self_copy->U($other);
}

sub X {
	my ($self, $other) = @_;
	foreach my $chrom (keys %$self) {
		next unless defined $other->{$chrom};
		$self->{$chrom}->X( $other->{$chrom} );
	}
	foreach my $chrom (keys %$other) {
		if ( !defined $self->{$chrom} ) {
			$self->{$chrom} = dclone $other->{$chrom};
		}
	}
	return $self;
}

sub xor {
	my ($self, $other) = @_;
	my $self_copy = $self->clone();
	return $self_copy->X($other);
}

sub D {
	my ($self, $other) = @_;
	my @remove;
	foreach my $chrom (keys %$self) {
		next unless $other->{$chrom};
		$self->{$chrom}->D( $other->{$chrom} );
	}
	return $self;
}

sub diff {
	my ($self, $other) = @_;
	my $self_copy = $self->clone();
	return $self_copy->D($other);
}

sub C {
	my ($self) = @_;
	foreach my $chrom (keys %$self) {
		$self->{$chrom}->C();
	}
	return $self;
}

sub complement {
	my ($self) = @_;
	my $self_copy = $self->clone();
	return $self_copy->C();
}

=head2 AUTOLOADed Operators

  Comparison  : equal, equivalent, superset, subset
  Generator   : cover and holes, inset(trim) and pad
  Membership  : member, insert, remove

=cut

sub AUTOLOAD {
	my @comaprison = qw(equal equivalent superset subset);
	my @generator  = qw(trim pad inset cover holes);
	my @membership = qw(member insert remove);
	my @cardi_all  = qw(empty finite);
	our $AUTOLOAD;
	no strict 'refs';
	(my $method = $AUTOLOAD) =~ s/.*:://s; # remove package name
  	return if $method eq 'DESTROY'; # ignore DESTORY method

  	if ( $method eq 'equal' ) {
  		*{$AUTOLOAD} = sub {
  			my $self  = shift @_;
  			my $other = shift @_;
  			my @self = sort keys %$self;
  			my @other = sort keys %$other;
  			if (@self == @other &&
  				( all { $self[$_] eq $other[$_] } 0..$#self ) &&
  				( all { $self->{$_}->equal($other->{$_}) } @self )
  				) {
  				return 1;
  			}
  			else {
  				return 0;
  			}
  		};
  	}
  	elsif ($method eq 'superset') {
  		*{$AUTOLOAD} = sub {
  			my $self  = shift @_;
  			my $other = shift @_;
  			my @other = sort keys %$other;
  			if (( all { defined $self->{$_} } @other ) &&
  				( all { $self->{$_}->superset($other->{$_}) } @other )
  				) {
  				return 1;
  			}
  			else {
  				return 0;
  			}
  		};
  	}
  	elsif ($method eq 'subset') {
  		*{$AUTOLOAD} = sub {
  			my $self  = shift @_;
  			my $other = shift @_;
  			return (all { defined $other->{$_} && $self->{$_}->$method($other->{$_}) } keys %$self);
  		};
  	}
  	elsif ( grep {$method eq $_} @generator ) {
  		*{$AUTOLOAD} = sub {
  			my $self  = shift @_;
  			my %newset;
  			foreach my $chrom (keys %$self) {
  				$newset{$chrom} = $self->{$chrom}->$method(@_);
  			}
  			return bless \%newset, ref $self;
  		};
  	}
  	elsif ( grep {$method eq $_ } @membership ) {
  		*{$AUTOLOAD} = sub {
  			my $self  = shift @_;
  			my $chrom = shift @_;
  			return undef unless defined $self->{$chrom};
  			return $self->{$chrom}->$method(@_);
  		};
  	}
  	elsif ( grep {$method eq $_ } @cardi_all ) {
  		*{$AUTOLOAD} = sub {
  			my $self  = shift @_;
  			return (all { $self->{$_}->$method() } keys %$self);
  		};
  	}
  	else {
  		croak "I does not understand $method";
  	}
  	goto &{$AUTOLOAD};
}


=head1 AUTHOR

Xueya Zhou, C<< <xueyazhou at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-genome at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Genome>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Genome::Ranges::IntSet


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

1; # End of Genome::Ranges::IntSet

package Utils::List;

use strict;
use warnings;
use Carp;
use Data::Dumper;
use List::Util qw|first product|;
use List::MoreUtils qw|all|;
use Set::IntSpan;


use base qw|Exporter|;
our @EXPORT_OK = qw|which_min which_max run_list all_pairs all_combs replace_first replace_all
					insert_before insert_after parse_fields split_list split_list_fixed|;
our %EXPORT_TAGS = ( 'all' => \@EXPORT_OK );

=head1 NAME

Utils::List - Simple array utilities

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Implement some simple useful functions for list which are missed in other packages.

See also C<List::Util> and C<List::MoreUtils>

	$ii = which_min(@data);
	$maxval = $data[$ii];

=head1 EXPORT

=over 5

=item * which_max

=item * which_min

=item * run_list

=back

=head1 SUBROUTINES/METHODS

=head2 which_max/min DATA

Return the index of max/min element.

In scalar context, return the B<first> index; in array context, return all index.

=cut

sub which_min {
  my @data =  @_;
  return unless @data > 0;
  my $min_ii = 0;
  for(my $ii = 1; $ii < @data; $ii ++) {
	if ($data[$ii] < $data[$min_ii]) {
	  $min_ii = $ii;
	}
  }
  if (wantarray) {
	my @min_ii = ();
	for(my $ii = 0; $ii < @data; $ii ++) {
	  if ($min_ii == $ii) {
		push @min_ii, $ii;
	  }
	  else {
		if ($data[$ii] == $data[$min_ii]) {
		  push @min_ii, $ii;
		}
	  }
	}
	return @min_ii;
  }
  else {
	return $min_ii;
  }
}

sub which_max {
  my @data = @_;
  return unless @data > 0;
  my $max_ii = 0;
  for(my $ii = 1; $ii < @data; $ii ++) {
	if ($data[$ii] > $data[$max_ii]) {
	  $max_ii = $ii;
	}
  }
  if (wantarray) {
	my @max_ii = ();
	for(my $ii = 0; $ii < @data; $ii ++) {
	  if ($max_ii == $ii) {
		push @max_ii, $ii;
	  }
	  else {
		if ($data[$ii] == $data[$max_ii]) {
		  push @max_ii, $ii;
		}
	  }
	}
	return @max_ii;
  }
  else {
	return $max_ii;
  }
}


=head2 run_list INTEGERS...

Return run list for an array of integers.

=cut

sub run_list {
	my @nints = @_;
	my $set;
	if (@nints == 1 && ref $nints[0] eq 'ARRAY') {
		$set = Set::IntSpan->new(@{$nints[0]});
	}
	else {
		$set = Set::IntSpan->new(@nints);
	}
	my $runs = $set->run_list;
	if (wantarray) {
		 my @ranges = split(',', $runs);
		 return @ranges;
	}
	else {
		return $runs;
	}
}

=head2 all_pairs

Return all pairswise combination of different elements in a list.

=cut

sub all_pairs {
	my @list = @_;
	my @pairs;
	for(my $ii = 0; $ii < @list; $ii ++) {
		for(my $jj = $ii + 1; $jj < @list; $jj ++) {
			push @pairs, [$list[$ii], $list[$jj]];
		}
	}
	if (wantarray) {
		return @pairs;
	}
	else {
		return \@pairs;
	}
}

=head2 all_combs

Given multiple arrays of different sizes, return the full list of all combinations
that have one element from each array.

Example: @a=qw(A B C); @b=qw(1 2 3)
then all_combs(\@a, \@b) returnes [A,1], [B,1], [C,1], [A,2], [B,2], [C,2], [A,3], [B,3], [C,3] 

=cut

sub all_combs {
	my (@arrays, @labels);
	if(all { ref $_ eq 'ARRAY'  } @_) {
		@arrays = @_;
	}
	else {
		my %data = @_;
		if ( (all { ref $_ eq 'ARRAY' } values %data) && 
			 (all { ref $_ eq '' } keys %data) ) {
			@labels = sort keys %data;
			@arrays = @data{@labels};
		}
		else {
			die "Incorrect arguments to all_combs: should be an array of arrayref or hash of arraref";
		}
	}

	my @sizes = map { scalar(@$_) } @arrays;
	my @combs;
	for(my $ii = 0; $ii < product(@sizes); $ii ++) {
		my @mind = _index_sing2multi($ii, @sizes);
		if (@labels) {
			push @combs, { map { ($labels[$_] => $arrays[$_][$mind[$_]]) } 0..$#arrays };
		}
		else {
			push @combs, [map { $arrays[$_][$mind[$_]] } 0..$#arrays];
		}
	}
	if (wantarray) {
		return @combs;
	}
	else {
		return \@combs;
	}
}

sub _index_sing2multi {
	my $single = shift @_;
	my @sizes = @_;
	my $bound = product(@sizes);
	if ($single >= $bound) {
		die "Single index is out of bound: $bound";
	}
	my @multi;
	my ($divide, $remain);
	while(@sizes > 0) {
		my $size = pop @sizes;
		if (@sizes > 0) {
			my $prod = product(@sizes);
			$divide = int($single/$prod);
			$remain = $single % $prod;
			unshift @multi, $divide;
		}
		else {
			$remain = $single;
			if ($remain >= $size) {
				die "Remain is larger than size: $remain >= $size";
			}
			unshift @multi, $remain;
		}
		$single = $remain;
	}
	return @multi;
}


=head2 replace_first/replace_all

Replace the element at given mark in an array. When multiple marks are found
the replac_first only replace the first occurence, whereas replace_all will
replace all marks found. Return the number of replaced marks.

=cut

sub replace_first {
	my ($listref, $mark, $content) = @_;
	my $ii = first { $listref->[$_] eq $mark } 0..$#$listref;
	unless (defined $ii) {
		carp "Cannot find $mark in list";
		return 0;
	}
	$listref->[$ii] = $content;
	return 1;
}

sub replace_all {
	my ($listref, $mark, $content) = @_;
	my @iis = grep { $listref->[$_] eq $mark } 0..$#$listref;
	unless(@iis) {
		carp "Cannot find $mark in list";
		return 0;
	}
	foreach my $ii (@iis) {
		$listref->[$ii] = $content;
	}
	return scalar(@iis);
}


=head2 insert_before/after LISTREF, POS, CONTENTS...

Insert elements before or after certain mark in an array. When multiple marks are found,
the first one will be used. Return the index of inserted element in the new array.

=cut

sub insert_before {
	my ($listref, $mark, @contents) = @_;
	my $ii = first { $listref->[$_] eq $mark } 0..$#$listref;
	unless (defined $ii) {
		carp "Cannot find $mark in list";
		return;
	}
	if ($ii == 0) {
		unshift @$listref, @contents;
	}
	else {
		splice @$listref, $ii, 0, @contents; 
	}
	return $ii;
}

sub insert_after {
	my ($listref, $mark, @contents) = @_;
	my $ii = first { $listref->[$_] eq $mark } 0..$#$listref;
	unless (defined $ii) {
		carp "Cannot find $mark in list";
		return;
	}
	if ($ii == $#$listref) {
		push @$listref, @contents;
	}
	else {
		splice @$listref, $ii+1, 0, @contents; 
	}
	return $ii+1;
}

=head2 parse_fields

Parse a comma separated list of field names, e.g. Chr,Pos,St,Ed.
If comma is part of field names, escape is allowed, e.g. Chr,Pos\,hg19,Start,End.

=cut

sub parse_fields {
	my ($fstr) = @_;
	# Parse fields
	my @fields = split(/(?<!\\),/, $fstr);
	# The look for escaped commas
	foreach my $field (@fields) {
		if ($field =~ /(?<!\\)\\,/) {
			$field =~ s/(?<!\\)\\,/,/g;
		} 
		if ($field =~ /(?<=\\)\\/) {
			$field =~ s/\\\\/\\/g;
		}
	}
	return @fields;
}

=head2 split_list ARRAREF, SIZE

Split a list into roughly similar sized segments.

The SIZE parameter is only a suggest, the actual size of split segments may 
be different from this value in order to keep the size of all segements as
similar as possible.

=cut

sub split_list {
	my ($list, $size) = @_;
	my @segs;
	my $nsplit = int(scalar(@$list)/$size+0.5);
	if ($nsplit > 1) {
		my $actsize = int(scalar(@$list)/$nsplit+0.5);
		for(my $ii = 0; $ii < $nsplit; $ii ++) {
			my $st = $ii * $actsize;
			my $ed;
			if ($ii == $nsplit - 1) {
				$ed = @$list - 1;
			}
			else {
				$ed = ($ii+1) * $actsize - 1;
			}
			push @segs, [@{$list}[$st..$ed]];
		}
	}
	else {
		print STDERR "Nsplit=$nsplit: no need to split?\n";
		@segs = ($list);
	}
	if (wantarray) {
		return @segs;
	} 
	else {
		return \@segs;
	}
}

# Second version is to split list into fixed number of segments
sub split_list_fixed {
	my ($list, $nsplit) = @_;
	my @segs;
	if ($nsplit > 1) {
		my $size = int(scalar(@$list)/$nsplit+0.5);
		for(my $ii = 0; $ii < $nsplit; $ii ++) {
			my $st = $ii * $size;
			my $ed;
			if ($ii == $nsplit - 1) {
				$ed = @$list - 1;
			}
			else {
				$ed = ($ii+1) * $size - 1;
			}
			push @segs, [@{$list}[$st..$ed]];
		}
	}
	else {
		print STDERR "Nsplit=$nsplit: no need to split?\n";
		@segs = ($list);
	}
	if (wantarray) {
		return @segs;
	} 
	else {
		return \@segs;
	}
}



=head1 AUTHOR

Xueya Zhou, C<< <xueyazhou at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-utils at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Utils>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Utils::List


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Utils>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Utils>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Utils>

=item * Search CPAN

L<http://search.cpan.org/dist/Utils/>

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

1; # End of Utils::List

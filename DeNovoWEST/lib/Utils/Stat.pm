package Utils::Stat;

use strict;
use warnings;
use Carp;
use List::Util qw|min max sum|;
use POSIX qw|floor ceil|;
use Iterator::Simple qw|iterator|;
use Utils::Hash qw|merge_opts|;

use base qw|Exporter|;
our @EXPORT_OK = qw|mean weightedmean median mad var std sum_sq corr chisq_2x2
					index_hist quantile percent
					sample rand_perm permute|;
our %EXPORT_TAGS = ( 'all' => \@EXPORT_OK );


=head1 NAME

Utils::Stat - Simple statistical utilities.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Functions for doing basic statistical analysis.
See also: C<Statistics::Descriptive>, C<Statistics::Distributions>

=head1 EXPORT

=over 5

=item * mean

=item * median

=item * mad

=item * var

=item * std

=item * sum_sq

=item * corr

=item * chisq_2x2

=item * quantile

=item * percent

=item * index_hist

=item * sample

=item * permute

=item * rand_permute

=back

=head1 SUBROUTINES/METHODS

By default the input for this module is array or array ref with no missing data.

Alternatively, when data are integers, it can also be represented as a hash where
the number of appearance (counts) of each integer is stored as key-value pairs.

=head2 mean

Calculate mean for an array of numbers.
The input can be array, array ref, or hist (hash ref).

=cut

sub mean {
	if (@_) {
		if (@_ == 1) {
			if (ref $_[0] eq 'ARRAY') {
				return sum(@{$_[0]})/scalar(@{$_[0]});
			}
			elsif (ref $_[0] eq 'HASH') {
				my ($num, $den);
				while(my ($val, $count) = each %{$_[0]}) {
					$num += $val*$count;
					$den += $count;
				}
				unless ($den) {
					return undef;
				}
				else {
					return $num/$den;
				}
			}
			else {
				unless(ref $_[0]) {
					return $_[0];
				}
				else {
					croak "Unknown object type : ". ref $_[0];
				}
			}
		}
		else {
			return sum(@_)/scalar(@_);
		}
	}
	else {
		return undef;
	}
}

=head2 weightedmean

Calculate weighted mean for an array of numbers
The input should be two arrayrefs of the same length, first for values, second for weights.

=cut

sub weightedmean {
	my ($values, $weights) = @_;
	unless(ref $values eq 'ARRAY' && ref $weights eq 'ARRAY') {
		croak "weigthedmean: must provide two array refs for values and weights";
	}
	unless(scalar(@$values) == scalar(@$weights)) {
		croak "values and weights must be of the same length!";
	}
	my ($num, $den);
	for(my $ii = 0; $ii < @$values; $ii ++) {
		$num += $values->[$ii]*$weights->[$ii];
		$den += $weights->[$ii];
	}
	unless ($den) {
		return undef;
	}
	else {
		return $num/$den;
	}
}

=head2 median

Calculate median (50% percentile) for an array of numbers.

=cut

sub median {
	if (@_) {
		my @x;
		if (@_ == 1) {
			if (ref $_[0] eq 'ARRAY') {
				@x = sort { $a<=>$b } @{$_[0]};
				my $mind = int(scalar(@x)/2);
				if (@x % 2 == 1) {
					return $x[$mind];
				}
				else {
					return 0.5*($x[$mind-1]+$x[$mind]);
				}
			}
			elsif (ref $_[0] eq 'HASH') {
				return quantile($_[0], 0.5);
			}
			else {
				croak "Unknown object type : ". ref $_[0];
			}
		}
		else {
			@x = sort { $a<=>$b } @_;
			my $mind = int(scalar(@x)/2);
			if (@x % 2 == 1) {
				return $x[$mind];
			}
			else {
				return 0.5*($x[$mind-1]+$x[$mind]);
			}
		}
	}
	else {
		return undef;
	}
}

=head2 mad

Median absolute deviation

=cut
	
sub mad {
	my $med = median(@_);
	if (defined $med) {
		my @dat;
		if (@_ == 1) {
			if (ref $_[0] eq 'ARRAY') {
				@dat = sort { $a<=>$b } map { abs($_-$med) } @{$_[0]};
			}
			elsif (ref $_[0] eq 'HASH') {
				@dat = sort { $a<=>$b } map { abs($_-$med) } values %{$_[0]};
			}
			else {
				croak "Unknown object type : " . ref $_[0];
			}
		}
		else {
			@dat = sort { $a<=>$b } map { abs($_-$med) } @_;
		}
		return median(\@dat);
	}
	else {
		return undef;
	}
}


=head2 var

Calculate sample variance for an array of numbers.

=cut

sub var {
	if (@_) {
		my @x;
		if (@_ == 1) {
			if (ref $_[0] eq 'ARRAY') {
				return sum_sq($_[0], $_[0])/(scalar(@{$_[0]})-1);
			}
			elsif (ref $_[0] eq 'HASH') {
				my $mean_h = mean($_[0]);
				return undef unless defined $mean_h;
				my ($num, $den);
				while(my ($val, $count) = each %{$_->[0]}) {
					$num += (($val-$mean_h)**2)*$count;
					$den += $count;
				}
				if ($den > 1) {
					return $num/($den-1);
				}
				else {
					return undef;
				}
			}
			else {
				croak "Unknown object type : ". ref $_[0];
			}
		}
		else {
			return sum_sq(\@_, \@_)/(scalar(@_)-1);
		}
	}
	else {
		return undef;
	}
}

=head2 std

Calculate sample standard deviation.

=cut

sub std {
	if (@_) {
		my $var = var(@_);
  		if (defined $var && $var >= 0) {
			return sqrt($var);
  		}
  		else {
			return undef;
  		}
  	}
}


=head2 sum_sq X, Y

Sum of product deviation to mean. When X==Y, it is sum squared deviation.
X and Y should be array ref.


=cut

sub sum_sq {
	my ($x, $y) = @_;
	croak "Input must be array ref" unless ref $x eq 'ARRAY' && ref $y eq 'ARRAY';
	my $mean_x = mean(@$x);
	my $mean_y = mean(@$y);
	my $sum = 0;
	for (my $ii = 0; $ii < @$x; $ii ++) {
		$sum += ($x->[$ii]-$mean_x)*($y->[$ii]-$mean_y);
	}
	return $sum;
}

=head2 corr_r X, Y

Pearson correlation coefficient (r) for X and Y. X, Y should be array ref.

NOTE: r2 is also known as coefficient of determination.

=cut

sub corr {
	my ($x, $y) = @_;
	croak "Input must be array ref" unless ref $x eq 'ARRAY' && ref $y eq 'ARRAY';
	my $ssxx = sum_sq($x, $x);
	my $ssyy = sum_sq($y, $y);
	my $ssxy = sum_sq($x, $y);
	if ($ssxx > 0 && $ssyy> 0) {
		return $ssxy/sqrt($ssxx*$ssyy);
	}
	else {
		return undef;
	}
}

=head2 chisq_2x2 A, B, C, D

Chi-squared statistics from 2x2 table. It can be used together with chi-squared
distribution to find p-values.

=cut

sub chisq_2by2 {
  my ($a,$b,$c,$d) = @_;
  return (($a*$d-$b*$c) ** 2 * ($a+$b+$c+$d)) / (($a+$b)*($c+$d)*($b+$d)*($a+$c))
}

=head2 quantile DATA, P [,OPTIONS]

Estimate quantile for a given sample. The input data array must be sorted.

=head3 Options

The function assume that input data is not sorted and will do the sorting 
automatically. However, if the data is already sorted and performance is an issue,
you can disable sorting by C<sort=0>. 

DATA can be either array, arrayref, or hashref (as histogram). Sorting is
irrelevant for histogram.

P can be scalar or an arraref to provide multiple percentiles. In case when P
is arrayref, the output will be array or arrayref depending on context.

You can also change the quantile weighting medthod by the C<method> option 
(see tech details below).

=head3 Details

All sample quantiles are defined as weighted averages of
consecutive order statistics. Sample quantiles of type i are
defined by:
           Q(p) = (1 - gamma) x[j] + gamma x[j+1],
  where (j-m)/n <= p < (j-m+1)/n, x[j] is the jth order statistic,
  n is the sample size, the value of gamma is a function of
  j = floor(np + m) and g = np + m - j, and m is a constant
We used m = 1-p is the default used by R and S (method=7).

=cut 

sub _quantile_helper {
	my ($n, $p, $type) = @_;
	$type = 7 unless defined $type;
	my ($m, $j, $g, $gamma);
	# Discontnous sample quantile type 1~3
	if ($type >= 1 && $type <= 3) {
		if ($type == 1 || $type == 2) {
			$m = 0;
		}
		elsif ($type == 3) {
			$m = -0.5;
		}
		else {
			croak "Unsupported type $type";
		}
		$j = floor($n*$p+$m);
		$g = $n*$p+$m-$j;
		if ($type == 1) {
			$gamma = $g == 0 ? 0 : 1;
		}
		elsif ($type == 2) {
			$gamma = $g == 0 ? 0.5 : 1;
		}
		else {
			$gamma = ($g == 0 && $j % 2 == 0) ? 0 : 1;
		}
	}
	elsif ($type >= 4 && $type <= 9) {
		if ($type == 4) {
			$m = 0;
		}
		elsif ($type == 5) {
			$m = 0.5;
		}
		elsif ($type == 6) {
			$m = $p;
		}
		elsif ($type == 7) {
			$m = 1-$p;
		}
		elsif ($type == 8) {
			$m = ($p+1)/3.0;
		}
		elsif ($type == 9) {
			$m = $p/4.0+0.375;
		}
		else {
			croak "Unsupported type $type";
		}
		$j = floor($n*$p+$m);
		$gamma = $n*$p+$m-$j;
	}
	else {
		croak "Unsupported type $type";
	}
	return ($j, $gamma);
}

sub quantile {
	my ($data, $p, $argref) = @_;
	my $arg = merge_opts($argref, sort => 1, method => 7);
	if (ref $data eq 'ARRAY') {
		my $x;
		if ($arg->{sort}) {
			my @x = sort { $a<=>$b } @$data;
			$x = \@x;
		}
		else {
			$x = $data;
		}
		my $n = scalar(@$x);
		if (ref $p eq 'ARRAY') {
			my @res = map { 
				my ($j, $gamma) = _quantile_helper($n, $_, $arg->{method});
				((1-$gamma)*$x->[$j-1]+$gamma*$x->[$j])
			} @$p;
			if (wantarray) {
				return @res;
			}
			else {
				return \@res;
			}
		}
		else {
			my ($j, $gamma) = _quantile_helper($n, $p, $arg->{method});
			return ((1-$gamma)*$x->[$j-1]+$gamma*$x->[$j]);
		}	
	}
	elsif (ref $data eq 'HASH') {
		my $n = sum(values %$data);
		if (ref $p eq 'ARRAY') {
			my @res = map {
					my ($j, $gamma) = _quantile_helper($n, $_, $arg->{method});
					(1-$gamma)*index_hist($data,$j)+$gamma*index_hist($data,$j+1)
				} @$p;
		}
		else {
 	 		my ($j, $gamma) = _quantile_helper($n, $p, $arg->{method});
  			return((1-$gamma)*index_hist($data,$j)+$gamma*index_hist($data,$j+1));
		}
	}
	else {
		croak "Unknown object type : ".ref $data;
	}
}

=head2 percent DATA

Input data can be array ref or hash ref (histogram).

It calculates rank percentiles for each value and returns arrayref or hashref
for each values in the original data. 

By default, larger percentile corrrsponds to larger value.
If you want larger value to ranked at top (i.e, smaller rank), use C<rev=1> option.

Similar to C<quantile>, C<percent> also assume data are not sorted (from smallest to
largest). You can disable sorting by sort=0. 

By default is rank percentile RANK/SIZE. You can add psuedo counts: (RANK+P)/(SIZE+P).
You can give P by C<pseudo> option.


=cut 

# In R:
# perc.rank <- function(x, xo) length(x[x <= xo])/length(x)*100

sub percent {
	my ($data, $argref) = @_;
	my $arg = merge_opts($argref, rev => 0, sort => 1, pseudo => 0);
	if (ref $data eq 'ARRAY') {
		my $x;
		if ($arg->{sort}) {
			my @x = sort { $a<=>$b } @$data;
			$x = \@x;
		}
		else {
			$x = $data;
		}
		my %rnk;
		# Note: ties are handled automatically
		if ($arg->{rev}) {
			for(reverse(0..$#$x)) {
				$rnk{$x->[$_]} = (@$x-$_+$arg->{pseudo})/(@$x+$arg->{pseudo}); 
			}
		}
		else {
			for(0..$#$x) {
				$rnk{$x->[$_]} = ($_+1+$arg->{pseudo})/(@$x+$arg->{pseudo}); 
			}
		}
		if (wantarray) {
			return map { $rnk{$_} } @$data;
		}
		else {
			return [map { $rnk{$_} } @$data];
		}
	}
	elsif (ref $data eq 'HASH') {
		my @num = keys %$data;
		my $sum = sum(values %$data);
		if (!defined $sum || $sum <= 0) {
			return undef;
		}
		my %acc;
		if ($arg->{rev}) {
			foreach my $val (@num) {
				$acc{$val} = sum(map { $data->{$_} } grep { $_ >= $val } @num);
			}
		}
		else {
			foreach my $val (@num) {
				$acc{$val} = sum(map { $data->{$_} } grep { $_ <= $val } @num);
			}
		}
		if (wantarray) {
			return map { $_ => ($acc{$_}+$arg->{pseudo})/($sum+$arg->{pseudo}) } @num;
		}
		else {
			return { map { $_ => ($acc{$_}+$arg->{pseudo})/($sum+$arg->{pseudo}) } @num };
		}
	}
	else {
		croak "Unknown object type : ".ref $data;
	}
}

=head2 index_hist HIST, II

Helper function: the ii-th sorted values from the hist hash ref.

=cut

sub index_hist {
	my ($hist, $ii) = @_;
	my $acc = 0;
	foreach my $val (sort { $a <=> $b } keys %$hist) {
		$acc += $hist->{$val};
		if ($acc >= $ii) {
			return $val;
		}
	}
}

=head2 sample LIST [, OPTIONS]

Random sample elements with or without replacement from LIST. 

* C<nsub>: use this option to sample a subset from original list. If undefined,
then all items will be sampled.

* C<replace>: sample with (1) or without (0 or undef) replacement, default 1.

When sampling without replacement, nsub must be provided and less than total.
If in this case nsub == total, you should use rand_perm for random permutation.

NOTE: It's a good practie to set seeds by srand before running this function.

=cut

# random sample with replacement
sub sample {
	my ($list, $argref) = @_;
	croak "Input list should be array ref" unless ref $list eq 'ARRAY';
	my $arg = merge_opts($argref, nsub => undef, replace => 1);
	if ($arg->{replace}) {
		if (wantarray) {
			return @{$list}[_sample_wr_ind(scalar @$list, $arg->{nsub})];
		}
		else {
			return [ @{$list}[_sample_wr_ind(scalar @$list, $arg->{nsub})] ];
		}
	}
	else {
		if (wantarray) {
			return @{$list}[_sample_ind(scalar @$list, $arg->{nsub})];
		}
		else {
			return [ @{$list}[_sample_ind(scalar @$list, $arg->{nsub})] ];	
		}
	}
}

# helper function
sub _sample_wr_ind {
	my ($n, $sub) = @_;
	$sub = $n unless defined $sub;
	croak "Subcount should be no more than total" unless $n >= $sub && $sub > 0;
	my @result;
	for my $i (0..$sub-1) {
		$result[$i] = int(rand($n));
	}
	return @result;
}

sub _sample_ind {
	my ($n, $k) = @_;
	$k = 1 unless defined $k;
	croak "Subcount should be less than total" unless $n > $k && $k > 0;
	my @result;
	if ($k*6 > $n) {
		my @pool = (0..$n-1);
		for my $i (0..$k-1) {
			my $j = int(rand($n-$i));
			$result[$i] = $pool[$j];
			$pool[$j] = $pool[$n-$i-1];
		}
	}
	else {
		my %selected;
		for my $i (0..$k-1) {  # invariant:  non-selected at [0,n-i)
			my $j = int(rand($n));
			while (defined $selected{$j}) {
				$j = int(rand($n));
			}
			$result[$i]=$selected{$j}=$j;
		}
	}
	return @result;
}

=head2 rand_perm

Return a random permutation of list.

Set seed before call this function.

=cut

sub rand_perm {
	if (@_ == 1) {
		if (ref $_[0] eq 'ARRAY') {
			return sort { int(rand 2)*2-1 } @{$_[0]};	
		}
		else {
			croak "Unknown object type : ". ref $_[0];
		}
	}
	else {
		if (wantarray) {
			return sort { int(rand 2)*2-1 } @_;	
		}
		else {
			return [ sort { int(rand 2)*2-1 } @_ ];	
		}
	}
}

=head2 permute

The permutation of a given list. Return a iterator that generate
one permutation at one time.

Usage:

	my $iter = permute(@b)
	while(@a = iter->()) {
	... do things with @a ...
	}

=cut

sub _n_to_pat {
	my @odometer;
	my ($n, $length) = @_;
	for my $i (1 .. $length) {
		unshift @odometer, $n % $i;
		$n = int($n/$i);
	}
	return $n ? () : @odometer;
}

sub _pattern_to_permutation {
	my $pattern = shift;
	my @items = @{shift()};
	my @r;
	for (@$pattern) {
		push @r, splice(@items, $_, 1);
	}
	@r;
}

sub permute {
	my @items = @_;
	my $n = 0;
	return iterator {
		my @pattern = _n_to_pat($n, scalar(@items));
		my @result = _pattern_to_permutation(\@pattern, \@items);
		$n++;
		if (grep { $_ } @result) {
			if (wantarray) {
				return @result;
			}
			else {
				return \@result;
			}
		}
		else {
			return undef;
		}
		#return @result;
	};
}



=head1 AUTHOR

Xueya Zhou, C<< <xueyazhou at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-utils at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Utils>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Utils::Stat


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

1; # End of Utils::Stat

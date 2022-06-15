package Utils::Stat::Running;

use strict;
use warnings;
use Carp;

# Ref: https://www.johndcook.com/blog/skewness_kurtosis

sub new {
	my ($class) = @_;
	return bless { N => 0, M1 => 0.0, M2 => 0.0, M3 => 0.0, M4 => 0.0 }, $class;
}

sub clear {
	my $self = shift;
	foreach my $field (qw|N M1 M2 M3 M4|) {
		$self->{$field} = 0;
	}
	return $self;
}

=head2 push NUM...

The algorithm in C:
    n++;
    delta = x - M1;
    delta_n = delta / n;
    delta_n2 = delta_n * delta_n;
    term1 = delta * delta_n * n1;
    M1 += delta_n;
    M4 += term1 * delta_n2 * (n*n - 3*n + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;
    M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2;
    M2 += term1;

=cut

sub push {
	my $self = shift;
	foreach my $x (@_) {
		my ($delta, $delta_n, $delta_n2, $term1);
 		my $n1 = $self->{N};
 		$self->{N} ++;
    	$delta = $x - $self->{M1};
    	$delta_n = $delta / $self->{N};
    	$delta_n2 = $delta_n * $delta_n;
    	$term1 = $delta * $delta_n * $n1;
    	$self->{M1} += $delta_n;
    	$self->{M4} += $term1 * $delta_n2 * ($self->{N}**2 - 3*$self->{N} + 3) + 6 * $delta_n2 * $self->{M2} - 4 * $delta_n * $self->{M3};
    	$self->{M3} += $term1 * $delta_n * ($self->{N} - 2) - 3 * $delta_n * $self->{M2};
    	$self->{M2} += $term1;
	}
	return $self;
}

sub npts {
	my $self = shift;
	return $self->{N};
}

sub mean {
	my $self = shift;
	return $self->{M1};
}

sub var {
	my $self = shift;
	unless($self->{N} > 1) {
		die "Cannot calculate variance for less than two data points";
	}
	return $self->{M2}/($self->{N}-1);
}

sub std {
	my $self = shift;
	return sqrt($self->var());
}

sub skew {
	my $self = shift;
	return sqrt($self->{N}*$self->{M3}/$self->{M2}**1.5);
}

sub kurtosis {
	my $self = shift;
	return $self->{N}*$self->{M4}/$self->{M2}**2-3.0;
}


1;


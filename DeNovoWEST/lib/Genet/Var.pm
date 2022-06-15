package Genet::Var;

use strict;
use warnings;
use Carp;
use Data::Dumper;
use List::MoreUtils qw|all uniq|;
use Iterator::Simple qw|iterator|;

use base qw|Exporter|;

our @EXPORT_OK = qw|is_normalized is_cpg normalize left_trim right_trim_or_left_extend|;

=head2 is_normalized VAR

Check if a variant is normalized.
Note: this function does not check if the variant is left aligned.

=cut

sub _unpack_alt {
	my ($var) = @_;
	unless (all { defined $var->{$_} } qw|CHROM POS REF ALT|) {
		print Dumper $var;
		croak "A variant must contain CHROM/POS/REF/ALT fields";
	}
	my @alleles;
	unless (ref $var->{ALT}) {
		@alleles = ($var->{REF}, $var->{ALT});
	}
	else {
		unless (ref $var->{ALT} eq 'ARRAY') {
			croak "Multiple alt alleles should be stored in an array";
		}
		@alleles = ($var->{REF}, @{$var->{ALT}});
	}
	return @alleles;
}

sub is_normalized {
	my ($var) = @_;

	my @alleles = _unpack_alt($var);
	
	my ($first_base, $last_base, $exist_len_one_allele);
	my ($first_base_same, $last_base_same, $same) = (1,1,1);

	for(my $ii = 0; $ii < @alleles; $ii ++) {
		next if $alleles[$ii] eq '*';
		my $len = length($alleles[$ii]);
		$exist_len_one_allele = 1 if $len == 1;
		if ($ii) {
			if ($first_base ne uc(substr($alleles[$ii], 0, 1))) {
				$first_base_same = 0;
			}
			if ($last_base ne uc(substr($alleles[$ii], -1))) {
				$last_base_same = 0;
			}
			$same = $same & uc($alleles[$ii]) eq uc($alleles[0]);
		}
		else {
			$first_base = uc(substr($alleles[0], 0, 1));
			$last_base = uc(substr($alleles[0], -1));
		}
	}
	if ($same) {
		return 1;
	}
	if ($last_base_same || (!$exist_len_one_allele && $first_base_same)) {
		return 0;
	}
	return 1;
}


=head2 is_cpg VAR, SEQ

Check is a SNV is in the context of CpG hypermutable site.

CpG => TpG  

=cut

sub is_cpg {
	my ($var, $sq) = @_;
	if ($var->{REF} eq 'C' && $var->{ALT} eq 'T') {
		my $flankbase = $sq->get_slice($var->{CHROM}, $var->{POS}+1, $var->{POS}+1);
		if ($flankbase eq 'G') {
			return 1;
		}
		else {
			return 0;
		}
	}
	elsif ($var->{REF} eq 'G' && $var->{ALT} eq 'A') {
		my $flankbase = $sq->get_slice($var->{CHROM}, $var->{POS}-1, $var->{POS}-1);
		if ($flankbase eq 'C') {
			return 1;
		}
		else {
			return 0;
		}
	}
	else {
		return 0;
	}
}

=head2 norm_var VAR

Normalize variant by applying the following two functions

=cut

sub normalize {
	my ($var, $sq, $debug) = @_;
	return left_trim(right_trim_or_left_extend($var, $sq, $debug), $debug);
}


=head2 left_trim VAR

Left trim a variant.

=cut 

sub left_trim {
	my ($var, $debug) = @_;

	my @alleles = _unpack_alt($var);

	my $to_left_trim = 1;
	my $left_trimmed = 0;

	while($to_left_trim) {
		for(my $ii = 0; $ii < @alleles; $ii ++) {
			next if $alleles[$ii] eq '*';
			if (length($alleles[$ii]) == 1 || 
				substr($alleles[$ii], 0, 1) ne substr($alleles[0], 0, 1)) {
				$to_left_trim = 0;
				last;
			}
		}
		if ($to_left_trim) {
			for (my $ii = 0; $ii < @alleles; $ii ++) {
				substr($alleles[$ii], 0, 1, "");
			}
			$left_trimmed ++;
		}
		if ($debug) {
			print STDERR join("\t", $to_left_trim, $left_trimmed, @alleles), "\n";	
		}
	}

	if ($left_trimmed) {
		$var->{POS} += $left_trimmed;
		$var->{REF} = shift @alleles;
		$var->{ALT} = ref $var->{ALT} eq 'ARRAY' ? \@alleles : $alleles[0];
	}
	return $var;
}

=head2 right_trim_or_left_extend VAR

Right trims ot left extend a variant.

=cut

sub right_trim_or_left_extend {
	my ($var, $sq, $debug) = @_;

	#$debug = 1;
	my @alleles = _unpack_alt($var);

	my ($to_right_trim, $to_left_extend) = (1, 0);
	my ($right_trimmed, $left_extended) = (0, 0);

	while ($to_right_trim || $to_left_extend) {
		
		($to_right_trim, $to_left_extend) = (1, 0);

		for(my $ii = 0; $ii < @alleles; $ii ++) {
			# Skip the star allele
			next if $alleles[$ii] eq '*';
			if ($alleles[$ii] ne '') {
				if (substr($alleles[$ii], -1) ne substr($alleles[0], -1)) {
					$to_right_trim = 0;
					# do not break here, you need to check for empty alleles that might exist
				}

				#left_extend will be taken care of later
				#if (substr($alleles[$ii], 0, 1) ne substr($alleles[0], 0, 1)) {
				#	$to_left_extend = 1;
				#}

				if ($var->{POS} == 1 && length($alleles[$ii]) == 1) {
					$to_right_trim = 0;
					last;
				}
			}
			else {
				$to_right_trim = 0;
				$to_left_extend = 1;
				last;
			}
		}

		if ($to_right_trim) {
			for(my $ii = 0; $ii < @alleles; $ii ++) {
				next if $alleles[$ii] eq '*';
				substr($alleles[$ii], -1, 1, "");
			}
			$right_trimmed ++;
		}

		if ($to_left_extend) {
			$left_extended ++;
			my $pos = $var->{POS}-$left_extended;
			my $refbase;
			if ($sq->exists($var->{CHROM}, $pos)){
				$refbase = $sq->get_slice($var->{CHROM}, $pos, $pos);
			}
			else {
				croak "Cannot fetch sequence at $var->{CHROM}:$pos";
			}
			for(my $ii = 0; $ii < @alleles; $ii ++) {
				$alleles[$ii] = $refbase.$alleles[$ii];
			}
		}
		if ($debug) {
			print STDERR join("\t", $to_right_trim, $to_left_extend, $right_trimmed, $left_extended, @alleles), "\n";
		}
	}

	if ($left_extended || $right_trimmed) {
		$var->{POS} -= $left_extended;
		$var->{REF} = shift @alleles;
		$var->{ALT} = ref $var->{ALT} eq 'ARRAY' ? \@alleles : $alleles[0];
	}
	return $var;
}




1;
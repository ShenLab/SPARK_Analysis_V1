package Variants;

use strict;
use warnings;
use Carp;
use Data::Dumper;
use List::Util qw|sum min max|;
use List::MoreUtils qw|all any uniq|;
use Perl6::Slurp;
use Tie::IxHash;
use Set::IntSpan;
use Config::Std;
use Regexp::Common qw|number|;
use Sort::Versions;
use Utils::Number qw|commafy|;
use Utils::Hash qw|merge_opts str2hash peek_hash|;
use Utils::Stat qw|mean median std mad weightedmean|;
use Utils::List qw|all_combs insert_after insert_before|;
use Utils::Seq qw|rev_comp extract_fasta|;
use Utils::File qw|open_file|;
use Utils::File::Iter qw|iter_file|;
use Genet::Ped;
use Genet::Var qw|normalize|;
use Genet::File::VCF qw|parse_vcf_header|;
use FaSlice;
use Genome::UCSC qw|hg_chrom|;
use Genome::Ranges qw|validate_elem spec_to_range|;
use Genome::Ranges::IntSpan;
use Genome::UCSC::TwoBit;

use base qw|Exporter|;

use FindBin qw|$Bin|;
use lib "$Bin/../lib";
use Shared qw|find_proberng parse_fstr|;


our @EXPORT_OK = qw|var_type is_cpg var_dist parse_nearby cluster_vars sort_vars
					get_vcfped_samps get_fieldsinfo get_fields 
					expand_site flatten_geno arrange_genodat expand_fam
					cnv_type cnv_id cnvid_to_dat cnvr_region cnv_flanklen cnv_flank_region expand_cnvfam
					%eff xfactor mutrater parse_operations anno_fields
					|;

# Functional consequences and their ranks are given in this web page
# https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences
# In the ranks below, we modified default ranking so the synonymous has higher priority than splice_region.

our %eff = (
	transcript_ablation =>	['HIGH',	1],
	splice_acceptor 	=>	['HIGH',	2],
	splice_donor 		=> 	['HIGH',	3],
	stop_gained			=> 	['HIGH',	4],
	frameshift 			=>	['HIGH',	5],
	stop_lost 			=>	['HIGH',	6],
	start_lost 			=> 	['HIGH',	7],
	transcript_amplification =>	['HIGH',	8],	
	inframe_insertion 	=> 	['MODERATE',	9],
	inframe_deletion 	=>	['MODERATE',	10],
	missense			=>	['MODERATE',	11],
	protein_altering 	=>	['MODERATE',	12],	
	synonymous			=>	['LOW',		13],
	incomplete_terminal_codon	=> ['LOW',	14],	
	start_retained		=> 	['LOW',	15,	'Silent'],
	stop_retained		=>  ['LOW',	16,	'Silent'], 
	splice_region		=>  ['LOW',	17],	
	coding_sequence		=> 	['MODIFIER',	18],	 ## <-- default cutoff, but can be adjusted
	mature_miRNA		=> 	['MODIFIER',	19],	
	'5_prime_UTR'		=>	['MODIFIER',	20],	
	'3_prime_UTR'		=>	['MODIFIER',	21],	
	non_coding_transcript_exon	=>	['MODIFIER',	22],	
	intron				=>	['MODIFIER',	23],	
	NMD_transcript		=> 	['MODIFIER',	24],	
	non_coding_transcript  => ['MODIFIER',	25],	
	upstream_gene		=> 	['MODIFIER',	26],	
	downstream_gene		=>	['MODIFIER',	27],	
	TFBS_ablation		=> 	['MODIFIER',	28],	
	TFBS_amplification	=>	['MODIFIER',	29],	
	TF_binding_site		=>	['MODIFIER',	30],	
	regulatory_region_ablation	=>	['MODERATE',	31],	
	regulatory_region_amplification	=>	['MODIFIER',	32],	
	feature_elongation	=> 	['MODIFIER',	33],	
	regulatory_region 	=>	['MODIFIER',	34],	
	feature_truncation	=>	['MODIFIER',	35],	
	intergenic			=>	['MODIFIER',	36]	);

=head1 SYNOPSIS

Utility functions related to analyze variants.

=head2 var_type REF, ALT

Determine the variant type. SNV or Indel (including MNV).

=cut

# Determine the variant type
sub var_type {
	my ($ref, $alt, $nostrict) = @_;
	if ($ref =~ /^[acgtACGT]$/ && $alt =~ /^[acgtACGT]$/) {
		return 'SNV';
	}
	elsif ($ref =~ /^[acgtACGT]+$/ && $alt =~ /^[acgtACGT]+$/) {
		if (length($ref) == length($alt) && length($ref) > 1) {
			return 'MNV';
		}
		else {
			return 'Indel';
		}	
	}
	else {
		if ($nostrict) {
			warn "Unknown variant type: $ref > $alt";
			return 'Unknown';
		}
		else {
			croak "Unknown variant type: $ref > $alt";
		}
	}
}


sub is_cpg {
	my ($context) = @_;
	if ($context =~ /CG([ACTG])>CA\1/ || $context =~ /([ACGT])CG>\1TG/) {
		return 1;
	}
	else {
		return 0;
	}
}

=head2 var_dist POS1,REF1,POS2,REF2

Calculate the distance between two variants.
Distance can be <=0 if two variants overlap each other.

=cut

sub var_dist {
	my ($pos1, $ref1, $pos2, $ref2) = @_;
	if ($pos1 > $pos2) {
		($pos1, $pos2) = ($pos2, $pos1);
		($ref1, $ref2) = ($ref2, $ref1);
	}
	my $end1 = $pos1+length($ref1)-1;
	my $st2 = $pos2;
	return $st2 - $end1;
}


=head2 parse_nearby STRING, QUAL

Parse nearby variants from de_recall. Format: "12_124352_AG_A[Q40]".

=cut

sub parse_nearby {
	my ($string, $qthres) = @_;
	my @vars;
	foreach my $var (split(',', $string)) {
		my ($chr, $pos, $ref, $alt, $qual) = ($var =~ /^((?:chr)?[0-9XY]+)_(\d+)_([ACGT]+)_([ACGT]+)\[Q(\d+)\]$/);
		#print join("\t", $chr, $pos, $ref, $alt, $qual), "\n";
		unless(defined $chr && defined $pos && defined $ref && defined $alt && defined $qual) {
			warn "Cannot parse variant $var";
			next;
		}
		if (defined $qthres) {
			if ($qual >= $qthres) {
				push @vars, [$chr, $pos, $ref, $alt];
			}
		}
		else {
			push @vars, [$chr, $pos, $ref, $alt, $qual];
		}
	}
	return @vars;
}


=head1 merge_vars REFSEQ,VAR1,VAR2,...

Merge nearby variant cluster.

Assume variants have no overlap with each other.

Variants will be presented by [CRHOM,POS,REF,ALT] array.

=cut

# Merge nearby variant clusters to MNVs
sub cluster_vars {
	# We will first order variants by chrom and position
	# extract the left-most and right-most positions
	# then substitute the reference sequence with alt alleles from right to left
	my $sq = shift @_;
	my @vars = @_;
	if (@vars <= 1) {
		warn "Less than two variants in the cluster!";
		return $vars[0];
	}
	else {
		# All variants must be on the same chromosome
		unless(all { $vars[0][0] eq $_->[0] } @vars) {
			die "Not all variants in the cluster are on the same chromosome";
		}
		my @varsord = sort { $a->[1] <=> $b->[1] } @vars;
		my $refal = $sq->get_slice($varsord[0][0], $varsord[0][1], $varsord[-1][1]+length($varsord[-1][2])-1);
		my %mnv = (CHROM => $varsord[0][0], POS => $varsord[0][1], REF => $refal);
		foreach my $var (reverse @varsord) {
			my $relpos = $var->[1]-$mnv{POS};
			substr($refal, $relpos, length($var->[2]), $var->[3]);
		}
		$mnv{ALT} = $refal;
		# normalize if necessary
		normalize(\%mnv, $sq);
		if (wantarray) {
			return @mnv{qw|CHROM POS REF ALT|};
		}
		else {
			return \%mnv;
		}	
	}
}



=head1 sort_vars VAR1, VAR2 ...

Sort variants by chrom, position, ref, allele.

Variant list can be an array of array refs with [Chr, Pos, Ref, Alt].
It can also be an array of variant IDs "Chr:Pos:Ref:Alt".

=cut

sub sort_vars {
	my @vars;
	if (ref $_[0] eq 'ARRAY') {
		unless(all { scalar(@$_) == 4 } @_) {
			die "Not all variants have four fields: chr, pos, ref, alt!";
		}
		@vars = sort { versioncmp($a->[0], $b->[0]) || $a->[1] <=> $b->[1] || 
					   uc($a->[2]) cmp uc($b->[2]) || uc($a->[3]) cmp uc($b->[3]) } @_;
	}
	else {
		unless(all { $_ =~ /^\w+:\d+:\w+:\w+$/ } @_) {
			die "Not all variants are chr:pos:ref:alt!";
		}
		@vars = sort { versioncmp($a->[0], $b->[0]) || $a->[1] <=> $b->[1] || 
					   uc($a->[2]) cmp uc($b->[2]) || uc($a->[3]) cmp uc($b->[3]) } map { [split(':')] } @_;
	}
	return @vars;
}


=head2 get_fieldsinfo VCF, OVERRIDE

Parse INFO and FORMAT number string for each field from VCF header.

=cut

# get site and geno fields from VCF header
sub get_fieldsinfo {
	my ($vcf, $site, $geno) = @_;
	my ($header, $info, $sampids) = parse_vcf_header($vcf);

	#my %sfields = map { $_ => $info->{INFO}{$_}{Number} } keys %{$info->{INFO}};
	tie my %sfields, 'Tie::IxHash';
	$sfields{CHROM} = 1;
	$sfields{POS} = 1;
	$sfields{REF} = 1;
	$sfields{ID} = '.';
	$sfields{ALT} = 'A';
	$sfields{QUAL} = '.';
	$sfields{FILTER} = '.';
	foreach my $field (keys %{$info->{INFO}}) {
		$sfields{$field} = $info->{INFO}{$field}{Number};
	}

	if (defined $site) {
		print STDERR "Appending field info for Site\n";
		my %stype = str2hash($site, { psep => ',', kvsep => ':'});
		foreach my $field (keys %stype) {
			my $origfield = (split(q|\.|, $field))[0];
			$sfields{$origfield} = $stype{$origfield} // do croak "Cannot find field type for $origfield";
		}
	}

	#my %gfields = map { $_ => $info->{FORMAT}{$_}{Number} } keys %{$info->{FORMAT}};
	tie my %gfields, 'Tie::IxHash';
	foreach my $field (keys %{$info->{FORMAT}}) {
		$gfields{$field} = $info->{FORMAT}{$field}{Number};
	}

	if (defined $geno) {
		print STDERR "Appending field info for Geno\n";
		my %gtype = str2hash($geno, { psep => ',', kvsep => ':'});
		foreach my $field (keys %gtype) {
			my $origfield = (split(q|\.|, $field))[0];
			$gfields{$origfield} = $gtype{$origfield} // do croak "Cannot find field type for $origfield";
		}
	}

	# Validity check
	while(my ($field, $numchar) = each %sfields) {
		unless ($numchar =~ /^[\.ARG]$/ || $numchar =~ /^\d+$/ ) {
			die "Cannot recognize number char $numchar for INFO field $field\n";
		}
	}
	while(my ($field, $numchar) = each %gfields) {
		unless ($numchar =~ /^[\.ARG]$/ || $numchar =~ /^\d+$/ ) {
			die "Cannot recognize number char $numchar for GENO field $field\n";
		}
	}

	return (\%sfields, \%gfields);
}


=head2 get_vcfped_samps VCF, PED, [PATTERN]

Get a list of samples that are common to both VCF and PED. For each sample
we include information from six standard columns from PED file plus optional
duplication/MZ twins information. 

Under the "trios" mode, it will output offspring within a list of complete trios.
Under this mode, duplicates in parents will be lost.

=cut

sub get_vcfped_samps {
	my ($vcffile, $pedfile, $argref) = @_;
	my $arg = merge_opts($argref, trios => undef, twins => undef, ignore => qr/_Re(\d*)$/);

	my $vcf = Genet::File::VCF->new($vcffile);

	my %samps;
	my $ped;
	my $pattern;
	if (defined $arg->{ignore}) {
		my $str = $arg->{ignore};
		$str =~ s/^['"]//; $str =~ s/['"]$//;
		$pattern = qr/$str/;
	}
	$ped = Genet::Ped->new($pedfile, { ignore => $pattern });
	
	my %vcfsamp = map { $_ => 1 } $vcf->get_sampids();
	if ($arg->{trios}) {
		foreach my $famid ($ped->get_famids()) {
			foreach my $trio ($ped->get_trios($famid)) {
				if (all { defined $vcfsamp{$_} } @$trio) {
					my $sex = $ped->get_sex($famid, $trio->[0]) // '0';
					my $phe = $ped->get_aff($famid, $trio->[0]) // '0';
					my $gender = $sex eq '1' ? 'Male' : $sex eq '2' ? 'Female' : 'NA';
					my $pheno = $phe eq '1' ? 'Unaffected' : $phe eq '2' ? 'Affected' : 'NA';
					$samps{$trio->[0]} = { FID => $famid, IID => $trio->[0],
										   DAD => $trio->[1], MOM => $trio->[2],
										   GENDER => $gender, PHENO => $pheno };
				}
			}
		}
	}
	else {
		foreach my $famid ($ped->get_famids()) {
			foreach my $iid ($ped->get_members($famid)) {
				next unless defined $vcfsamp{$iid};
				my $sex = $ped->get_sex($famid, $iid) // '0';
				my $phe = $ped->get_aff($famid, $iid) // '0';
				my $gender = $sex eq '1' ? 'Male' : $sex eq '2' ? 'Female' : 'NA';
				my $pheno = $phe eq '1' ? 'Unaffected' : $phe eq '2' ? 'Affected' : 'NA';
				$samps{$iid} = { FID => $famid, IID => $iid, 
								 DAD => $ped->get_father($famid, $iid) // '0',
								 MOM => $ped->get_mother($famid, $iid) // '0',
								 GENDER => $gender, PHENO => $pheno };
			}
		}
	}
	
	if (defined $pattern) {
		foreach my $iid (map { (split)[1] } slurp $pedfile) {
			next unless defined $vcfsamp{$iid};
			if ($iid =~ /$pattern/) {
				(my $origid = $iid) =~ s/$pattern//;
				if (defined $samps{$origid}) {
					push @{$samps{$origid}{DUP}}, $iid;
				}
			}
		}
	}
	if (defined $arg->{twins}) {
		my %twins;
		open my $fin, $arg->{twins} or die "Cannot open twins list: $arg->{twins}";
		while(<$fin>) {
			my @pair = split;
				unless(@pair == 2) {
				die "Incorrect number of samples in line: $_"
			}
			if (defined $twins{$pair[0],$pair[1]}) {
				warn "Twin pair $pair[0],$pair[1] has already been defined!";
				next;
			}
			$twins{$pair[0],$pair[1]} = 1;
			$twins{$pair[1],$pair[0]} = 1;
			# Only add to DUP entry for exsiting samples
			if (defined $samps{$pair[0]} && defined $samps{$pair[1]}) {
				push @{$samps{$pair[0]}{DUP}}, $pair[1];
				push @{$samps{$pair[1]}{DUP}}, $pair[0];
			}
		}
	}
	return \%samps;
}


=head2 expand_site SITE, FIELDS

  "Expand" the "site" data for each biallelic combination of ref and alt allele.

  FIELDS should be hash that specifies field name and a "number character" parsed from
  the VCF file header (1, A, R, G, .). For G type, we will assume ploidy=2 (polyploidy not implemented).

  Array should store information one allele per element
  Fields with missing data should be filled in with NA string '.'
  Return the regularized site data for all non-* alleles in a site.
  FILTER field will not be altered.

=cut

sub expand_site {
	my ($site, $sfields, $argref) = @_;
	my $arg = merge_opts($argref, nastring => '.', ploidy => 2);

	# Number of alternative alleles
	my $nalt = @{$site->{ALT}};

	# Move info field to upper level
	while(my ($field, $type) = each %$sfields) {
		if (defined $site->{INFO}{$field}) {
			if (defined $site->{$field}) {
				die "INFO field $field is in conflict with other builtin field names" ;
			}
			$site->{$field} = $site->{INFO}{$field};
		}
	}
	_validate_fields($site, $sfields, $nalt);
	
	# Now expand site to sites
	my @vars;
	for(my $ii = 0; $ii < $nalt; $ii ++) {
		# * allele should not be skipped
		#next if $site->{ALT}[$ii] eq '*';
		push @vars => _expand_allele_ii($site, $sfields, $ii, $nalt, $arg->{nastring});
	}
	return @vars;
}

# Helper function
sub _validate_fields {
	my ($site, $fields, $nalt, $extra) = @_;
	while(my ($field, $type) = each %$fields) {
		next unless exists $site->{$field} && defined $site->{$field};
		if ($type eq 'A') {
			if ($nalt > 1) {
				unless (ref $site->{$field} eq 'ARRAY' && @{$site->{$field}} == $nalt) {
					print $extra if defined $extra;
					print Dumper $site;
					croak "Field $field is of type A but does not have $nalt elements for a multi-allic site";
				}
			}
			else {
				if (ref $site->{$field} eq 'ARRAY') {
					unless (@{$site->{$field}} == 1) {
						print $extra if defined $extra;
						print Dumper $site;
						croak "Field $field is of type A but does not have 1 element for a bi-allic site";
					}	
				}
			}
		}
		elsif ($type eq 'R') {
			unless (ref $site->{$field} eq 'ARRAY' && @{$site->{$field}} == $nalt + 1) {
				print $extra if defined $extra;
				print Dumper $site;
				croak "Field $field is of type R but does not have $nalt+1 elements";
			}
		}
		elsif ($type eq 'G') {
			unless (ref $site->{$field} eq 'ARRAY' && @{$site->{$field}} == ($nalt+1)*($nalt+2)/2) {
				print $extra if defined $extra;
				print Dumper $site;
				croak "Field $field is of type G but does not have ($nalt+1)*($nalt+2)/2 elements";
			}
		}
		elsif ($type =~ /^(\d+)$/) {
			my $nel = $1;
			if ($nel == 1) {
				if (ref $site->{$field} eq 'ARRAY') {
					print $extra if defined $extra;
					print Dumper $site;
					croak "Field $field is expected to be a scalar";
				}	
			}
			elsif ($nel > 1) {
				unless (ref $site->{$field} eq 'ARRAY' && @{$site->{$field}} == $nel) {
					print $extra if defined $extra;
					print Dumper $site;
					croak "Field $field is of type $nel but does not have $nel elements";
				}
			}
			elsif ($nel == 0) {
				unless (exists $site->{$field} && defined $site->{$field}) {
					print Dumper $site;
					croak "Field $field is of type $nel but not consistent with this type";
				}
			}
			else {
				croak "Incorrect field type for $field: $nel";
			}
		}
	}
	return $site;
}

# Assume site data are validated
# index ii in this function is 0-based
sub _expand_allele_ii {
	my ($site, $fields, $ii, $nalt, $nastring) = @_;
	my %dat;
	while(my ($field, $type) = each %$fields) {
		if ($type eq 'A') {
			if (defined $site->{$field}) {
				if (ref $site->{$field} eq 'ARRAY') {
					$dat{$field} = $site->{$field}[$ii];
				}
				else {
					$dat{$field} = $site->{$field}
				}
			}
			else {
				$dat{$field} = $nastring;
			}
		}
		elsif ($type eq 'R') {
			$dat{$field.".REF"} = defined $site->{$field} ? $site->{$field}[0] : $nastring;
			$dat{$field.".ALT"} = defined $site->{$field} ? $site->{$field}[$ii+1] : $nastring;
		}
		elsif ($type =~ /^(\d+)$/) {
			my $nel = $1;
			if ($nel == 1) {
				$dat{$field} = $site->{$field} // $nastring;
			}
			elsif ($nel > 1) {
				for(my $jj = 1; $jj <= $nel; $jj ++) {
					if (defined $site->{$field}) {
						$dat{$field.".".$jj} = $site->{$field}[$jj-1];
					}
					else {
						$dat{$field.".".$jj} = $nastring;
					}
				}
			}
			elsif ($nel == 0) {
				$dat{$field} = exists $site->{$field} ? 1 : 0; 
			}
			else {
				die "Incorrect number of element $nel for site field $field";
			}
		}
		elsif ($type eq 'G') {
			$dat{$field.".HOMREF"} = defined $site->{$field} ? $site->{$field}[0] : $nastring;
			# for P=2, the index of the genotype "a/b" where a<=b is b(b+1)/2+a
			# E.g. for P=2 and N=2, the ordering is 00,01,11,02,12,22	
			# b = $ii + 1, a = 0 or $ii + 1
			$dat{$field.".HET"} = defined $site->{$field} ? $site->{$field}[($ii+1)*($ii+2)/2] : $nastring;
			$dat{$field.".HOMALT"} = defined $site->{$field} ? $site->{$field}[($ii+1)*($ii+2)/2+($ii+1)] : $nastring;
		}
		else {
			if (defined $site->{$field}) {
				if (ref $site->{$field} eq 'ARRAY') {
					if ($field eq 'ID' || $field eq 'FILTER') {
						$dat{$field} = join(";", @{$site->{$field}});
					}
					else {
						$dat{$field} = join(",", @{$site->{$field}});
					}
				}
				else {
					$dat{$field} = $site->{$field};
				}
			}
			else {
				$dat{$field} = $nastring;
			}	
		}
	}
	return \%dat;
}

=head2 flatten_geno GDAT, FIELDS, JJ, NALT

  Extract allele specific information for genotype
  NOTE: JJ is JJ-th alternative allele, 1-based.

=cut

sub flatten_geno {
	my ($gdat, $gfields, $jj, $nalt, $argref) = @_;
	my $arg = merge_opts($argref, nastring => '.', ploidy => 2, extra => undef);
	_validate_fields($gdat, $gfields, $nalt, $arg->{extra});
	return _expand_allele_ii($gdat, $gfields, $jj-1, $nalt, $arg->{nastring});
}


# Help function in reformating genotype output

sub arrange_genodat {
	my ($sampdat, $fields, $nachar) = @_;
	my @output;
	foreach my $gf (keys %$fields) {
		my @gfdat;
		foreach my $gdat (@$sampdat) {
			push @gfdat, $gdat->{$gf} // $nachar;
		}
		push @output, \@gfdat;
	}
	return @output;
}

#
# Expand genotypes of family members when parsing var_filter output
#
# We will also pack "field" from original sample to "_field"
sub expand_fam {
	my ($dat, $famrm) = @_;
	# Determine the fields for expansion	
	my @fields = grep { /\.FamMembers$/ } sort keys %$dat;
	# Expand the data, and test for equal length
	my %expand;
	foreach my $field (qw|FamMembers Relations Phenotypes|, @fields) {
		unless(defined $dat->{$field}) {
			die "Cannot find field $field when expanding family data";
		}
		if ($dat->{$field} eq "") { # <-- fix the empty field to avoid zero expansion
			$dat->{$field} = '.';
		}
		$expand{$field} = [split(',', $dat->{$field})];
	}
	my $n_samp = scalar(@{$expand{FamMembers}});
	unless(all { scalar(@$_) == $n_samp } values %expand) {
		print STDERR Dumper $dat;
		die "Not all expanded fields have the same length as the number of family members: $n_samp";
	}
	# Re-pack the expanded data
	my @data;
	for(my $ii = 0; $ii < $n_samp; $ii ++) {
		# Removing family members that appear in the family exclusion list
		next if defined $famrm && defined $famrm->{$expand{FamMembers}[$ii]};
		my %data = (IID => $expand{FamMembers}[$ii],
					Relation => $expand{Relations}[$ii],
					Pheno => $expand{Phenotypes}[$ii] );
		foreach my $field (@fields) {
			(my $newfd = $field) =~ s/\.FamMembers$//;
			$data{$newfd} = $expand{$field}[$ii];
		}
		# Adding fields from original sample
		while(my ($key, $val) = each %$dat) {
			next if defined $expand{$key};
			if (defined $data{$key}) {
				if (defined $data{"_$key"}) {
					print Dumper $dat;
					die "Hidden field _$key has been defined!";
				}
				$data{"_$key"} = $val;
			}
			else {
				$data{$key} = $val;
			}
		}
		push @data, \%data;
	}
	return @data;
}

# Expand family info for CNV data
sub expand_cnvfam {
	my ($dat, $prefix, $famrm) = @_;
	$prefix = "FamMember" unless defined $prefix;
	# Determine the fields for expansion	
	my @fields = grep { /^${prefix}_/ && !/^${prefix}_Carrier/ } sort keys %$dat;
	# Expand the data, and test for equal length
	my %expand;
	foreach my $field (@fields) {
		unless(defined $dat->{$field}) {
			die "Cannot find field $field when expanding family data";
		}
		if ($dat->{$field} eq "") {
			$dat->{$field} = ".";
		}
		$expand{$field} = [split(',', $dat->{$field})];
	}
	my $n_samp = scalar(@{$expand{$fields[0]}});
	unless(all { scalar(@$_) == $n_samp } values %expand) {
		die "Not all expanded fields have the same length as the number of family members";
	}
	# Re-pack the expanded data
	my @data;
	for(my $ii = 0; $ii < $n_samp; $ii ++) {
		# Removing family members that appear in the family exclusion list
		next if defined $famrm && defined $famrm->{$expand{"${prefix}_IDs"}[$ii]};
		my %data;
		foreach my $field (@fields) {
			(my $newfd = $field) =~ s/^${prefix}_//;
			$newfd =~ s/s$//; # Also strip off the plural suffix
			$data{$newfd} = $expand{$field}[$ii];
		}
		# Adding fields from original sample
		while(my ($key, $val) = each %$dat) {
			next if defined $expand{$key};
			if (defined $data{$key}) {
				if (defined $data{"_$key"}) {
					print Dumper $dat;
					die "Hidden field _$key has been defined!";
				}
				$data{"_$key"} = $val;
			}
			else {
				$data{$key} = $val;
			}
		}
		push @data, \%data;
	}
	return @data;
}


=head2 cnv_type CN, OPTIONS

Return standardized CNV type: Dup/Del/Normal

We do not deal with absolute copy number, but only classify CNV type. 
CN=2 will be treated specically with the following rules
1. When no other option is give, CN=2 will be taken as duplication. This is based on assumption
  that we are deal with male chrX CNVs and normal copy number should not be reported in autosomes.
2. When maleX option is on, CN=2 will only be duplication whereas CN=1 will be taken as deletion.
3. For CNs that cannot be classified as del or dup, default is set them to normal.
  If asis is set, they will appear in original format.

There is an further option strict, when switched on, Normal type CNV is not allowed.

=cut

sub cnv_type {
	my ($cn, $argref) = @_;
	my $arg = merge_opts($argref, strict => undef, maleX => undef, asis => undef);
	my $type;
	if ($cn =~ /^(dup|gain|amp)/i || ($cn =~ /^\d+$/ && $cn >= 2)) {
		if ($cn eq '2') {
			if ($arg->{maleX}) {
				$type = 'Dup';
			}
			else {
				if ($arg->{strict}) {
					die "Normal copy is not allowed under strict mode of cnv_type!";
				}
				else {
					#unless ($arg->{asis}) {
						$type = 'Normal';
					#}
				}
			}
		}
		else {
			$type = 'Dup';
		}
	} 
	elsif ($cn =~ /^(del|loss)/i || ($cn =~ /^\d+$/ && $cn < 2)) {
		if ($cn eq '1') {
			if ($arg->{maleX}) {
				if ($arg->{strict}) {
					die "Normal copy of male chrX is not allowed under strict mode of cnv_type";
				}
				else {
					#unless ($arg->{asis}) {
						$type = 'Normal';
					#}
				}
			}
			else {
				$type = 'Del';
			}
		}
		else {
			$type = 'Del';
		}	
	}
	else {
		if ($arg->{strict}) {
			die "Cannot determine CNV type: $cn";
		}
		else {
			unless ($arg->{asis}) {
				warn "Cannot determine CNV type: $cn... will assume \"Normal\"";
				$type = 'Normal';
			}
			
		}
	}
	return $type;
}


=head2 cnv_id CHR, START, END, CN

Create CNV ID: ChromRegionSpec_Type

=cut

sub cnv_id {
	my ($chr, $start, $end, $cn, $argref) = @_;
	return $chr.":".commafy($start)."-".commafy($end)."_".cnv_type($cn, $argref);
}

# Unpack CNV ID to region and type
# Note that we do not explicitly check for CNType, so it can also be use for other SVs.
sub cnvid_to_dat {
	my ($cnvid) = @_;
	my ($spec, $cntype) = split('_', $cnvid);
	my ($chr, $start, $end) = spec_to_range($spec);
	unless(defined $cntype && defined $chr && defined $start && defined $end) {
		die "Cannot parse CNVID: $cnvid";
	}
	unless($start < $end) {
		die "Incorrect range for CNV: $cnvid";
	}
	return ($chr, $start, $end, $cntype);
}

# Create CNVR from a list of overlapping CNVs 
sub cnvr_region {
	my @cnvids = @_;
	my ($chrom, @spans, @cntypes);
	foreach my $cnvid (@cnvids) {
		my ($chr, $start, $end, $cntype) = cnvid_to_dat($cnvid);
		unless(defined $chrom) {
			$chrom = $chr;
		}
		else {
			unless($chrom eq $chr) {
				die "CNVs in one CNVR are foun on different chromosomes";
			}
		}
		push @spans, [$start, $end];
		push @cntypes, $cntype;
	}

	# CNVR is represented as an Intspan object
	my $region = Set::IntSpan->new(\@spans);
	my @ranges = $region->spans();
	if (@ranges > 1) {
		print Dumper \@ranges;
		die "cnvr_region: CNVR is disjoint";
	}
	elsif (@ranges == 0) {
		die "cnvr_region: CNVR is empty";
	}

	my @cntype = uniq sort @cntypes;
	if (@cntype == 1) {
		return ($chrom, $ranges[0][0], $ranges[0][1], $cntype[0]);
	}
	else {
		return ($chrom, $ranges[0][0], $ranges[0][1], "Err");
	}
}



=head2 cnv_flanklen 

Determine the default flanking region length for CNV plotting.

=cut

sub cnv_flanklen {
	my ($length) = @_;
	if ($length >= 500_000) {
		return int(0.5*$length);
  	}
  	elsif ($length <= 50_000) {
  		return int(5*$length);
	}
  	else {
  		# 50~500k, ratio of 500k/length should be 1~10
		return int((0.5*500_000/$length)*$length);
  	}
}

=head2 cnv_flanklen_region

Determine the flanking region for CNV plotting.
This is more flexible than the above convenient function, allowing the use of probe file
and specify parameters to lower and upper bound of sizes, and fractions of flanking sizes.
If probes are provided, size will be measured in number of probes.

If a CNV has no probe coverage, it will return an empty region.

=cut

sub cnv_flank_region {
	my ($chrom, $start, $end, $argref) = @_;
	my $arg = merge_opts($argref, probes => undef, 
						sizeupper => 500_000, sizelower => 50_000,
						ratioupper => 0.5, ratiolower => 5);

	unless($arg->{sizeupper} > $arg->{sizelower} && $arg->{sizelower} > 0 &&
		   $arg->{ratioupper} < $arg->{ratiolower} &&  $arg->{ratioupper} > 0) {
		die "Incorrect size range and ratios at cnv_flank_region";
	}

	my ($length, $flank, $chrom_probe, $lo_pos, $hi_pos);
	if ($arg->{probes}) {
		my $firstpbchr = peek_hash($arg->{probes});
		my $probechr;
		if ($firstpbchr =~ /^chr/) {
			$probechr = 1;
		}
		else {
			$probechr = 0;
		}
		
		$chrom_probe = $chrom;
		if ($chrom =~ /^chr/ && $probechr == 0) {
			$chrom_probe =~ s/^chr//; 
		}
		elsif ($chrom !~ /^chr/ && $probechr == 1) {
			$chrom_probe = hg_chrom($chrom_probe);
		}

		($lo_pos, $hi_pos) = find_proberng($arg->{probes}, $chrom_probe, $start, $end);
		unless(defined $lo_pos && defined $hi_pos) {
			return;
		}
		$length = $hi_pos - $lo_pos + 1;
	}
	else {
		$length = $end - $start + 1;
	}

	# Determine flank length
	if ($length >= $arg->{sizeupper}) {
		$flank = int($length*$arg->{ratioupper});
	}
	elsif ($length <= $arg->{sizelower}) {
		$flank = int($length*$arg->{ratiolower});
	}
	else {
		my $ratio = $arg->{ratioupper} *  ($arg->{sizeupper}/$length) *
			($arg->{ratiolower}/$arg->{ratioupper})/($arg->{sizeupper}/$arg->{sizelower});
		$flank = int($length*$ratio);
	}
	
	# Determine the flanking region
	if ($arg->{probes}) {
		my $st_pos = max(0, $lo_pos-$flank);
		my $ed_pos = min($hi_pos+$flank, scalar(@{$arg->{probes}{$chrom_probe}{Ends}})-1);
		my @range = ($chrom, $arg->{probes}{$chrom_probe}{Starts}[$st_pos], $arg->{probes}{$chrom_probe}{Ends}[$ed_pos]);
		validate_elem(@range);
		
		if (wantarray) {
			return @range;
		}
		else {
			return sprintf("%s:%d-%d", @range);
		}

	}
	else {
		my @range = ($chrom, max(1,$start-$flank), $end+$flank);
		validate_elem(@range);
		if (wantarray) {
			return @range;
		}
		else {
			return sprintf("%s:%d-%d", @range);
		}
	}
}


=head2 xfactor

Adjustment factor for chrX mutation rate.

=cut

sub xfactor {
	my ($pmale, $alpha) = @_;
	unless($pmale >= 0 && $pmale <= 1) {
		die "Incorrect proportion of male sample: $pmale";
	}
	$alpha = 3.4 unless defined $alpha; # <- this is the DDD recommended alpha
	my $pfemale = 1-$pmale;
	my $autosomal = 2;
	my $female_transmissions = 1;
	my $male_transmissions = $pfemale;
	my $male_factor = 2 / (1 + (1 / $alpha));
	my $female_factor = 2 / (1 + $alpha);
	my $x_factor = (($male_transmissions * $male_factor) +
				    ($female_transmissions * $female_factor)) / $autosomal;
	return $x_factor;
}

=head2 mutrater 

Create a subroutine to calculate point mutation rate at each genomic position.
The subrountine will take input from a line of annotated variant, and return
a relative or absolute mutation rate for that position.

Arguments to the subroutine creater include:
	* Method -- Method used to calculate relative mutation rate (details below)
				It also determines the way to parse the lookup table.
	* Lookup -- Lookup table of context dependent mutation rate
	* Fasta -- Genome reference sequence in fasta format (or 2bit format)
	* Weight -- Optional, weight at each single base pair (e.g. to correct for coverage)
	* Scale -- the scaling factor(s) for the mutation rate (constant or bedgraph file). 
	* Chunk -- Fast mode, gives a chunk size to store sequence and weights in RAM.

The lookup table can either be local (3mer or 7mer context rate) or exome/genome 
(mutation rate associated with each Chrom,Position,Ref,Alt)

Currently, we have implemented two versions of local mutation rate based on sequence context.
The first method was tri-nucleotide sequence context model, we obtained local rates
from DenovoNear package ("3merDenovoNear"). The second method was 7mer context model
in which rates were derived from extremely rare variants from WGS data. We obtained
local rates from Mr-Eel package (Carlson et al. 2019). 

When using local context to look for mutation rates, genome reference sequence is used to look
for sequence context or we will make use of "Context" column to determine the sequence context.

More complex methods are available. For example, Francioli et al. 2015 and Carlson et al. 2019's 
predicted mutation rate based on their models. And DDD 2019's per-bp mutation rates at each 
position of the exome. To accomondate those complex models, we can also lookup per-bp mutation
rates from a tabix indexed database table. The DB-table gives position and allele specific mutation 
rates and should have five columns: Chrom,Position,Ref,Alt,Rate. We will use the rate from the 5th column.
If the position and allele specific mutation rate cannot be found in the lookup table, 
it will fall back to use local rate.

Scale can be either a constant or a bed graph file. For the latter case, we can get region-specific 
scaling factors. If some region does not exist in bedgraph, we will check if other regions from
the same chromosome exist. If so, we will fetch the factor from the nearest interval; otherwise
we will fall back to mean scaling factor across all the regions.

Weight file is also position-specific. It should be tabix index with three columns: Chrom, Position and Weight.
If a position cannot be found in the weight file, then it will have zero weight (e.g., no coverage).

The subroutine query context, mutation rate and weight for each SNV at any given genomic position 
one at a time. Under faster mode, we will assume all enumeration of SNVs from the input table are ordered by 
genomic positions. So we can slurp the sequence and weights for a genomic region and keep them in a lookup table 
for faster access.

Note: chromosome nomenclautures in fasta, rate/weight table must be consistent with SNV enumeration table!
No warning will be issued when chromosome nomencalutures do not match!

=cut

sub mutrater {
	my ($argref) = @_;
	my $arg = merge_opts($argref, Method => undef, Lookup => undef, Fasta => undef, Rate => undef,
								  Weight => undef, Scale => 1, Chunk  => 1000000, Count => undef); # 1000000

	unless(defined $arg->{Method}) {
		die "Must provide Method to determine motif length!";
	}
	unless(defined $arg->{Lookup} && -f $arg->{Lookup}) {
		die "Cannot find lookup table: $arg->{Lookup}!";
	}

	my $lookup = _slurp_lookup($arg->{Method}, $arg->{Lookup});
	if (defined $arg->{Count}) {
		if (defined $arg->{Rate}) {
			die "Per-base mutation rate is not used under Count mode!";
		}
		# Reset all rate to 1 under count mode
		foreach my $context (sort keys %$lookup) {
			$lookup->{$context} = 1;
		}
	}

	my ($motiflen) = ($arg->{Method} =~ /(\d+)mer/);
	unless(defined $motiflen) {
		warn "Cannot determine motif length from Method, fallback to default 3mer";
		$motiflen = 3;
	}
	unless($motiflen % 2 == 1) {
		die "Motif length must be an odd number";
	}
	my $adjlen = ($motiflen-1)/2;

	# Reference sequence slicer
	my $sq;
	if (defined $arg->{Fasta}) {
		if ($arg->{Fasta} =~ /\.2bit$/) {
			$sq = Genome::UCSC::TwoBit->new($arg->{Fasta});
		}
		else {
			$sq = FaSlice->new(file => $arg->{Fasta});
		}
	}

	# Weight and mutation rate at each positions
	my ($tabwt, $tabrate);
	if (defined $arg->{Weight} || defined $arg->{Rate}) {
		require Bio::DB::HTS::Tabix;
		if ($arg->{Weight}) {
			$tabwt = Bio::DB::HTS::Tabix->new(filename => $arg->{Weight});
		}
		if ($arg->{Rate}) {
			$tabrate = Bio::DB::HTS::Tabix->new(filename => $arg->{Rate});
		}
	}

	# Global or region-specific scaling factor
	# Average scale is used when target location is outside the bedgraph region 
	my ($scaler, $avgscale);
	if ($arg->{Scale} =~ /^$RE{num}{real}$/) {
		$avgscale = $arg->{Scale};
	}
	else {
		unless(-f $arg->{Scale}) {
			die "Cannot find region-specific scaling factor bedgraph: $arg->{Scale}";
		}
		my ($sum, $weight);
		my $fin = open_file($arg->{Scale});
		while(<$fin>) {
			my @dat = split;
			unless(@dat == 4) {
				print STDERR $_;
				die "Incorrect number of columns in the scale file";
			}
			my $len = $dat[2]-$dat[1];
			if ($len <= 0) {
				print STDERR $_;
				die "Incorrect region length in the scale file!";
			}
			$weight += $len;
			$sum += $dat[3] * $len;
		}
		$avgscale = $sum/$weight;
		my %opt;
		if ($arg->{Scale} =~ /\.(bed|bedgraph)$/ || $arg->{Scale} =~ /\.(bed|bedgraph)\.gz$/) {
			$opt{bed} = 1;
		}
		$scaler = Genome::Ranges::IntSpan->new($arg->{Scale}, \%opt);
	}

	# Define subroutine to calculate mutation rate for each SNV
	# If refseq is not provided, we will obtain a sequence context from annotation.
	# If sequence context contain N, we will calculate the average rate from
	# all possible context
	my $rater;

	my ($prev_chr, $prev_pos, $refseq) = ("", 1, "");
	my ($prev_name, $prev_loc, $weights, $rates) = ("", undef, {}, {});

	$rater = sub {
		unless($_[0]{Ref} =~ /^[ACGT]$/ && $_[0]{Alt} =~ /^[ACGT]$/) {
			die "Incorrect Ref/Alt allele: $_[0]{Ref}/$_[0]{Alt}";
		}

		if ($arg->{Fasta}) {
			if ($arg->{Chunk}) {
			# When Chunk is provided, we fetch sequence from [ Pos-AdjLen ~ Pos-AdjLen+(Chunk-1) ] <- 1-based, inclusive
			# $prev_pos is set to the start point of the range where sequence was collected
				if ($_[0]{Chrom} ne $prev_chr ||
					$_[0]{Position}+$adjlen > $prev_pos+$arg->{Chunk}-1) {
					my $start = max($_[0]{Position}-$adjlen, 1);
					my $end = min($_[0]{Position}-$adjlen+$arg->{Chunk}-1, $sq->{SEQLEN}{$_[0]{Chrom}});
					$refseq = $sq->get_slice($_[0]{Chrom}, $start, $end);
					$prev_chr = $_[0]{Chrom};
					$prev_pos = $_[0]{Position}-$adjlen;
				}
			}
			else {
			# When Chunk is not specified, we will fetch sequence from [ 1 ~ ChromSize ]
			#	if ($_[0]{Chrom} ne $prev_chr) {
			#		$refseq = extract_fasta($arg->{Fasta}, $_[0]{Chrom});
			#		$prev_chr = $_[0]{Chrom};
			#		$prev_pos = 1;
			#	}
			#}
				# When chunk is not specified, we will fetch local motif sequence
				my $start = $_[0]{Position}-$adjlen;
				my $end = $_[0]{Position}+$adjlen;
				# Start and end should be within chromosome range, we will not check this
				$refseq = $sq->get_slice($_[0]{Chrom}, $start, $end);
				$prev_chr = $_[0]{Chrom};
				$prev_pos = $_[0]{Position}-$adjlen;
			}
		}
		else {
			die "Cannot find sequence context" unless defined $_[0]{Context};
		}

		if ($arg->{Weight} || $arg->{Rate}) {
			if ($arg->{Chunk}) {		
				if ($_[0]{Chrom} ne $prev_name ||
					$_[0]{Position} > $prev_loc+$arg->{Chunk}-1) {
					if ($arg->{Weight}) {
						$weights = {};
						my $iter = $tabwt->query(sprintf("%s:%d-%d", $_[0]{Chrom}, $_[0]{Position}, $_[0]{Position}+$arg->{Chunk}-1));
						if (defined $iter) {
							while(my $line = $iter->next) {
								chomp($line);
								my ($chr, $pos, $wt) = split(/\t+/, $line);
								$weights->{$pos} = $wt;
							}
						}
					}
					if ($arg->{Rate}) {
						$rates = {};
						my $iter = $tabrate->query(sprintf("%s:%d-%d", $_[0]{Chrom}, $_[0]{Position}, $_[0]{Position}+$arg->{Chunk}-1));
						if (defined $iter) {
							while(my $line = $iter->next) {
								chomp($line);
								my ($chr, $pos, $ref, $alt, $rt) = split(/\t+/, $line);
								my $varid = join(":", $chr, $pos, $ref, $alt);
								$rates->{$varid} = $rt;
							}
						}
					}
					$prev_name = $_[0]{Chrom};
					$prev_loc = $_[0]{Position};
				}
			}
			else {
				# Without Chunk, collect weight/rate at a single position
				if ($arg->{Weight}) {
					$weights = {};
					my $iter = $tabwt->query(sprintf("%s:%d-%d", $_[0]{Chrom}, $_[0]{Position}, $_[0]{Position}));
					while(my $line = $iter->next) {
						chomp($line);
						my ($chr, $pos, $wt) = split(/\t+/, $line);
						$weights->{$pos} = $wt;
					}
				}
				if ($arg->{Rate}) {
					$rates = {};
					my $iter = $tabrate->query(sprintf("%s:%d-%d", $_[0]{Chrom}, $_[0]{Position}, $_[0]{Position}));
					while(my $line = $iter->next) {
						chomp($line);
						my ($chr, $pos, $ref, $alt, $rt) = split(/\t+/, $line);
						my $varid = join(":", $chr, $pos, $ref, $alt);
						$rates->{$varid} = $rt;
					}
				}
			}
		}

		my ($from, $to);
		if ($arg->{Fasta}) {
			$from = uc(substr($refseq, $_[0]{Position}-$prev_pos-$adjlen, $motiflen));
			$to = $from; substr($to, $adjlen, 1, $_[0]{Alt});
			# Context field will be updated if RefSeq is provided <-- don't update context in place!
			# $_[0]{Context} = "$from>$to";
		}
		else {
			($from, $to) = split('>', $_[0]{Context});
			# Validate sequence context and alleles
			unless(defined $from && length($from) == $motiflen &&
				   defined $to   && length($to)   == $motiflen) {
				die "Incorrect length of sequence context for $arg->{Method}: $_[0]{Context}"
			}
		}
		unless(substr($from, $adjlen, 1) eq $_[0]{Ref} && 
			substr($to, $adjlen, 1) eq $_[0]{Alt}) {
			die "Ref or Alt allele does not match sequence context!";
		}

		# Find out the final weight (scale factor * weight)
		my $fscale;
		if (defined $scaler) {
			$fscale = $scaler->lookup_nearest($_[0]{Chrom}, $_[0]{Position});
			$fscale = $avgscale unless defined $fscale;
		}
		else {
			$fscale = $avgscale;
		}
		my $weight;
		if ($arg->{Weight}) {
			# Set weight to 0 if not found in the weight file
			my $wt = $weights->{$_[0]{Position}} // 0;
			$weight = $wt * $fscale;
		}
		else {
			$weight = $fscale;
		}

		my $varid = join(":", @{$_[0]}{qw|Chrom Position Ref Alt|});
		if (defined $rates->{$varid}) {
			if (wantarray) {
				return ($rates->{$varid}*$weight, "$from>$to");
			}
			else {
				return $rates->{$varid}*$weight;
			}
		}

		if ($from =~ /^[ACGT]+$/ && $to =~ /^[ACGT]+$/) {
			my $rate = $lookup->{$from,$to} // do { die "Cannot find rate for $from>$to" };
			die "Cannot mutation rate for $from>$to" unless defined $rate;
			if (wantarray) {
				return ($rate*$weight, "$from>$to");
			}
			else {
				return $rate*$weight;
			}
		}
		elsif ($from =~ /^[ACGTN]+$/ && $to =~ /^[ACGTN]+$/) {
			# In the presence of N bases, we will average out all possible combinations in the context
			my (@froms, @tos);
			for(my $ii = 0; $ii < length($from); $ii ++) {
				my $currbp = substr($from, $ii, 1);
				if($currbp eq 'N') {
					push @froms, [qw|A C G T|];
				}
				else {
					push @froms, [$currbp];
				}
				if ($ii == $adjlen) {
					push @tos, [$_[0]{Alt}];
				}
				else {
					push @tos, $froms[$ii];
				}
			}
			my @comb_froms = all_combs(@froms);
			my @comb_tos = all_combs(@tos);
			my @rates;
			for(my $jj = 0; $jj < @comb_froms; $jj ++) {
				my $from_jj = join("", @{$comb_froms[$jj]});
				my $to_jj = join("", @{$comb_tos[$jj]});
				push @rates, $lookup->{$from_jj,$to_jj} // do { die "Cannot find rate for $from>$to" };
			}
			#print STDERR "Averaging over ambiguous allele\n";
			if (wantarray) {
				return (mean(@rates)*$weight, "$from>$to");
			}
			else {
				return mean(@rates)*$weight;
			}
		}
		else {
			die "Incorrect base pair in sequence context: $from>$to";
		}
	};

	return $rater;
}

# Helper function: slurp lookup table

sub _slurp_lookup {
	my ($method, $tabfile) = @_;
	my %lookup;
	if ($method eq '3merDenovoNear') {
		my $it = iter_file($tabfile);
		while(my $dat = $it->()) {
			$lookup{$dat->{from},$dat->{to}} = $dat->{mu_snp};
		}
	} 
	elsif ($method eq '7merMrEel') {
		# Mr-eel recommended scaling factor
		my $factor = 1.2e-8*3e9/36087319;
		my $it = iter_file($tabfile);
		while(my $dat = $it->()) {
			my @context = ($dat->{Sequence} =~ /([ACGT]{7})\(([ACGT]{7})\)/);
			unless(@context == 2) {
				die "Cannot find sequence context from both strands: $dat->{Sequence}!"
			}
			foreach my $from (@context) {
				my $ref = substr($from, 3, 1);
				foreach my $alt (qw|A G C T|) {
					next if $alt eq $ref;
					my $key;
					if ($ref eq 'A' || $ref eq 'G') {
						$key = $ref.rev_comp($ref)."_".$alt.rev_comp($alt);
					}
					else {
						$key = rev_comp($ref).$ref."_".rev_comp($alt).$alt;
					}
					if ($dat->{$key} == 0) {
						die "Zero mutation rate for $ref>$alt in $from context";
					}
					my $to = $from; 
					substr($to, 3, 1, $alt);
					$lookup{$from,$to} = $dat->{$key} * $factor;
				}
			}
		}
	}
	else {
		die "Cannot recognize method: $method";
	}
	return \%lookup;
}

=head2 parse_operations

Parse the operation string used to "group by" a list of data, and return operation call backs.

Example: count_distinct(FamID,IID):CarrierFreq;distinct(Syndrome):KnownPathoCNVs

Group by operations are defined in BEDtools package: https://bedtools.readthedocs.io/en/latest/content/tools/groupby.html

Note: Some operations like first/last depends on the sorting of input intervals which should
be treated separately.

=cut

sub parse_operations {
	my ($fstr, $argref) = @_;
	my $arg = merge_opts($argref, dbfields => undef, orsep => ',', dfsep => '|');
	unless(defined $arg->{dbfields}) {
		die "Must provide dbfields";
	}
	my $dbfields;
	if (ref $arg->{dbfields} eq 'HASH') {
		$dbfields = $arg->{dbfields};
	}
	elsif (ref $arg->{dbfields} eq 'ARRAY') {
		my %dbfds;
		@dbfds{@{$arg->{dbfields}}} = @{$arg->{dbfields}};
		$dbfields = \%dbfds;
	}
	else {
		$dbfields = parse_fstr($arg->{dbfields});
	}

	tie my %ops, 'Tie::IxHash';
	foreach my $opstr (split(';', $fstr)) {
		my ($opfunc, $opname) = split(':', $opstr);
		unless(defined $opfunc && defined $opname) {
			die "Cannot find operation name: $opstr";
		}
		if (defined $ops{$opname}) {
			die "Operation column $opname has already been defined!";
		}
		if ($opfunc =~ /^(\w+)\((\S+)\)$/) {
			my ($funcname, $params) = ($1, $2);
			my @params = split(',', $params);
			foreach my $param (@params) {
				unless(grep { $_ eq $param } values %$dbfields) {
					die "Cannot find field $param in dbfile!";
				}
			}
			$ops{$opname} = _create_operator($funcname, \@params, $arg);
		}
		else {
			die "Incorrect specification of operation function: $opfunc"
		}
	}
	return \%ops;
}

# Helper function
sub _create_operator {
	my ($funcname, $params, $arg) = @_;
	my $operator;
	if ($funcname eq 'sum' || $funcname eq 'max' || $funcname eq 'min' ||
		$funcname eq 'mean' || $funcname eq 'median' || $funcname eq 'stdev' || $funcname eq 'mad') {
		unless (@$params == 1) {
			die "Function $funcname only take one parameter";
		}
		# For numerical functions, non-number will be ignored automatically based on regex
		if ($funcname eq 'sum') {
			$operator = sub {
				return sum(grep { $_ =~ /^$RE{num}{real}$/ } map { $_->{$params->[0]} } @_);
			};
		} 
		elsif ($funcname eq 'min') {
			$operator = sub {
				return min(grep { $_ =~ /^$RE{num}{real}$/ } map { $_->{$params->[0]} } @_);
			};
		}
		elsif ($funcname eq 'max' ) {
			$operator = sub {
				return max(grep { $_ =~ /^$RE{num}{real}$/ } map { $_->{$params->[0]} } @_);
			};
		}
		elsif ($funcname eq 'mean' ) {
			$operator = sub {
				return mean(grep { $_ =~ /^$RE{num}{real}$/ } map { $_->{$params->[0]} } @_);
			};
		}
		elsif ($funcname eq 'median' ) {
			$operator = sub {
				return median(grep { $_ =~ /^$RE{num}{real}$/ } map { $_->{$params->[0]} } @_);
			};
		}
		elsif ($funcname eq 'stdev') {
			$operator = sub {
				return std(grep { $_ =~ /^$RE{num}{real}$/ } map { $_->{$params->[0]} } @_);
			};
		}
		elsif ($funcname eq 'mad') {
			$operator = sub {
				return mad(grep { $_ =~ /^$RE{num}{real}$/ } map { $_->{$params->[0]} } @_);
			};
		}
	}
	elsif ($funcname eq 'weightedsum' || $funcname eq 'weightedmean') {
		unless(@$params == 2) {
			die "Function $funcname can only take two parameters";
		}
		# We will use the calculate sum/mean/median of the first parameter using second parameter as weight
		# It requires both fields to be numerical.
		if ($funcname eq 'weightedsum') {
			$operator = sub {
				my @valwts = grep { $_->[0] =~ /^$RE{num}{real}$/ && $_->[1] =~ /^$RE{num}{real}$/ } map { [$_->{$params->[0]},$_->{$params->[1]}] } @_;
				return sum(map { $_->[0]*$_->[1] } @valwts);
			};
		}
		elsif ($funcname eq 'weightedmean') {
			$operator = sub {
				my @valwts = grep { $_->[0] =~ /^$RE{num}{real}$/ && $_->[1] =~ /^$RE{num}{real}$/ } map { [$_->{$params->[0]},$_->{$params->[1]}] } @_;
				my @vals =  map { $_->[0] } @valwts;
				my @wts  =  map { $_->[1] } @valwts;
				return weightedmean(\@vals, \@wts);
			}
		}	
	}
	# The following functions applies to both numeric and character variables
	elsif ($funcname eq 'count') {
		$operator = sub { return scalar(@_) };
	}
	elsif ($funcname eq 'count_distinct') {
		$operator = sub { 
			return scalar(uniq sort map { join($arg->{dfsep}, @{$_}{@$params}) } @_); 
		};
	}
	elsif ($funcname eq 'collapse') {
		$operator = sub {
			return join($arg->{orsep}, map { join($arg->{dfsep}, @{$_}{@$params}) } @_);
		};
	}
	elsif ($funcname eq 'first') {
		$operator = sub {
			my $dat = shift @_;
			return join($arg->{dfsep}, @{$dat}{@$params}); 
		};
	}
	elsif ($funcname eq 'last') {
		$operator = sub {
			my $dat = pop @_;
			return join($arg->{dfsep}, @{$dat}{@$params}); 
		};
	}
	elsif ($funcname eq 'distinct') {
		$operator = sub {
			return join($arg->{orsep}, uniq sort map { join($arg->{dfsep}, @{$_}{@$params}) } @_);
		};
	}
	elsif ($funcname eq 'freqdesc') {
		$operator = sub {
			my %hist;
			foreach my $elem (map { join($arg->{dfsep}, @{$_}{@$params}) } @_) {
				$hist{$elem} ++;
			}
			return join($arg->{orsep}, map { "$_($hist{$_})" } sort { $hist{$b} <=> $hist{$a} || $a cmp $b } keys %hist);
		};
	}
	elsif ($funcname eq 'mode') {
		$operator = sub {
			my %hist;
			foreach my $elem (map { join($arg->{dfsep}, @{$_}{@$params}) } @_) {
				$hist{$elem} ++;
			}
			return (sort { $hist{$b} <=> $hist{$a} } keys %hist)[0];
		};
	}
	elsif ($funcname eq 'freqasc') {
		$operator = sub {
			my %hist;
			foreach my $elem (map { join($arg->{dfsep}, @{$_}{@$params}) } @_) {
				$hist{$elem} ++;
			}
			return join($arg->{orsep}, map { "$_($hist{$_})" } sort { $hist{$a} <=> $hist{$b} || $a cmp $b } keys %hist);
		};
	}
	elsif ($funcname eq 'antimode') {
		$operator = sub {
			my %hist;
			foreach my $elem (map { join($arg->{dfsep}, @{$_}{@$params}) } @_) {
				$hist{$elem} ++;
			}
			return (sort { $hist{$a} <=> $hist{$b} } keys %hist)[0];
		};
	}
	else {
		die "Cannot recognize func name: $funcname";
	}
	return $operator;
}

=head2 anno_fields

Determine the fields in the annotation output from config.

=cut


sub anno_fields {
	my ($conf, $noxtra) = @_;
	my %config;
	if (-f $conf) {
		read_config $conf => %config;
	}
	elsif (ref $conf) {
		%config = %$conf;
	}
	else {
		die "Cannot recognize config data or file: $conf";
	}

	# 
	# Fields from VEP annotations
	#
	my @fields = qw(VarID Chrom Position Ref Alt GeneID GeneEff TransCount TransIDs TransEffs CodonChg AAChg);
	my $veparg;
	if (!defined $config{VEP}{Option} || $config{VEP}{Option} =~ /^\s*$/) {
		# Default options:
		$veparg = "--gencode_basic --hgvs --symbol --transcript_version --biotype --numbers";
	}
	else {
		$veparg = $config{VEP}{Option};
	}
	if ($veparg =~ /\-\-symbol/) {
		insert_before(\@fields, "GeneID", "Symbol");
	}
	if ($veparg =~ /\-\-biotype/) {
		insert_after(\@fields, "TransIDs", "TransBiotypes");
	}
	if ($veparg =~ /\-\-hgvs/) {
		insert_before(\@fields, "CodonChg", "cDNAChg");
	}
	if ($veparg =~ /\-\-canon/) {
		insert_before(\@fields, "TransIDs", "TransCanon");
	}
	if (defined $config{VEP}{Plugin}) {
		unless (ref $config{VEP}{Plugin}) {
			$config{VEP}{Plugin} = [$config{VEP}{Plugin}];
		}
		if (grep { /^Context/ } @{$config{VEP}{Plugin}}) {
			insert_after(\@fields, "Alt", "Context");
		}
		if (grep { /^Ancestral/ } @{$config{VEP}{Plugin}}) {
			insert_after(\@fields, "Alt", "Anc");
		}
		if (grep { /^SpliceRegion/ } @{$config{VEP}{Plugin}}) {
			insert_after(\@fields, "TransEffs", "SpliceReg");
		}
	}
	if (defined $config{VEP}{Gene_Fields}) {
		my $gfields = parse_fstr($config{VEP}{Gene_Fields}, 1);
		insert_before(\@fields, "TransIDs", values %$gfields);
	}
	if (defined $config{VEP}{Trans_Fields}) {
		my $tfields = parse_fstr($config{VEP}{Trans_Fields}, 1);
		push @fields => values %$tfields;
	}
	if (defined $config{VEP}{Custom_Fields}) {
		unless(ref $config{VEP}{Custom_Fields}) {
			$config{VEP}{Custom_Fields} = [$config{VEP}{Custom_Fields}];
		}
		my $vepxtra = join(",", @{$config{VEP}{Custom_Fields}});
		my $vfields = parse_fstr($vepxtra, 1);
		push @fields => values %$vfields;
	}

	#
	# Additional gene and transcript level information
	#
	if (defined $config{Gene}{GXref_Fields}) {
		my @gxrefs;
		unless (ref $config{Gene}{GXref_Fields}) {
			$config{Gene}{GXref_Fields} = [$config{Gene}{GXref_Fields}];
		}
		foreach my $gstr (@{$config{Gene}{GXref_Fields}}) {
			my $gfields = parse_fstr($gstr, 1);
			my @gfields = values %$gfields;
			shift @gfields;
			unless(@gfields > 0) {
				warn "No extra gene info field found: $gstr";
			}
			else {
				push @gxrefs, @gfields;
			}
		}
		insert_after(\@fields, "GeneID", @gxrefs);
	}
	if (defined $config{Transcript}{TXref_Fields}) {
		my @txrefs;
		unless (ref $config{Transcript}{TXref_Fields}) {
			$config{Transcript}{TXref_Fields} = [$config{Transcript}{TXref_Fields}];
		}
		foreach my $tstr (@{$config{Transcript}{TXref_Fields}}) {
			my $tfields = parse_fstr($tstr, 1);
			my @tfields = values %$tfields;
			shift @tfields;
			unless(@tfields > 0) {
				warn "No extra trans info field found: $tstr";
			}
			else {
				push @txrefs, @tfields;
			}
		}
		insert_after(\@fields, "TransIDs", @txrefs);
	}

	#
	# Inhouse D-mis annotations
	#
	if (defined $config{Missense}{DBFile_Fields}) {
		my @dmis;
		unless (ref $config{Missense}{DBFile_Fields}) {
			$config{Missense}{DBFile_Fields} = [$config{Missense}{DBFile_Fields}];
		}
		my @stdfds = qw(Chrom Position Ref Alt GeneID AAChg);
		my %stdfds; @stdfds{@stdfds} = @stdfds;
		foreach my $vstr (@{$config{Missense}{DBFile_Fields}}) {
			my $vfields = parse_fstr($vstr, 1);
			my @vfields = grep { !defined $stdfds{$_} } values %$vfields;
			unless(@vfields > 0) {
				warn "No extra D-mis field found: $vstr";
			}
			else {
				push @dmis, @vfields;
			}
		}
		push @fields, @dmis;
	}

	#
	# ANNOVAR filter/region-based fields
	#
	if (defined $config{ANNOVAR}) {
		my @fstrs;
		if (defined $config{ANNOVAR}{Filter_Fields}) {
			if (ref $config{ANNOVAR}{Filter_Fields} eq 'ARRAY') {
				push @fstrs, @{$config{ANNOVAR}{Filter_Fields}};				
			}
			else {
				push @fstrs, $config{ANNOVAR}{Filter_Fields};
			}
		}
		if (defined $config{ANNOVAR}{Region}) {
			unshift @fstrs, $config{ANNOVAR}{Region};
		}
		unless(@fstrs > 0) {
			warn "No ANNOVAR fields found!";
		}
		my $afields = parse_fstr(join(",", @fstrs), 1);
		push @fields, values %$afields;
	}

	#
	# Extra fields from input and sample level fields
	#
	unless ($noxtra) {
		if (defined $config{Sample} && defined $config{Sample}{IDField}) {
			my $xfields = parse_fstr($config{Sample}{IDField}, 1);
			unshift @fields, values %$xfields;
		}
		if (defined $config{Input}{XtraFields}) {
			my $xfields = parse_fstr($config{Input}{XtraFields}, 1);
			insert_after(\@fields, "Alt", values %$xfields);
		}
		if ($config{Sample}{IDField} && defined $config{Sample}{IXref_Fields}) {
			unless(ref $config{Sample}{IXref_Fields}) {
				$config{Sample}{IXref_Fields} = [$config{Sample}{IXref_Fields}];
			}
			my @sampxrefs;
			foreach my $xstr (@{$config{Sample}{IXref_Fields}}) {
				my $xfields = parse_fstr($xstr, 1);
				my @xfields = values %$xfields;
				shift @xfields;
				unless(@xfields > 0) {
					warn "No extra sample info field found: $xstr";
				}
				else {
					push @sampxrefs, @xfields;
				}
			}
			splice @fields, 1, 0, @sampxrefs;
		}
	}

	my @uqfields = uniq sort @fields;
	unless(scalar(@uqfields) == scalar(@fields)) {
		my %count;
		foreach my $field (@fields) {
			$count{$field} ++;
		}
		my @dupfds;
		foreach my $field (grep { $count{$_} > 1 } sort keys %count) {
			push @dupfds, $field;
		}
		print join(",", @dupfds), "\n";
		die "Duplicated field names found!";
	}

	return @fields;
}


1;
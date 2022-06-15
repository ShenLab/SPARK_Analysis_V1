package Genet::File::BPed;

use strict;
use warnings;
use Carp;
use IO::File;
use Perl6::Slurp;
use List::Util qw|sum|;
use List::MoreUtils qw|all none natatime|;
use Iterator::Simple qw|iterator|;
use POSIX qw|ceil|;
use Utils::Hash qw|merge_opts|;
use Utils::File::Iter qw|slurp_file|;

=head1 NAME

Genet::File::BPed - Read and write PLINK BPed format.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

    use Genet::File::BPed;

    my @fullsamp = map { [ (split)[0..5] ] } slurp "CEU.fam";
    my @subsamp;
    foreach my $ii (0..$#fullsamp) {
        push @subsamp, $fullsamp[$ii] unless $ii % 2;
    }

    my $reader = Genet::File::BPed->new("$ENV{HOME}/Temp/CEU",
        {fiid => 1, gstr => 1, samp => [ map { $_->[0].'-'.$_->[1] } @subsamp ]});

    my $writer = Genet::File::BPed->new("$ENV{HOME}/Temp/CEUSubset",
        {samp => \@subsamp, fiid => 1, write => 1});

    while(my @dat = $reader->read_next) {
        $reader->recalc_major(@dat);
        last if $dat[0]{CHROM} > 1;
        $writer->write_next(@dat);
    }


=head1 SUBROUTINES/METHODS

=head2 $class->new PREFIX [, OPTIONS]

Create an object for read or write genotypes from PLINK BPed file.

When fam/bim/bed files exist, the constructed object is by default used for reading
unless forced to overrite the existing files.

When none of such files exist, the constructed object will be used for writing.
In this case, sample level information must be provided. In its minimal, it can be a
list of sample IDs. The full information should be the same standard 6-column
FAM file used by plink. Each individual can be stored as arrayref or hashref.

=head3 Options

=over 5

=item C<samp> 

For reading, this is the subset of samples for reading, default is to read all.
For writing, this is the full sample list to write into fam file.

=item C<fiid>

The keys to genotype hash is IID by default, use this option to switch to "FID-IID".

=item C<gstr>

The values in the genotype hash will be represented as strings like "minor/major".
Default is undef. Switch off this option to use the number of major alleles.

=item C<write>

Force writing, without check the existing files.

=back


=cut

sub new {
	my ($class, $prefix, $argref) = @_;
	my $arg = merge_opts($argref, samp => undef, fiid => undef,
		 gstr => undef, write => undef);


	my @FAMFIELDS = qw|FID IID DAD MOM SEX PHE|;
	my $MAGICBITS = '011011000001101100000001';

	if (!$arg->{write} && (all { -f "$prefix.$_" } qw|fam bim bed|)) {
		my @samp = slurp_file("$prefix.fam",  {header => 0, fields => [@FAMFIELDS]});

		open my $fbim, "$prefix.bim" or croak "Cannot open bim file : $prefix.bim";
		open my $fbed, "<:raw", "$prefix.bed" or croak "Cannot open bed file : $prefix.bed";

		my $header;	read($fbed, $header, 3);
		croak "Not a SNP-major ordered bed file" unless
			unpack("B24",$header) eq $MAGICBITS;


		my (@all_sampids);
		if ($arg->{fiid}) {
			@all_sampids = map { $_->{FID}.'-'.$_->{IID} } @samp;
		}
		else {
			@all_sampids = map { $_->{IID} } @samp;
		}

		my (@sampind);
		if ($arg->{samp}) {
			croak "Subset of samples should be provided as an arrayref"
				unless ref $arg->{samp} eq 'ARRAY';
			foreach my $iid (@{$arg->{samp}}) {
				my @ind = grep { $iid eq $all_sampids[$_] } 0..$#all_sampids;
				croak "Cannot determine sample $iid in fam file" unless @ind == 1;
				push @sampind, $ind[0];
			}
		}
		else {
			@sampind = 0..$#all_sampids;
		}

		return bless { samp => \@samp, nsamp => scalar(@samp), sampids => \@all_sampids, 
			prefix => $prefix, fhbim => $fbim, fhbed => $fbed, nbyte => ceil(scalar(@samp)/4),
			subind => \@sampind, gstr => $arg->{gstr}, write => 0 }, $class;
	}
	elsif ($arg->{write} || (none { -f "$prefix.$_" } qw|fam bim bed|)) {
		croak "Must provide samp list for writing" 
			unless defined $arg->{samp} && ref $arg->{samp} eq 'ARRAY';
		my $samplst = $arg->{samp};
		my @samp;
		if (ref $samplst->[0] eq 'HASH') {
			foreach my $indiv (@$samplst) {
				croak "Each individual must have : FID IID DAD MOM SEX PHE"
					unless (all { defined $indiv->{$_} } @FAMFIELDS);
				push @samp, $indiv;
			}
		}
		elsif (ref $samplst->[0 eq 'ARRAY']) {
			foreach my $indiv (@$samplst) {
				croak "Each individual should have six columns" unless @$indiv == 6;
				my %dat;
				@dat{@FAMFIELDS} = @$indiv;
				push @samp, \%dat;
			}
		}
		else {
			foreach my $indiv (@$samplst) {
				push @samp, { FID => $indiv, IID => $indiv, DAD => '0', MOM => '0', SEX => '-9', PHE => '-9' };
			}
		}

		my @all_sampids;
		if ($arg->{fiid}) {
			@all_sampids = map { $_->{FID}.'-'.$_->{IID} } @samp;
		}
		else {
			@all_sampids = map { $_->{IID} } @samp;
		}

		# Write the sample file
		my $fout = IO::File->new("$prefix.fam", "w");
		foreach my $indiv (@samp) {
			print $fout join("\t", @{$indiv}{@FAMFIELDS}), "\n";
		}
		$fout->close;

		#open my $fbim, ">$prefix.bim" or croak "Cannot write to bim file : $prefix.bim";
		#open my $fbed, ">:raw", "$prefix.bed" or croak "Cannot write to bed file : $prefix.bed";
		my $fbim = IO::File->new("$prefix.bim", "w") or croak "Cannot write to bim file : $prefix.bim";
		my $fbed = IO::File->new("$prefix.bed", "w") or croak "Cannot write to bed file : $prefix.bed";
		binmode($fbed);

		# Write the magic bits
		print $fbed pack ("B24", "011011000001101100000001");

		return bless { samp => \@samp, nsamp => scalar(@samp), sampids => \@all_sampids, 
			prefix => $prefix, fhbim => $fbim, fhbed => $fbed, nbyte => ceil(scalar(@samp)/4),
			write => 1 }, $class;

	}
	elsif (!$arg->{write} && !(all { -f "$prefix.$_" } qw|fam bim bed|)) {
		croak "Not all fam/bim/bed files exist for reading";
	}
	else {
		croak "Some of fam/bim/bed files for $prefix already exist, please check to avoid clobbering";
	}
}


=head2 Sample Accessors

$self->get_samp: Return sample list, each sample is hashref.

$self->get_sampids : Return sample IDs, either IID or FID-IID depending on the init option.

=cut


foreach my $attr (qw|samp sampids|) {
	no strict 'refs';
	*$attr = sub { my $self = shift; if (wantarray) { @{$self->{$attr}} } else { return $self->{$attr} } };
}



=head2 $self->read_next

Read genotypes for the next variant. It returns \%site and \%geno, in array or arrayref
depending on the context.

The site level information will contain following fields as hash or array.
CHROM, ID, POS, GPOS, MAJOR, MINOR

The genotype will be a hash, indexed by the compound sample name IID or FID-IID.
Values will be the genotype string representation or number of major alleles.

Samples with missing genotypes should exist in the hash, with undef values.

=cut

my %GCODE = ("00" => 0, "01" => 1, "11" => 2, "10" => undef);

sub read_next {
	my ($self) = @_;
	croak "Cannot read a writer object!" if $self->{write};
		
	my $siteline = $self->{fhbim}->getline();
	my $nbyteread = read($self->{fhbed}, my $genobytes, $self->{nbyte},0);

	if ( defined $siteline && $nbyteread != $self->{nbyte} ) {
		croak "bim and bed file are not consistent: $genobytes vs $self->{nbyte}!";
	}
	
	if (defined $siteline && defined $genobytes) {
		my ($chrom, $snp, $genet, $pos, $minor, $major) = split(/\s+/, $siteline);
		my %site = (CHROM => $chrom, POS => $pos, GPOS => $genet, 
				 	ID => $snp, MAJOR => $major, MINOR => $minor);

		my (%geno, %af);
	 	my @geno = map { my $s = unpack("b8",$_);
					  	 (substr($s,0,2), substr($s,2,2),
					  	  substr($s,4,2), substr($s,6,2));
						} split('',$genobytes);
		unless ($self->{gstr}) {
			foreach my $ii (@{$self->{subind}}) {
				$geno{$self->{sampids}[$ii]} = $GCODE{$geno[$ii]};
			}
		}
		else {
			my %GSTR = ("00" => "$minor/$minor", "01" => "$minor/$major", 
						"10" => "./.", "11" => "$major/$major");
			foreach my $ii (@{$self->{subind}}) {
				$geno{$self->{sampids}[$ii]} = $GSTR{$geno[$ii]};
			}
		}
		if (wantarray) {
			return (\%site, \%geno);
		}
		else {
			return [\%site, \%geno];
		}
	}
	else {
		return;
	}
}


=head2 Allele Frequencies.

allele_count SITE, GENO [, SAMP]: calculate allele counts. Returns a hash/hashref.

allele_freq SITE, GENO [, SAMP]: calculate allele frequencies. Returns a hash/hashref.


recalc_major SITE, GENO [, SAMP]: re-calculate the major allele in place


After subsetting the whole sample, allele frequency may change, so major/minor allele
may also change. 

B<NOTE>: by default, all samples were used for allele frequency calculation.
Optionally, users can provide a list of samples (e.g. founders) on whom allele
frequencies are calculated.

=cut

sub allele_count {
	my ($self, $site, $geno, $samp) = @_;

	my @genos;
	if (defined $samp) {
		@genos = map { $geno->{$_} } grep { defined $geno->{$_} } @$samp;
	}
	else {
		@genos = grep { defined $_ } values %$geno;
	}

	my %af = ($site->{MAJOR} => 0, $site->{MINOR} => 0);

	if ($self->{gstr}) {
		foreach my $geno (@genos) {
			foreach my $al (split(q|/|, $geno)) {
				croak "Allele $al is not defined in site" unless defined $af{$al};
				$af{$al} ++;
			}
		}
	}
	else {
		foreach my $geno (@genos) {
			$af{$site->{MAJOR}} += $geno;
			$af{$site->{MINOR}} += 2-$geno;
		}
	}

	if (wantarray) {
		return %af;
	}
	else {
		return \%af;
	}
}

sub allele_freq {
	my ($self, $site, $geno, $samp) = @_;
	my %af = $self->allele_count($site, $geno, $samp);
	my $ntot = sum(values %af);
	if ($ntot > 0) {
		foreach my $al (keys %af) {
			$af{$al} /= $ntot;
		}
	}
	if (wantarray) {
		return %af;
	}
	else {
		return \%af;
	}
}

sub recalc_major {
	my ($self, $site, $geno, $samp) = @_;

	my %af = $self->($site, $geno, $samp);

	if ($af{$site->{MINOR}} > $af{$site->{MAJOR}}) {
		($site->{MAJOR}, $site->{MINOR}) = ($site->{MINOR}, $site->{MAJOR});
		unless ($self->{gstr}) {
			foreach my $iid (keys %$geno) {
				if (defined $geno->{$iid}) {
					$geno->{$iid} = 2 - $geno->{$iid};
				}
			}
		}
	}
	return $self;
}

=head2 Genotyping Missingness

geno_miss : the number of missing genotypes

geno_miss_rate : the missing genotype rate

geno_nomiss : the number of non-missing genotyes

geno_nomiss_rate : the non-missing genotype rate

=cut

sub geno_miss {
	my ($self, $geno) = @_;
	my $miss = grep { !defined $_ } values %$geno;
	return $miss;
}

sub geno_nomiss {
	my ($self, $geno) = @_;
	return @{$self->{subind}} - $self->geno_miss($geno);
}

sub geno_miss_rate {
	my ($self, $geno) = @_;
	return $self->geno_miss($geno)/@{$self->{subind}};
}

sub geno_nomiss_rate {
	my ($self, $geno) = @_;
	return 1-$self->geno_miss_rate($geno);
}


=head2 $self->write_next

Return a bped file writer.

Each call of the write require site and geno as input, which are the same as returned
from the tierator. It will consecutively write to the stored file handler.

The required input is site and genotype level information, the same as returned from
the iterator. Genotypes will be the number of major alleles.


=cut

my %BCODE = ('0' => "00", '1' => "01", '2' => "11", "NA" => "10");
my @BIMFIELDS = qw|CHROM ID GPOS POS MINOR MAJOR|;

sub write_next {
	my ($self, $site, $geno) = @_;
	croak "Cannot write a reader object!" unless $self->{write};

	#*FHBIM = $self->{fhbim};
	#*FHBED = $self->{fhbed};
	#print FHBIM join("\t", @{$site}{@BIMFIELDS}), "\n";
	my ($minor, $major);
	if (ref $site eq 'ARRAY') {
		$self->{fhbim}->print(join("\t", @$site), "\n");
		($minor, $major) = ($site->[4], $site->[5]);
	}
	elsif (ref $site eq 'HASH') {
		print STDERR "Fields required: ", join(", ", @BIMFIELDS), "\n";
		croak "Not all fields can be found in site";
		$self->{fhbim}->print(join("\t", @{$site}{@BIMFIELDS}), "\n");
		($minor, $major) = @{$site}{qw|MINOR MAJOR|};
	}
	else {
		croak "Incorrect site info";
	}

	my $it = natatime 4, @{$self->{sampids}};
	while(my @samps = $it->()) {
		my $gbyte = join "", map { 
			if (!defined $_) {
				$BCODE{'NA'}
			}
			elsif ($_ eq '0' || $_ eq "$minor/$minor") {
				$BCODE{'0'}
			}
			elsif ($_ eq '1' || $_ eq "$major/$minor" || $_ eq "$minor/$major") {
				$BCODE{'1'}
			}
			elsif ($_ eq '2' || $_ eq "$major/$major") {
				$BCODE{'2'}
			}
			else {
				$BCODE{'NA'}
			}
		 }  @{$geno}{@samps};
		 #print FHBED pack("b8", $gbyte)
		 $self->{fhbed}->print(pack("b8", $gbyte));
	}
	return 1;
}


=head1 AUTHOR

Xueya Zhou, C<< <xueyazhou at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-genet at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Genet>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Genet::File::BPed


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Genet>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Genet>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Genet>

=item * Search CPAN

L<http://search.cpan.org/dist/Genet/>

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

1; # End of Genet::File::BPed

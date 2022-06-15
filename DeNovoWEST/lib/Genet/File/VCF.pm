package Genet::File::VCF;

use strict;
use warnings;
use Carp;
use IO::Detect;
use List::Util qw|sum min max|;
use List::MoreUtils qw|all uniq|;
use Iterator::Simple qw|iterator|;
use IO::Detect;
use File::Which;
use Text::ParseWords qw|parse_line|;
use Tie::IxHash;
use Utils::Hash qw|merge_opts array2hash|;
use Utils::File qw|open_file|;
use Genome::UCSC qw|hg_par|;
use Data::Dumper;

use base qw|Exporter|;

our @EXPORT_OK = qw|parse_vcf_header slurp_vcf_header check_vcf_chrpref allele_count|;
our %EXPORT_TAGS = ( 'all' => \@EXPORT_OK );


=head1 NAME

Genet::File::VCF - Parsing VCF file.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

    use Genet::File::VCF;
    use Perl6::Slurp;

    my %pedline = map { my @row = split; $row[1] => [@row] } slurp "CDH_MGH.ped";

    my $vcf = Genet::File::VCF->new("$ENV{HOME}/temp/CDH_MGH_VQSR.vcf.gz");

    my @samps = grep { /^\d+$/ } $vcf->get_sampids;

    my $vit = $vcf->iter({strict => 1, samp => \@samps, region => "20:43000000-63000000"});

    my $writer = Genet::File::BPed->new("$ENV{HOME}/Temp/DH_MGH_VQSR",
                  {samp => [map { $pedline{$_} } @samps], write => 1});

    while(my ($site, $geno) = $vit->()) {

        my $varid = "$site->{CHROM}_$site->{POS}_$site->{REF}_$site->{ALT}[0]";
        $site->{ID} = $varid;

        next unless scalar(@{$site->{ALT}}) == 1 && $site->{FILTER}[0] eq 'PASS';

        my %geno; 
        while(my ($iid, $gdat) = each %$geno) {
            if (defined $gdat->{GT}) {
                my $gct = sum( split(q|/|, $gdat->{GT}) );
                croak "Incorrect number of ref allele count" unless $gct >= 0 && $gct <= 2;
                $geno{$iid} = 2-$gct;
        }
        else {
             $geno{$iid} = undef;
        }
    }

    $writer->write_next([$site->{CHROM},$site->{ID},0,$site->{POS},$site->{ALT}[0],$site->{REF}],\%geno);
}

=head1 DESCRIPTION

The main functionality of this module is to read and write (limited) the VCF file.
The performance is slow because it unpack all the information stored in aline.

See also L<Vcf> distributed with VCFtools for more functionailities.


=head1 SUBROUTINES/METHODS

=head2 new FILE [, OPTIONS]

Create a new VCF file reader object.

FILE can be a file or file handle. In the latter case, it only supports sequential reading
from the beginning.

=cut

sub new {
	my ($class, $file) = @_;

	my $ftype = is_filehandle($file);
	my ($header, $info, $sampids) = parse_vcf_header($file);
	#print Dumper $info;

	return bless { file => $file, ftype => $ftype, header => $header, 
		sampids => $sampids, %$info }, $class;

}

=head2 Helper Functions

=head3 parse_vcf_header FILE

  Return header, sample names, and field information from VCF file or file handle.
  The field information will extract from FILTER/FORMAT/INFO.

=cut

my @STDFIELDS = qw|CHROM POS ID REF ALT QUAL FILTER INFO FORMAT|;
$STDFIELDS[0] = '#'.$STDFIELDS[0];

sub parse_vcf_header {
	my ($file) = @_;
	my $fin;
	if (is_filehandle($file)) {
		$fin = $file;
	}
	else {
		 $fin = open_file($file);
	}

	my ($header, %info, @sampids);
	#%info = (FILTER => {}, INFO => {}, FORMAT => {});
	foreach my $Class (qw|FILTER INFO FORMAT|) {
		tie my %classdat, 'Tie::IxHash';
		$info{$Class} = \%classdat;
	}
	my $line;
	while ($line = <$fin>) {
		croak "No header found in the input: $line" unless $line =~ /^#/;
		last if $line =~ /^#CHROM/;

		$header .= $line;
		if ($line =~ /^##([^=<]+)=\<(.+)\>$/) {
			my ($label, $content) = ($1, $2);
			next unless defined $info{$label};
			my %dat = map { $_ =~ /^([^=]+)=(.+)$/; ($1 => $2) } parse_line(',', 1, $content);
			# print $label, "\t", join(" ", keys %dat), "\n";
			$info{$label}{$dat{ID}} = \%dat;
  		}
	}

	chomp($line);
	my @fields = split(/\t/, $line);
	for(my $ii = 0; $ii < min(scalar(@fields), 8); $ii ++) {
		croak "Field name unmatch: $fields[$ii] <> $STDFIELDS[$ii]" 
			unless $fields[$ii] eq $STDFIELDS[$ii];
	}
	if (@fields > 9) {
		for(my $ii = 9; $ii < @fields; $ii ++) {
			push @sampids, $fields[$ii];
		}
	}
	my @siduq = uniq sort @sampids;
	croak "Duplicated sample names exist in VCF header" unless @sampids == @siduq;

	if (wantarray) {
		return ($header, \%info, \@sampids);
	}
	else {
		return [$header, \%info, \@sampids];
	}
}

=head3 slurp_vcf_header FILE

  Return header, and a list of sample names from VCF file or file handle. 
  It did not parse header information fields.

=cut

sub slurp_vcf_header {
	my ($file) = @_;
	my $fin = open_file($file);

	my ($header, @sampids);
	my $line;
	while ($line = <$fin>) {
		croak "No header found in the input" unless $line =~ /^#/;
		last if $line =~ /^#CHROM/;
		$header .= $line;
	}

	chomp($line);
	my @fields = split(/\t/, $line);
	my $fupper;
	if (@fields > 9) {
		$fupper = 9;
	}
	else {
		$fupper = @fields;
	}
	for(my $ii = 0; $ii < $fupper; $ii ++) {
		croak "Field name unmatch: $fields[$ii] <> $STDFIELDS[$ii]" 
			unless $fields[$ii] eq $STDFIELDS[$ii];
	}
	for(my $ii = 9; $ii < @fields; $ii ++) {
		push @sampids, $fields[$ii];
	}
	my @siduq = uniq sort @sampids;
	croak "Duplicated sample names exist in VCF header" unless @sampids == @siduq;

	if (wantarray) {
		return ($header, \@sampids);
	}
	else {
		return [$header, \@sampids];
	}
}

# Check if VCF use chr prefix in chromosome names
sub check_vcf_chrpref {
	my ($vcfpath) = @_;
	my ($vcfchr, $fvc);
	if ($vcfpath =~ /\.vcf\.gz$/) {
		$fvc = IO::Uncompress::Gunzip->new($vcfpath, MultiStream => 1)
			or die "Cannot open zipped VCF: $vcfpath";
	}
	else {
		open $fvc, $vcfpath or die "Cannot open vcf: $vcfpath:";
	}
	my $chr;
	while(<$fvc>) {
		unless (/^#/) {
			$chr = (split)[0];
			last;
		}
	}
	die "Cannot determine vcf chr name" unless defined $chr;
	if ($chr =~ /^chr/) {
		$vcfchr = 1;
	}
	else {
		$vcfchr = 0;
	}
	return $vcfchr;
}

=head2 list_vqsr_filter {TYPE => THRESHOLD}

   Given threshold and variant type, list all GATK's VQSR filters that PASS (<=) the threshold.
   When VQSR filter flag is absent, it will output only PASS.
   Note: the returned list will be filters that will KEEP the variant.

=cut

sub list_vqsr_filter {
	my ($self, $argref) = @_;
	my $arg = merge_opts($argref, SNP => 99.9, INDEL => 99.0);
	# parse all VQSR filters
	my %filters;
	foreach my $filt (keys %{$self->{FILTER}}) {
		if ($filt =~ /^VQSRTranche([A-Z]+)([0-9\.]+)to([0-9\.]+)$/) {
			my ($vartype, $lower, $upper) = ($1, $2, $3);
			$filters{$vartype}{$filt} = [$lower, $upper];
		}
	}
	my @filtout = qw|PASS|;
	# Now extract filters
	foreach my $vartype (qw|SNP INDEL|) {
		next unless defined $arg->{uc($vartype)};
		my $varfilt = $filters{$vartype};
		unless (defined $varfilt) {
			carp "Cannot find VQSR filter for $vartype";
		}
		else {
			push @filtout, grep { $varfilt->{$_}[1] <= $arg->{uc($vartype)} } keys %{$varfilt};
		}
	}
	if (wantarray) {
		return @filtout;
	}
	else {
		return [@filtout];
	}
}


=head2 Accessors

get_sampids : get a list of sample IDs

get_info : get all INFO fields

get_filter : get all FILTER labels

get_format : get all FORMAT tags

=cut

sub get_sampids {
	my $self = shift @_;
	if (wantarray) {
		return @{$self->{sampids}};
	}
	else {
		return $self->{sampids};
	}
}

foreach my $attr (qw|info filter format|) {
	no strict 'refs';
	*$attr = sub {
		my $self = shift;
		if (wantarray) { 
			return %{$self->{uc($attr)}} 
		}
		else { 
			return $self->{uc($attr)} 
		} 
	};
}

=head2 get_vcf_fields

Get standard VCF fields plus sample IDs.

It is possible to provide a subset of samples (including empty list).

=cut

sub get_vcf_fields {
	my ($self, $sampids) = shift;
	my @vcfsamps = $self->get_sampids();
	my %allsamps = array2hash(@vcfsamps);
	if (defined $sampids) {
		croak "Must provide list of sample IDs with arrayref" unless ref $sampids eq 'ARRAY';
		# verify that sampleids are subset of existing sammples.
		croak "Provided samples must be a subset of existing samples" unless (all { defined $allsamps{$_} } @$sampids);
		return (@STDFIELDS, @$sampids);
	}
	else {
		return (@STDFIELDS, @vcfsamps);
	}
}

=head2 $self->iter [OPTIONS]

Return an iterator for reading VCF file one line at a time.

The return value of calling the iterator is SITE, GENO and INFO.
If file handle instead of file name is provided to construct object,
then iterator will start from the current position of stored file handle.

The following fields will be included in site-level information:
C<CHROM, POS, ID, REF, ALT, QUAL>, and C<FILTER>, where C<ID, ALT> and C<FILTER> will 
be stored as arrayref to accomondate possible multiple values.

For genotypes, the returned hash are index by sample IDs. Each value is a hashref
with field names from FORMAT definition.


=head3 Options

=over 5

=item C<samp>

Fetch genotypes for a subset of samples. This only affect the returned genotypes.
Values in the INFO will not be recalculated, unless C<recalc> is turned on.

=item C<recalc>

To recalculate AC/AN in the INFO field. If AC/AN is not defined in the original VCF, 
then new field will be added. 

NOTE: In case the number of observed alleles is less than that in the original VCF,
it will NOT change ALT field, because many other fields depend on the number of ALT
alleles, e.g. AC. When this happen, some alleles will have count 0 in AC field.

=item C<strict>

By default, we do not check against field/info/filter definitions at the header.
It will further slow down the process, but can be enabled by turn on C<strict> option.

=item C<region>

User can also provide a region specification ("chr:st-ed"). Only valid for bgzip'ed VCF file.

=item C<order>

Keep original order of INFO fields and FORMAT tags. This may be useful when write
the information back to VCF line. 

=item C<keepline>

Also output original line.

=cut

sub iter {	
	my ($self, $argref)  = @_;
	my $arg = merge_opts($argref, samp => undef, recalc => undef,
		strict => 1, region => undef, order => undef, keepline => undef);

	# index for the provided list of subsamps
	my (@sampind);
	my $nsamp = @{$self->{sampids}};
	if (defined $arg->{samp}) {
		foreach my $iid (@{$arg->{samp}}) {
			my @ind = grep { $iid eq $self->{sampids}[$_] } 0..$nsamp-1;			
			croak "$iid is not found or found more than once" unless @ind == 1;
			push @sampind, $ind[0];
		}
	}
	else {
		@sampind = 0..$nsamp-1 if $nsamp > 0;
	}

	my $fin;
	unless ($arg->{region}) {
		if ($self->{ftype}) {
			$fin = $self->{file};
		}
		else {
			$fin = open_file($self->{file});
		}
	}
	else {
		croak "Region-based query only support bgzip'ed VCF" unless $self->{file} =~ /\.vcf\.gz$/;
		croak "Cannot find tabix index file" unless -f "$self->{file}.tbi";
		my $tabix = which("tabix");
		croak "Cannot find tabix executable" unless defined $tabix;
		open $fin, "$tabix -p vcf $self->{file} $arg->{region} |" or croak "Cannot open pipe";
	}

	if ($arg->{recalc}) {
		if ($arg->{strict}) {
			croak "Cannot find AC/AN in the orginal file" unless defined $self->{INFO}{AC} && defined $self->{INFO}{AN};
		}
		else {
			if (!defined $self->{INFO}{AC}) {
				$self->{INFO}{AC} = {ID => "AC", Number => "A", Type => "Integer", 
					Description => "Allele count in genotypes, for each ALT allele, in the same order as listed"
				};
			}
			if (!defined $self->{INFO}{AN}) {
				$self->{INFO}{AN} = {ID => "AC", Number => "A", Type => "Integer", 
					Description => " Total number of alleles in called genotypes"
				};
			}
		}
	}

	return iterator {
		while(my $line = <$fin>) {
			next if $line =~ /^#/;
			chomp($line);
			my @values = split(/\t/, $line);
			my %site;
			$site{CHROM} = $values[0];
			$site{POS}   = $values[1];
			# VCF spec allow multiple IDs
			$site{ID}    = $values[2] eq '.' ? undef : [ split(';',$values[2]) ];
			$site{REF}   = $values[3];
			# If ALT=., it is a mono-morphic site
			$site{ALT}   = $values[4] eq '.' ? undef : [ map { $_ } split(',',$values[4]) ];
			$site{QUAL}  = $values[5] eq '.' ? undef : $values[5];
			$site{FILTER}= $values[6] eq '.' ? undef : [ split(';',$values[6]) ];
			if (defined $site{FILTER}) {
				foreach my $filt (grep { $_ ne 'PASS' } @{$site{FILTER}}) {
					croak "FILETER flag $filt is not defined" unless defined $self->{FILTER}{$filt};
				}
			}
			my %info;
			tie %info, "Tie::IxHash" if $arg->{order};
			if ($values[7] ne '.') {
				foreach my $kvpair (split(';', $values[7])) {
					my ($key, $val) = split('=', $kvpair);
					if ($arg->{strict}) {
						croak "INFO field $key is not defined" unless defined $self->{INFO}{$key};
					}
					if (!defined $val) {
						$info{$key} = 1;
					}
					else {
						if ($val =~ /,/) {
							my @vals = split(',', $val);
							$info{$key} = \@vals;
						}
						else {
							$info{$key} = $val;
						}
					}
				}
			}
			$site{INFO} = \%info;


			my %geno; # genotype data required
			if (@sampind) {
				tie %geno, "Tie::IxHash" if $arg->{order};
				my @format = split(':', $values[8]);
				if ($arg->{strict}) {
					foreach my $fmt (@format) {
						croak "FORMAT tag $fmt is not defined" unless defined $self->{FORMAT}{$fmt};
					}
				}
				foreach my $jj (@sampind) {
					my %data;
					tie %data, "Tie::IxHash" if $arg->{order};
					my @vals = split(':', $values[$jj+9]);
					if ($arg->{strict}) {
						unless (@format == @vals) {
							print join("\t", @format), "\n", join("\t", @vals), "\n";
							croak "Number of fields in format does not match values";
						}
					} 
					for (my $ii = 0; $ii < @vals; $ii ++) {
						if ($vals[$ii] eq '.') {
							$data{$format[$ii]} = undef;
						}
						elsif ($vals[$ii] =~ /,/) {
							$data{$format[$ii]} = [ split(',', $vals[$ii]) ];
						}
						else {
							$data{$format[$ii]} = $vals[$ii];
						}
					}
					$geno{$self->{sampids}[$jj]} = \%data;
				}
			}
			if ($arg->{recalc}) {
				allele_count(\%site, \%geno);
			}
			if (wantarray) {
				if ($arg->{keepline}) {
					return ( \%site, \%geno, $line );
				}
				else {
					return ( \%site, \%geno );
				}
			}
			else {
				if ($arg->{keepline}) {
					return [ \%site, \%geno, $line ];
				}
				else {
					return [ \%site, \%geno ];
				}	
			}
		}
		return;
	};
}


=head2 Functions on Genotypes

allele_count SITE, GENO [, SAMP]: re-calculate AC/AN, and store it back to the INFO fields.

NOTE: the optional SAMP list can be provided when calculating allele count.

=cut

sub allele_count {
	my ($site, $geno, $samp) = @_;

	return 0 unless defined $site->{ALT};

	my @genos;
	if (defined $samp) {
		@genos = map { split(/\||\//, $geno->{$_}{GT}) } grep { defined $geno->{$_}{GT} } @$samp;
	}
	else {
		@genos = map { split(/\||\//, $_->{GT}) } values %$geno;
	}

	my @ac = (0) x (@{$site->{ALT}});
	my $nref = 0;
	foreach my $g (@genos) {
		next if $g eq '.';
		croak "Genotype code $g > number of alleles" if $g > @{$site->{ALT}};
		if ($g > 0) {
			$ac[$g-1] ++;
		}
		elsif ($g == 0) {
			$nref ++;
		}
		else {
			croak "Cannot recognize genotype code $g";
		}
	}

	# Update site level information
	$site->{INFO}{AC} = [@ac];
	$site->{INFO}{AN} = $nref + sum(@ac);

	return 1;
}

=head2 self->create_vcf_line SITE, GENO [, OPTIONS]

Return a VCF line for writing.

Optionally specify the format fields, info fields, and subset of samples.

=cut

sub create_vcf_line {
	my ($self, $site, $geno, $argref) = @_;
	my $args = merge_opts($argref, info => undef, format => undef, subset => undef);

	my @outsamps;
	if (defined $args->{subset}) {
		@outsamps = @{$args->{subset}};
		my %allsamps = array2hash($self->get_sampids());
		croak "Provided samples do not exist" unless (all { defined $allsamps{$_} } @outsamps);
	}
	else {
		@outsamps = $self->get_sampids();
	}

	my @row;
	foreach my $key (qw|CHROM POS ID REF ALT QUAL FILTER INFO|) {
		if (defined $site->{$key}) {
			my $info = $site->{$key};
			if ($key eq 'ID' || $key eq 'FILTER') {
				push @row, join(';', @$info);
			}
			elsif ($key eq 'ALT') {
				push @row, join(',', @$info);
			}
			elsif ($key eq 'INFO') {
				my @ifields;
				if (defined $args->{info}) {
					croak "Provided INFO fields do not all exit" unless (all { defined $self->{INFO}{$_} } @{$args->{info}});
					@ifields = @{$args->{info}};
				}
				else {
					@ifields = keys %$info;
				}
				push @row, join(";", map { ref $info->{$_} eq 'ARRAY' ?
						"$_=".join(",", @{$info->{$_}}) : "$_=$info->{$_}" } @ifields);	
			}
			else {
				push @row, $info;
			}
		} 
		else {
			push @row, ".";
		}
	 }

	 if (defined $geno) {
	 	my @format;
	 	if (defined $args->{format}) {
	 		croak "Provided FORMAT fields do not all exit" unless (all { defined $self->{FORMAT}{$_} } @{$args->{format}});
	 		@format = @{$args->{format}};
	 	}
	 	else {
	 		my $gt = (values %$geno)[0];
	 		@format = keys %$gt;
	 	}
	 	push @row, join(":", @format);
	 	foreach my $iid (@outsamps) {
	 		my @gt;
	 		foreach my $tag (@format) {
	 			if (defined $geno->{$iid}{$tag}) {
	 				my $val = $geno->{$iid}{$tag};
	 				if (ref $val eq 'ARRAY') {
	 					push @gt, join(',', @$val);
	 				}
	 				else {
	 					push @gt, $val;
	 				}
	 			}
	 			else {
	 				if ($tag eq 'GT') {
	 					push @gt, './.';
	 				}
	 				else {
	 					push @gt, '.'
	 				}
	 			}
	 		}
	 		push @row, join(":", @gt);
	 	}
	 }
	 return join("\t", @row);
} 




=head1 AUTHOR

Xueya Zhou, C<< <xueyazhou at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-genet at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Genet>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Genet::File::VCF


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

1; # End of Genet::File::VCF

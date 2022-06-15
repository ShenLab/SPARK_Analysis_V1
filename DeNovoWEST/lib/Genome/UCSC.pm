package Genome::UCSC;

use strict;
use warnings;
use Carp;
use Data::Dumper;
use List::MoreUtils qw(all minmax);
use Utils::Hash qw|merge_opts|;
use Utils::File::Iter qw|iter_file|;
use Genome::Ranges;
use Genome::UCSC::DB;
use Genome::UCSC::BinKeeper;

use base qw|Exporter|;
our @EXPORT_OK = qw|hg_par hg_chrom hg_chr is_hgchr cyto_band band_range
					%PAR iter_genetab slurp_genetab|;
our %EXPORT_TAGS = ( 'all' => \@EXPORT_OK );


=head1 NAME

Genome::UCSC - Utility functions to UCSC database!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

  use Genome::UCSC;

  my $queryer = cyto_band("hg19");
  my $band = $queryer($range);

  if (hg_par(X, 20000, "hg19") == 0) {
      print "Not in PAR region\n";
  }


=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 hg_par HGBUILD, CHROM, POS

Test if a genomic position on sex chromosome is within pseudo-autosomal region.
Return 0 if not in PAR, 1 if in PAR1, and 2 if in PAR2.

=cut

our %PAR = ( hg17 =>
			{ X1 => 2692881, X2 => 154494747, Y2 => 57372174 },
			 b35 => 
			{ X1 => 2692881, X2 => 154494747, Y2 => 57372174 }, 
			hg18 => 
			{ X1 => 2709520, X2 => 154584238, Y2 => 57443438 },
			b36 =>
			{ X1 => 2709520, X2 => 154584238, Y2 => 57443438 },
			hg19 => 
			{ X1 => 2699520, X2 => 154931044, Y2 => 59034050 },
			b37 =>
			{ X1 => 2699520, X2 => 154931044, Y2 => 59034050 },
			hg38 => 
			{ X1 => 2781479, X2 => 155701382, Y2 => 56887902 },
			b38 =>
			{ X1 => 2781479, X2 => 155701382, Y2 => 56887902 },
		); 

sub hg_par {
	my ($chrom, $pos, $hg) = @_;
	croak "Must provide human genome assembly version" unless defined $hg;
	$hg = lc($hg);
	croak "Genome assembly $hg is not supported" unless defined $PAR{$hg};
	if ($chrom eq "chrX" || $chrom eq 'X') {
		return 1 if $pos < $PAR{$hg}{X1};
		return 2 if $pos >= $PAR{$hg}{X2};
		return 0;
	}
	elsif ($chrom eq "chrY" || $chrom eq 'Y') {
		return 1 if $pos < $PAR{$hg}{X1};
		return 2 if $pos >= $PAR{$hg}{Y2};
		return 0;
	}
	else {
		return 0;
	}
}

=head2 hg_chrom ALIAS

Convert human chromsome names to UCSC format.

Supporting convertion from PLINK's internal representation and NCBI's chromsome code

Currently only consider chromsomes of primary hg assembly.

=cut

my %hg_chromnames = (XY => 'chrX', 23 => 'chrX', 24 => 'chrY', 25 => 'chrX', 26 => 'chrM', MT => 'chrM');
for my $chr (1..22, 'X', 'Y') {
	$hg_chromnames{"$chr"} = "chr".$chr;
}

sub hg_chrom {
	my ($alias) = @_;
	return $hg_chromnames{$alias} // $alias;
}


=head2 is_hgchr CHROM

Test if given name is canonical human chromosome name.

=cut

my %hg_chroms;
while(my ($chr, $name) = each %hg_chromnames) {
	$hg_chroms{$chr} = 1;
	$hg_chroms{$name} = 1;
}

sub is_hgchr {
	my ($chrom) = @_;
	return $hg_chroms{$chrom};
}

=head2 hg_chr CHROM, POS, HG

Convert human chromosome names to PLINK internal code

It will deal PAR regions automatically, and only consider chromsomes of primary hg assembly.

=cut


sub hg_chr {
	my ($chrom, $pos, $hg) = @_;
	$chrom =~ s/^chr//;
	if ($chrom =~ /^(\d+)$/) {
		if ($1 >= 1 && $1 <= 26) {
			return $1;
		}
		else {
			return 0;
		}
	}
	elsif ($chrom eq 'M' || $chrom eq 'MT') {
		return 26;
	}
	elsif ($chrom eq 'X' || $chrom eq 'Y') {
		if (defined $pos && defined $hg && hg_par($chrom, $pos, $hg) > 0) {
			return 25;
		}
		elsif ($chrom eq 'X') {
			return 23;
		}
		else {
			return 24;
		}
	}
	elsif ($chrom eq 'XY') {
		return 25;
	}
	else {
		return 0;
	}
}


=head2 cyto_band DB

Return a cytoband queryer, which is a subroutine that produce cytoband
given genomic interval.

=cut

sub cyto_band {
	my ($dbname) = @_;
	my $dbo = Genome::UCSC::DB->new($dbname);
	my $it = $dbo->iter("SELECT chrom, chromStart+1, chromEnd, name FROM cytoBand",
		{array => 1});
	my $bk = Genome::UCSC::BinKeeper->new();
	while(my $dat = $it->()) {
		$bk->add(@$dat);
	}
	my $queryer = sub  {
		my ($chrom, $chromStart, $chromEnd, $argref) = @_;
		my $arg = merge_opts($argref, qCover => 0);

		(my $chr = $chrom) =~s/^chr//;
		$chromEnd = $chromStart unless defined $chromEnd;
		if ($chromEnd < $chromStart) {
			croak "Incorrect interval: chromEnd < chromStart";
		}
		my @bands = map { $_->[2] } $bk->find_range($chrom, $chromStart, $chromEnd, {qCover => $arg->{qCover}});
		if ( @bands == 0) {
			return "NA";
		}
		elsif ( @bands == 1 ) {
			return $dbname =~ /^dm/ ? $bands[0] : $chr.$bands[0];
		}
		else {
			if ($dbname =~ /^dm/) {
				return $bands[0]."-".$bands[-1];
			}
			else {
				if ( all {$_ =~ /^p/} @bands ) {
					return $chr.$bands[-1]."-".$bands[0];
				}
				elsif ( all {$_ =~ /^q/} @bands ) {
					return $chr.$bands[0]."-".$bands[-1];
				}
				else {
					warn "Cytoband range crosses the centromere";
					return $chr.$bands[0]."-".$bands[-1];
				}
			}
		}
	};
}

=head2 band_range DB

Return a genomic range queryer, which is a subroutine that produce genomic positions
given cytoband.

=cut

sub band_range {
	my ($dbname) = @_;
	if ($dbname =~ /^dm/) {
		warn "No support for Dmel";
		return;
	}
	my $dbo = Genome::UCSC::DB->new($dbname);
	my $it = $dbo->iter("SELECT chrom, chromStart+1, chromEnd, name FROM cytoBand",
		{array => 1});
	my $rg = Genome::Ranges->new();
	while(my $dat = $it->()) {
		$rg->add(@$dat);
	}
	my $queryer = sub  {
		my ($bandrange) = @_;
		my ($start_band, $end_band) = split('\-', $bandrange);
  		my $chrom;
  		my @ranges;
  		if ($start_band =~ /^([1-9XY][0-9]*)pter/) {
  			$chrom = $1;
  			push @ranges, 1;
  		}
  		elsif ($start_band =~ /^([1-9XY][0-9]*)cen/) {
  			$chrom = $1;
  		}
  		elsif ($start_band =~ /^([1-9XY][0-9]*)([pq][0-9.]*)$/) {
  			$chrom = $1;
  			my $name = $2;
			push @ranges, map { ($_->[0], $_->[1]) }
	  			grep { $_->[2] =~ /^$name/ } @{$rg->{"chr$chrom"}};
  		}
  		if (defined $end_band && $end_band =~ /^([pq][0-9.]*)$/) {
			my $name = $1;
			push @ranges, map { ($_->[0], $_->[1]) }
	  			grep { $_->[2] =~ /^$name/ } @{$rg->{"chr$chrom"}};
  		}
  		elsif (defined $end_band && $end_band eq 'qter') {
  			# find the end of chromosome
  			my $nn = @{$rg->{"chr$chrom"}};
  			push @ranges, $rg->{"chr$chrom"}[$nn-1][1];
  		}

  		if (@ranges >= 1) {
			my ($start, $end) = minmax(@ranges);
			if (wantarray) {
				return ($chrom, $start, $end);
			}
			else {
				return [$chrom, $start, $end];
			}
 		}
 		else {
 			warn "Error band specification: $bandrange";
 			return undef;
 		}
	}
}

=head2 iter_gentab and slurp_gentab

Iteartor access to UCSC gene table, similar to C<iter_file> and C<slurp_file>

=cut

sub iter_genetab {
	my ($prefix, $argref) = @_;
	croak "Cannot find $prefix.sql" unless -f "$prefix.sql";
	croak "Cannot find $prefix.txt.gz" unless -f "$prefix.txt.gz";
	my @fields;
	my $fin = IO::File->new("$prefix.sql");
	# Parse SQL file to fetch field names.
	while(<$fin>) {
		last if /CREATE TABLE/;
	}
	while(<$fin>) {
		last if /KEY/;
		my $field = (split)[0];
		$field =~ s/`//g;
		push @fields, $field;
	}
	my %arg;
	if (defined $argref) {
		%arg = %$argref;
	}
	my $it = iter_file("$prefix.txt.gz", {header => 0, fields => \@fields, %arg});
	if (wantarray) {
		return ($it, \@fields);
	}
	else {
		return $it;
	}
}

sub slurp_genetab {
	my ($file, $argref) = @_;
	my $iter = iter_genetab($file, $argref);
	my @out;
	while(my $dat = $iter->()) {
		push @out, $dat;
	}
	if (wantarray) {
		return @out;
	}
	else {
		return \@out;
	}
}

=head1 AUTHOR

Xueya Zhou, C<< <xueyazhou at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-genome at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Genome>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Genome::UCSC


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


=head1 ACKNOWLEDGEMENTS


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

1; # End of Genome::UCSC

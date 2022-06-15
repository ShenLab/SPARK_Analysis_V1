package Genome::UCSC::GeneTable;

use strict;
use warnings;
use Carp;
use DBI;
use Array::IntSpan;
use List::Util qw|min max sum|;
use List::MoreUtils qw|all pairwise|;
use Genome::UCSC::DB;
use Utils::File::Iter;

use base qw|Exporter|;
our @EXPORT_OK = qw|iter_geneTab iter_genePred
					gdat_list_exons gdat_mrna_length gdat_list_introns
					gdat_is_coding gdat_list_coding_exons gdat_cds_length
					gdat_list_5prime_utrs gdat_list_3prime_utrs 
					gdat_list_utrs gdat_to_bed12_line |;
our %EXPORT_TAGS = ( 'all' => \@EXPORT_OK, 
					 'iter' => [grep { /^iter/ } @EXPORT_OK], 
					 'gdat' => [grep { /^gdat/ } @EXPORT_OK] );

=head1 NAME

Genome::UCSC::GeneTable - Access to UCSC's gene tables!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

	use Genome::Ranges;
	use Genome::Ranges::IntSet;
	use Genome::UCSC::Genes qw|iter_genePred gdat_list_coding_exons|;

	my $git = iter_genePred("$ENV{HOME}/temp/ensGene.txt.gz");

	my $rng = Genome::Ranges->new();
	while(my $gdat = $git->()) {
		foreach my $cds (gdat_list_coding_exons($gdat)) {
			$rng->add($gdat->{chrom}, $cds->[0]-2, $cds->[1]+2);
		}
	}
	my $exome = Genome::Ranges::IntSet->new($rng);
	print commafy($exome->size()), "\n";


=head1 EXPORT

=over 5

=item * gdat_list_exons

=item * gdat_mrna_length

=item * gdat_list_introns

=item * gdat_is_coding

=item * gdat_list_coding_exons 

=item * gdat_cds_length

=item * gdat_list_5prime_utrs 

=item * gdat_list_3prime_utrs 

=item * gdat_list_utrs 

=item * gdat_to_bed12_line

=back

=head1 SUBROUTINES/METHODS

=head2 GeneTable Iterators

=over 5

=item * iter_geneTab

=item * iter_genePred

=back

=cut

my @GENETAB_FIELDS = qw|name chrom strand txStart txEnd cdsStart cdsEnd
					    exonCount exonStarts exonEnds|;

sub iter_geneTab {
	my ($db, $table) = @_;
	my $dbo = Genome::UCSC::DB->new($db);
	croak "$db.$table is not a gene table" unless _is_genetab($dbo, $table);
	my $iter = $dbo->iter("SELECT * FROM $table");
	return $iter;
}

sub _is_genetab {
  my ($dbo, $table) = @_;
  my $fieldInfo = $dbo->getFieldInfo($table);
  if (all { defined $fieldInfo->{$_} } @GENETAB_FIELDS) {
	return 1;
  }
}

my @GENEPRED_FIELDS = qw|bin name chrom strand txStart txEnd cdsStart cdsEnd
						 exonCount exonStarts exonEnds score name2
						 cdsStartStat cdsEndStat exonFrames|;

sub iter_genePred {
	my ($file) = @_;
	return iter_file($file, { header => 0, fields => \@GENEPRED_FIELDS });
}


=head1 Operations on Gene Data Object

Gene data object is a hash reference to each transcript with mandatory
gentable fields.

B<Note>: when numbering exon/intron, strand is not considered.
So the first exon for genes on negative strand would appear as the last one,
5'UTR for genes on negative strand should indeed be 3'UTR. To fully account 
for the strand issue, use L<Genome::UCSC::Transcript> module.

The genomic intervals returned by the following functions are all 1-based,
inclusive, conforming to the inhouse standard.

=over 5

=item * gdat_list_exons

=item * gdat_mrna_length

=item * gdat_list_introns

=item * gdat_is_coding

=item * gdat_list_coding_exons 

=item * gdat_cds_length

=item * gdat_list_5prime_utrs 

=item * gdat_list_3prime_utrs 

=item * gdat_to_bed12_line

=back

=cut

# List coding exons for one transcript
# Remember all of them are 0-based half-open coordinates!
# Transform them into 1-base, inclusive intervals!
sub gdat_list_exons {
  my ($dat) = @_;
  my @exonStarts = split(q|,|, $dat->{exonStarts});
  my @exonEnds   = split(q|,|, $dat->{exonEnds});
  my @exonRanges;
  for (my $jj = 0; $jj < @exonStarts; $jj ++) {
	next unless $exonStarts[$jj] < $exonEnds[$jj];
	push @exonRanges, [ $exonStarts[$jj]+1, $exonEnds[$jj],
						sprintf("%s.e%d", $dat->{name}, $jj+1) ];
  }
  return @exonRanges;
}

sub gdat_mrna_length {
  my ($dat) = @_;
  my $length;
  my @exonStarts = split(q|,|, $dat->{exonStarts});
  my @exonEnds = split(q|,|, $dat->{exonEnds});
  for ( my $ii = 0; $ii < $dat->{exonCount}; $ii ++ ) {
	$length += $exonEnds[$ii]-$exonStarts[$ii];
  }
  return $length;
}

sub gdat_list_introns {
  my ($dat) = @_;
  if ($dat->{exonCount} > 1) {
	my @exonStarts = split(q|,|, $dat->{exonStarts}); shift @exonStarts;
	my @exonEnds   = split(q|,|, $dat->{exonEnds}); pop @exonEnds;
	my @intronRanges;
	for (my $jj = 0; $jj < @exonStarts; $jj ++) {
	  next unless $exonEnds[$jj] < $exonStarts[$jj];
	  push @intronRanges, [ $exonEnds[$jj]+1, $exonStarts[$jj],
							sprintf("%s.i%d", $dat->{name}, $jj+1) ];
	}
	return @intronRanges;
  }
  else {
	return ();
  }
}

sub gdat_is_coding {
	my ($dat) = @_;
	return $dat->{cdsStart} < $dat->{cdsEnd};
}

sub gdat_list_coding_exons {
  my ($dat) = @_;
  if ($dat->{cdsStart} < $dat->{cdsEnd}) {
	my $all_exons = Array::IntSpan->new( gdat_list_exons($dat) );
	my $coding_exons = $all_exons->get_range($dat->{cdsStart}+1, $dat->{cdsEnd});
	return @$coding_exons;
  }
  else {
	return ();
  }
}

# Calculate coding sequence length.
sub gdat_cds_length {
  my ($dat) = @_;
  my $length;
  if ($dat->{cdsStart} < $dat->{cdsEnd}) {
	my @exonStarts = split(q|,|, $dat->{exonStarts});
	my @exonEnds = split(q|,|, $dat->{exonEnds});
	for ( my $ii = 0; $ii < $dat->{exonCount}; $ii ++ ) {
	  if ($exonEnds[$ii] < $dat->{cdsStart} ||
		  $exonStarts[$ii] > $dat->{cdsEnd} ) {
		next;
	  }
	  else {
		$length += min($exonEnds[$ii], $dat->{cdsEnd}) -
		  max($exonStarts[$ii], $dat->{cdsStart});
	  }
	}
	return $length;
  }
  else {
	return 0;
  }
}

# UTRs
sub gdat_list_5prime_utrs {
  my ($dat) = @_;
  if ($dat->{cdsStart} > $dat->{txStart}) {
	my $all_exons = Array::IntSpan->new( gdat_list_exons($dat) );
	my $five_utrs = $all_exons->get_range($dat->{txStart}+1, $dat->{cdsStart});
	return @$five_utrs;
  }
  else {
	return ();
  }
}

sub gdat_list_3prime_utrs {
  my ($dat) = @_;
  if ($dat->{cdsEnd} < $dat->{txEnd}) {
	my $all_exons = Array::IntSpan->new( gdat_list_exons($dat) );
	my $three_utrs = $all_exons->get_range($dat->{cdsEnd}+1, $dat->{txEnd});
	return @$three_utrs;
  }
  else {
	return ();
  }
}

sub gdat_list_utrs {
  my ($dat) = @_;
  my @three_utrs = gdat_list_5prime_utrs($dat);
  my @five_utrs = gdat_list_3prime_utrs($dat);
  my @utrs = (@five_utrs, @three_utrs);
  return @utrs;
}

# Convert data to UCSC bed12 line
sub gdat_to_bed12_line {
  my $dat = shift;
  my @exonStarts = split(q|,|, $dat->{exonStarts});
  my @exonEnds   = split(q|,|, $dat->{exonEnds});
	croak "Incorrect number of exons" unless
	  scalar(@exonStarts) == $dat->{exonCount} &&
		scalar(@exonEnds) == $dat->{exonCount};
  my @blockSizes = map { $exonEnds[$_]-$exonStarts[$_] } 0..$dat->{exonCount}-1;
  my @blockStarts = map { $_-$dat->{txStart} } @exonStarts;
  return join("\t", $dat->{chrom}, $dat->{txStart}, $dat->{txEnd}, $dat->{name}, 0,
			  $dat->{strand}, $dat->{cdsStart}, $dat->{cdsEnd}, 0,
			  $dat->{exonCount}, join(",", @blockSizes).",", join(",", @blockStarts).",",
			 );
}

=head1 AUTHOR

Xueya Zhou, C<< <xueyazhou at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-genome at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Genome>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Genome::UCSC::Genes


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

1; # End of Genome::UCSC::Genes

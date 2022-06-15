package Genome::UCSC::DB;

use strict;
use warnings;
use Carp;
use DBI;
use List::MoreUtils qw(any all notall none mesh);
use Genome::UCSC::BinKeeper ":all";
use Utils::Hash qw|merge_opts|;

our @EXPORT=qw(hAddBinToQuery hAddBinToQueryGeneral);
use base qw(Utils::DB Exporter);

=head1 NAME

Genome::UCSC::DB - Perl API to UCSC Database.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 DESCRIPTION

Inherit the L<Utils::DB> module, and adding custom methods to UCSC DB.

=head1 SUBROUTINES/METHODS

=head2 $class->new DB [, OPTIONS]

Create a new object given database name, with default host/user/pass info.

=cut

sub new {
	my ($class, $db, $argref) = @_;
	croak "Must provide database name!" unless defined $db;
	my $arg = merge_opts($argref, host => "genome-mysql.cse.ucsc.edu",
		db => $db, user => "genomep", pass => "password", port => 3306);
  	return $class->SUPER::new($arg);
}

=head2 Chromosome Info

Glean information from chromInfo and cytoBand tables

Error messages will be displayed for local database.

=cut

# Return some sequence named in chromInfo from the given db,
# or undef if db has no chromInfo.
sub hDefaultChrom {
	return $_[0]->getQuickScalar("SELECT chrom FROM chromInfo limit 1");
}

# Get list of all chromosome names in database.
sub hAllChromNames {
	return $_[0]->getQuickArray("SELECT chrom FROM chromInfo");
}


=head2 Split tables

Accomodate split tables.

=cut

my $HDB_MAX_SEQS_FOR_SPLIT = 100;
# Return TRUE if split tables are allowed in database.
sub hCanHaveSplitTables {
	return undef unless $_[0]->tableExists("chromInfo");
	my $count = $_[0]->getTableSize("chromInfo");
	return 1 if $count > 0 && $count < $HDB_MAX_SEQS_FOR_SPLIT;
}

# Retrieve (or build if necessary) the cached hash of split-consolidated
# tables for current db.
sub getTableHash {
  my ($self) = @_;
  if (!defined $self->{tabHash}) {
	my @allTables = $self->listTables;
	if ($self->hCanHaveSplitTables) {
	  foreach my $tbl (@allTables) {
		my ($trackName, $chrom) = $self->hParseTableName($tbl);
		push @{$self->{tabHash}{$trackName}}, $tbl;
	  }
	}
	else {
	  foreach my $tbl (@allTables) {
		$self->{tabHash}{$tbl} = [$tbl];
	  }
	}
  }
  return $self->{tabHash};
}

# Return TRUE if a table exists in db.
sub hTableExists {
  my ($self, $table) = @_;
  my $tabHash = $self->getTableHash;
  my ($trackName, $chrom);
  if (defined $self->hCanHaveSplitTables) {
	($trackName, $chrom) = $self->hParseTableName($table);
  }
  else {
	$trackName = $table;
  }
  if ( defined $tabHash->{$trackName} &&
	  (any { $_ eq $table } @{$tabHash->{$trackName}})
	 ) {
	return 1;
  }
}

# Returns all split tables for rootName or just rootName if not split
# or undef if no such table. Receive results in list context
sub hSplitTableNames {
  my ($self, $rootName) = @_;
  my $tabHash = $self->getTableHash;
  if (defined $tabHash->{$rootName}) {
	return @{$tabHash->{$rootName}};
  }
}

# Return TRUE if track table (or split table) exists
sub hTableOrSplitExists {
  my ($self, $trackName) = @_;
  my $tabHash = $self->getTableHash;
  return defined $tabHash->{$trackName};
}

# Parse an actual table name like "chr17_random_blastzWhatever" into
# the track name (blastzWhatever) and chrom (chr17_random).
sub hParseTableName {
  my ($self, $table) = @_;
  my ($trackName, $chrom) = ($table, $self->hDefaultChrom);
  if ($table =~ /^chr/ || $table =~ /^Group/) {
	if ($table =~ /^(\w+)_([^_]+)$/) {
	  ($trackName, $chrom) = ($2, $1);
	}
  }
  return ($trackName, $chrom);
}


=head2 TableInfo

Table information will be used in constructing range queries.
hFindTableInfo will return the  TableInfo structure, containing
following fields:

rootName : Name without chrN_
isPos    : True if table is positional
isSplit  : True if table is split by chromosome
hasBin   : True if table starts with binField
chromField : Name of chromosome field
startField : Name of chromosome start field
endField   : Name of chromosome end field
nameField  : Name of item name field
scoreField : Name of score field
strandField : Name of strand field
cdsStartField : cds(thick) Start field
cdsEndField   : cds(thick) End field
countField    : exon(block)Count field
startsField   : Name of exon(block)Starts field
endsSizesField : Name of exon(block)Ends(Sizes) field
spanField : Name of span field (wiggle)
hasCDS  : True if it has cdsStart,cdsEnd fields
hasBlocks : True if it has count,starts,endsSizes
type : A guess at the trackDb type for this table

=cut

my @table_info_fields =
  qw|rootName isPos isSplit hasBin chromField startField endField
	 nameField scoreField strandField cdsStartField cdsEndField
	 countField startsField endsSizesField hasCDS hasBlock type|;

# Given a table return the fields corresponding to all the bed 12
# fields, if they exist.  Fields that don't exist in the given table
# will be set to "" or undef

sub _fitField {
  my ($hash, $fieldName, $retField) = @_;
  if (defined $hash->{$fieldName}) {
	$$retField = $fieldName;
	return 1;
  }
  else {
	$$retField = "";
	return 0;
  }
}

sub hFindBed12FieldsAndBin {
  my ($self, $table, $hti) = @_;
  my $sr = $self->{dbh}->selectall_hashref("DESC $table","Field");
  my $gotIt = 1;
  if (defined $sr->{bin}) {
	$hti->{hasBin} = 1;
  }
  # Look for bed-style or linkedFeatures names.
  if (defined $sr->{chrom} && defined $sr->{chromStart} && defined $sr->{chromEnd}) {
	$hti->{chromField} = "chrom";
	$hti->{startField} = "chromStart";
	$hti->{endField} = "chromEnd";
	_fitField($sr, "name", \$hti->{nameField}) ||
	  _fitField($sr, "id", \$hti->{nameField}) ||
		_fitField($sr, "acc", \$hti->{nameField})  ||
		  _fitField($sr, "frag", \$hti->{nameField}) ||
			_fitField($sr, "config", \$hti->{nameField});
	_fitField($sr, "score", \$hti->{scoreField});
	_fitField($sr, "strand", \$hti->{strandField});
	_fitField($sr, "thickStart", \$hti->{cdsStartField});
	_fitField($sr, "thickEnd", \$hti->{cdsEndField});
	_fitField($sr, "blockCount", \$hti->{countField}) ||
	  _fitField($sr, "lfCount", \$hti->{countField});
	_fitField($sr, "chromStarts", \$hti->{startsField}) ||
	  _fitField($sr, "blockStarts", \$hti->{startsField}) ||
		_fitField($sr, "lfStarts", \$hti->{startsField});
	_fitField($sr, "blockSizes", \$hti->{endsSizesField}) ||
	  _fitField($sr, "lfSizes", \$hti->{endsSizesField});
	_fitField($sr, "span", \$hti->{spanField});
  }
  # Look for names of psl and psl-like (chain, chainLink, net, altGraphX,
  # some older types). */
  elsif (defined $sr->{tName} && defined $sr->{tStart} && defined $sr->{tEnd}) {
	$hti->{chromField} = "tName";
	$hti->{startField} = "tStart";
	$hti->{endField} = "tEnd";
	_fitField($sr, "qName", \$hti->{nameField}) ||
	  _fitField($sr, "name", \$hti->{nameField}) ||
		_fitField($sr, "chainId", \$hti->{nameField});
	_fitField($sr, "strand", \$hti->{strandField});
	_fitField($sr, "blockCount", \$hti->{countField});
	_fitField($sr, "tStarts", \$hti->{startsField});
	_fitField($sr, "blockSizes", \$hti->{endsSizesField});
  }
  # Look for gene prediction names
  elsif (defined $sr->{chrom} && defined $sr->{txStart} && defined $sr->{txEnd}) {
	$hti->{chromField} = "chrom";
	$hti->{startField} = "txStart";
	$hti->{endField} = "txEnd";
	_fitField($sr, "geneName", \$hti->{nameField}) ||
	  _fitField($sr, "name", \$hti->{nameField});
	_fitField($sr, "score", \$hti->{scoreField});
	_fitField($sr, "strand", \$hti->{strandField});
	_fitField($sr, "cdsStart", \$hti->{cdsStartField});
	_fitField($sr, "cdsEnd", \$hti->{cdsEndField});
	_fitField($sr, "exonCount", \$hti->{countField});
	_fitField($sr, "exonStarts", \$hti->{startsField});
	_fitField($sr, "exonEnds", \$hti->{endsSizesField});
  }
  # Look for repeatMasker names
  elsif (defined $sr->{genoName} && defined $sr->{genoStart} && defined $sr->{genoEnd}) {
	_fitField($sr, "repName", \$hti->{nameField});
	_fitField($sr, "swScore", \$hti->{scoreField});
	_fitField($sr, "strand", \$hti->{strandField});
  }
  # chromGraph table
  elsif (defined $sr->{chrom} && defined $sr->{chromStart} ) {
	$hti->{chromField} = "chrom";
	_fitField($sr, "chromStart", \$hti->{startField}) ||
	  _fitField($sr, "pos", \$hti->{startField});
	$hti->{endField} = $hti->{startField};
	_fitField($sr, "id", \$hti->{nameField}) ||
	  _fitField($sr, "name", \$hti->{nameField});
	_fitField($sr, "strand", \$hti->{strandField});
  }
  elsif ($table =~ /^chr/ && $table =~ /_gl$/ && defined $sr->{start} && defined $sr->{end}) {
	$hti->{chromField} = "";
	$hti->{startField} = "start";
	$hti->{endField} = "end";
	_fitField($sr, "frag", \$hti->{nameField});
	_fitField($sr, "strand", \$hti->{strandField});
  }
  else {
	_fitField($sr, "acc", \$hti->{nameField}) ||
	  _fitField($sr, "id", \$hti->{nameField}) ||
		_fitField($sr, "name", \$hti->{nameField});
	$gotIt = 0;
  }
  return $gotIt;
}

# Find table information, return undef if no table
sub hFindTableInfo {
  my ($self, $rootTable, $skip) = @_;
  my $hti;
  @{$hti}{@table_info_fields} = ('') x @table_info_fields;
  my $fullName;

  unless ($skip) {
	if ( !$self->hTableOrSplitExists($rootTable)) {
	  return undef;
	}
	else {
	  $hti->{isSplit} = 0;
	  $fullName = $rootTable;
	  my $chrom = $self->hDefaultChrom;
	  if ($chrom) {
		if ($self->hTableExists("${chrom}_${rootTable}")) {
		  $hti->{isSplit} = 1;
		  $fullName = "${chrom}_${rootTable}";
		}
	  }
	}
  }
  else {
	$fullName = $rootTable;
  }

  $hti->{isPos}  = $self->hFindBed12FieldsAndBin($fullName, $hti);
  $hti->{hasCDS} = 1 if $hti->{cdStartField};
  $hti->{hasBlocks} = 1 if $hti->{startsField};
  if ($hti->{isPos}) {
	if ( $hti->{startsField} eq 'exonStarts' ) {
	  $hti->{type} = 'genePred';
	} elsif ( $hti->{startsField} eq 'chromStarts' ||
			  $hti->{startsField} eq 'blockStarts' ) {
	  $hti->{type} = 'bed 12';
	} elsif ( $hti->{startsField} eq 'lfStarts' ) {
	  $hti->{type} = 'linkedFeatures';
	} elsif ( $hti->{startsField} eq 'tStarts' ) {
	  $hti->{type} = 'psl';
	} elsif ( $hti->{cdsStartField} ) {
	  $hti->{type} = 'bed 8';
	} elsif ( $hti->{strandField} && !$hti->{chromField} ) {
	  $hti->{type} = 'gl';
	} elsif ( $hti->{spanField} ) {
	  $hti->{type} = 'wiggle';
	} elsif ( $hti->{nameField} ) {
	  $hti->{type} = 'bed 4';
	} elsif ( $hti->{endField} ) {
	  $hti->{type} = 'bed 3';
	} elsif ( $hti->{startField} eq $hti->{endField} ) {
	  $hti->{type} = 'chromGraph';
	} else {
	  croak "Unrecognized table type for $fullName !";
	}
  }
  return $hti;
}


=head2 BinField

BinField related range query.

=cut

# If table contains bin field, produce SQL clause that will restrict to
# relevant bins to query.
sub _hAddBinToQueryStandard {
  my ($binField, $start, $end, $selfContained) = @_;
  my ($bFirstShift, $bNextShift) = (binFirstShift, binNextShift);
  my ($startBin, $endBin) = ($start, $end-1);
  $startBin >>= $bFirstShift;
  $endBin   >>= $bFirstShift;
  my $levels = binLevels;
  my $clause = $selfContained ? "(" : "";
  for (my $ii = 0; $ii < $levels; $ii ++) {
	my $offset = binOffset($ii);
	if ($ii > 0) {
	  $clause .= " OR ";
	}
	if ($startBin == $endBin) {
	  $clause .= sprintf("%s=%u", $binField, $startBin+$offset);
	}
	else  {
	  $clause .= sprintf("%s>=%u AND %s<=%u",
						 $binField, $startBin+$offset,
						 $binField, $endBin+$offset);
	}
	$startBin >>= $bNextShift;
	$endBin   >>= $bNextShift;
  }
  if ($selfContained) {
	$clause .= sprintf(" OR %s=%u )", $binField, binOffsetOldToExtended);
	$clause .= " AND ";
  }
  return $clause;
}

sub _hAddBinToQueryExtended {
  my ($binField, $start, $end) = @_;
  my ($bFirstShift, $bNextShift) = (binFirstShift, binNextShift);
  my $BINRANGE_MAXEND_512M = BINRANGE_MAXEND_512M;
  my ($startBin, $endBin) = ($start, $end-1);
  $startBin >>= $bFirstShift;
  $endBin   >>= $bFirstShift;
  my $levels = binLevelsExtended;

  my $clause = "(";
  if ($start < $BINRANGE_MAXEND_512M ) {
	$clause .= _hAddBinToQueryStandard($binField, $start, $BINRANGE_MAXEND_512M, undef);
	$clause .= " OR ";
  }
  for (my $ii = 0; $ii < $levels; $ii ++) {
	my $offset = binOffsetExtended($ii);
	if ($ii > 0) {
	  $clause .= " OR ";
	}
	if ($startBin == $endBin) {
	  $clause .= sprintf("%s=%u", $binField, $startBin+$offset);
	}
	else {
	  $clause .= sprintf("%s>=%u AND %s<=%u",
						 $binField, $startBin+$offset,
						 $binField, $endBin+$offset);
	}
	$startBin >>= $bNextShift;
	$endBin   >>= $bNextShift;
  }
  $clause .= ") AND ";
  return $clause;
}

sub hAddBinToQueryGeneral {
  my ($binField, $start, $end) = @_;
  if ( $end <= BINRANGE_MAXEND_512M ) {
	return _hAddBinToQueryStandard($binField, $start, $end, 1);
  }
  else {
	return _hAddBinToQueryExtended($binField, $start, $end);
  }
}

sub hAddBinToQuery {
  my ($start, $end) = @_;
  return hAddBinToQueryGeneral("bin",$start,$end);
}



=head1 AUTHOR

Xueya Zhou, C<< <xueyazhou at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-genome at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Genome>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Genome::UCSC::DB


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

1; # End of Genome::UCSC::DB

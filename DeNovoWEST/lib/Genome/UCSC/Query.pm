package Genome::UCSC::Query;

use strict;
use warnings;
use Carp;
use DBI;
use Iterator::Simple qw|iterator|;
use List::MoreUtils qw|all|;
use Genome::Ranges qw|validate_elem|;
use Genome::UCSC::DB qw(hAddBinToQuery hAddBinToQueryGeneral);
use Utils::Hash qw|merge_opts|;
use Data::Dumper;

=head1 NAME

Genome::UCSC::Query - Range query in UCSC database table.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

    use Genome::UCSC::Query;

    my $qr = Genome::UCSC::Query->new("hg19", {table => 'refGene', skip => 1});

    my $it = $qr->iter(spec_to_range("chr11:61,373,776-61,705,075"));
    while(my $dat = $it->()) {
        print join(" ", @{$dat}{qw|name name2 chrom txStart txEnd|}), "\n";
    }

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 $class->new DB, TABLE [, OPTIONS]

Create a database table range query object, with the following options

   where, having, order, limit  : SQL subclause to rootTable
   fields  :  fields to be selected, should be comma separated, default '*'
   join    :  combine information from other tables (see below)
              [NOTE] join operations are not supported for split tables!
   skip    :  skip checking existance, can speed up the process
   verbose :  output debug information

The structure for joining
   join => "left join kgXref on rootTable.name = kgXref.kgID"

If two or three tables are joined, we assume that the rootTable
is the main table. If we want to have its record will be found in the
output, left join to the rootTable will be used. The range query
is based on the positional information of rootTable.
Note the fields of the auxiliary tables will not be rigorously checked.

We also assume that the fields names from different tables will not
conflict (that is the same name will mean the same thing). Otherwise
the name from main table will override the name from link table,
you can use alias in the fields to change this default behavior.
The dbrow selectecd by iterator is represented by a hashref.

=cut

sub new {
	my ($class, $db, $table, $argref) = @_;

	my $arg = merge_opts($argref,
				where => undef,
				having   => undef,
				order => 1,
				skip => 1, # skip checking existance, speed up the process
				limit => undef,
				join  => "",
				fields => '*', # comma separated field list
				verbose => 1, # output debug information
				user => undef, pass => undef, port => undef, # additional DB parameters
				);

	my %dbarg = map { $_ => $arg->{$_} } grep { defined $arg->{$_} } qw|user pass port|;
	my $hdb = Genome::UCSC::DB->new($db, \%dbarg);

	my $hti   = $hdb->hFindTableInfo($table, $arg->{skip});

	if (!defined $hti) {
		warn "Table $table does not exist or hFinTableInfo failed";
		return;
	}
	if (!$hti->{isPos}) {
		warn "Table $table is not positional!";
		return;
	}

	my ($SELECT, $FROM, $WHERE1, $WHERE2, $HAVING, $EXTRA, $ORDER);
	$SELECT = "SELECT $arg->{fields} ";
	$FROM   = "FROM " . ($hti->{isSplit} ? "SPLITTABLE AS MAIN " :
		"$table AS MAIN $arg->{join} ");
	$WHERE1 = ($hti->{isSplit} ? "" : "WHERE MAIN.chrom=CHROMOSOME ");
	$WHERE2 = ($WHERE1 eq "" ? "WHERE " : "AND ") .
	($hti->{hasBin} ? "ADDBINTOQUERY " : "") .
	"MAIN.$hti->{endField} >= CHROMSTART AND MAIN.$hti->{startField} < CHROMEND ";
	$EXTRA  = ($arg->{where} && $arg->{where} !~ /^order|limit/ ?
		"$arg->{where} " : "");
	$HAVING = $arg->{having} ? $arg->{having} : "";
	$ORDER  = ($arg->{order} ? "ORDER BY ".
		( $hti->{isSplit} ? "$hti->{startField} " :
			"$hti->{chromField}, $hti->{startField} " ) : "" ) .
	( $arg->{limit} ? "limit $arg->{limit}" : "" );

	my $self = bless { hdb => $hdb, split => $hti->{isSplit} ? "Split" : "",
	bin => $hti->{hasBin} ? "Bin" : "", verbose => $arg->{verbose} }, $class;
	if ( $self->isSplit ) {
		if ( !$self->hasBin ) {
			$WHERE2 =~ s/CHROMSTART/?/;
			$WHERE2 =~ s/CHROMEND/?/;
		}
		for my $tab ($hdb->hSplitTableNames($table)) {
			my (undef, $chrom) = $hdb->hParseTableName($tab);
			(my $from = $FROM) =~ s/SPLITTABLE/$tab/;
			$self->{query}{$chrom}   = "$SELECT $from $WHERE1 $WHERE2 ";
			$self->{querynr}{$chrom} = "$SELECT $from $WHERE1 ";
			if ($EXTRA) {
				if ($self->{query}{$chrom} =~ /WHERE/) {
					$self->{query}{$chrom} .= "AND $EXTRA ";
				}
				else {
					$self->{query}{$chrom} .= "WHERE $EXTRA ";
				}
				if ($self->{querynr}{$chrom} =~ /WHERE/) {
					$self->{querynr}{$chrom} .= "AND $EXTRA ";
				}
				else {
					$self->{querynr}{$chrom} .= "WHERE $EXTRA ";
				}
			}
			$self->{query}{$chrom} .= "$HAVING $ORDER";
			$self->{querynr}{$chrom} .= "$HAVING $ORDER";
		}
	}
	else {
		if (! $self->hasBin) {
			$WHERE1 =~ s/CHROMOSOME/?/;
			$WHERE2 =~ s/CHROMSTART/?/;
			$WHERE2 =~ s/CHROMEND/?/;
		}
		$self->{query} = "$SELECT $FROM $WHERE1 $WHERE2";
		$self->{querynr} = "$SELECT $FROM $WHERE1 ";
		if ($EXTRA) {
			if ($self->{query} =~ /WHERE/) {
				$self->{query} .= "AND $EXTRA ";
			}
			else {
				$self->{query} .= "WHERE $EXTRA ";
			}
			if ($self->{querynr} =~ /WHERE/) {
				$self->{querynr} .= "AND $EXTRA ";
			}
			else {
				$self->{querynr} .= "WHERE $EXTRA ";
			}
		}
		$self->{query} .= "$HAVING $ORDER";
		$self->{querynr} .= "$HAVING $ORDER";
	}

	return $self;
}

# Helper functions

sub hasBin  {  $_[0]{bin} ne ""   }

sub isSplit {  $_[0]{split} ne "" }

sub dbh     {  $_[0]{hdb}{dbh}    }

sub verbose { $_[0]{verbose} }

# Four internal query types:
# Split or not x Bin or not
sub type   {  $_[0]{split}.$_[0]{bin} }


=head2 $self->prepare CHROM, START, END

Prepare range query. The function accomdates different data table types.

When only chromosome is provided, then prepare bulk query for the chromosome.

Return a database statement handler.

The chromosome range will be 1-based and inclusive. It will internally converted
to UCSC's convention.

=cut

sub prepare {
	my ($self,$chrom,$start,$end) = @_;
	$end = $start if defined $start && !defined $end;
	validate_elem($chrom,$start,$end);
	my $pMethod = "_prepare".$self->type;
	my $sth = $self->$pMethod($chrom,$start,$end);
	$sth->execute();
	return $sth;
}

sub _prepareSplit {
	my ($self, $chrom, $start, $end) = @_;
	my $sth;
	if (defined $start && defined $end) {
		print STDERR $self->{query}{$chrom}, "\n" if $self->verbose;
		$sth = $self->dbh->prepare($self->{query}{$chrom});
		$sth->bind_param(1, $start);
		$sth->bind_param(2, $end  );
	}
	else {
		print STDERR $self->{querynr}{$chrom}, "\n" if $self->verbose;
		$sth = $self->dbh->prepare($self->{querynr}{$chrom});
	}
	return $sth;
}

sub _prepareSplitBin {
	my ($self, $chrom, $start, $end) = @_;
	my ($query, $sth);
	if (defined $start && defined $end) {
		$query = $self->{query}{$chrom};
		$query =~ s/CHROMSTART/$start/;
		$query =~ s/CHROMEND/$end/;
		$query =~ s/ADDBINTOQUERY/hAddBinToQueryGeneral("MAIN.bin", $start, $end)/e;
		$sth = $self->dbh->prepare($query);
	}
	else {
		$query = $self->{querynr}{$chrom};
	}
	print STDERR $query, "\n" if $self->verbose;
	$sth = $self->dbh->prepare($query);
	return $sth;
}

sub _prepare {
	my ($self, $chrom, $start, $end) = @_;
	my $sth;
	if (defined $start && defined $end) {
		print STDERR $self->{query}, "\n" if $self->verbose;
		$sth = $self->dbh->prepare($self->{query});
		$sth->bind_param(1, $chrom);
		$sth->bind_param(2, $start);
		$sth->bind_param(3, $end  );
	}
	else {
		print STDERR $self->{querynr}, "\n" if $self->verbose;
		$sth = $self->dbh->prepare($self->{querynr});
		$sth->bind_param(1, $chrom);
	}
	return $sth;
}

sub _prepareBin {
	my ($self, $chrom, $start, $end) = @_;
	my $query;
	if (defined $start && defined $end) {
		$query = $self->{query};
		$query =~ s/CHROMOSOME/'$chrom'/;
		$query =~ s/CHROMSTART/$start/;
		$query =~ s/CHROMEND/$end/;
		$query =~ s/ADDBINTOQUERY/hAddBinToQueryGeneral("MAIN.bin",$start, $end)/e;
	}
	else {
		$query = $self->{querynr};
		$query =~ s/CHROMOSOME/'$chrom'/;
	}
	print STDERR $query, "\n" if $self->verbose;
	my $sth = $self->dbh->prepare($query) or croak $DBI::errstr;
	return $sth
}

=head2 quick_query CHROM, START, END

Quick range query (chrom or chrom,start,end) : return an array of hashrefs

=cut

sub quick_query {
	my $self = shift @_;
	validate_elem(@_);
	my $sth = $self->prepare(@_);
	return $sth->fetchall_arrayref( {} );
}


=head2 iterator RANGES [, OPTIONS]

Given a list of ranges, execute range query for each given range.
Return an iterator to the results.

The argument is an arrayref of arrayrefs [chrom] or [chrom,start,end].

The return value of the iterator is hash/hashref (default) or array/arrayref,
depending on the context and can be controled by C<array> option.

=cut

sub iter {
	my ($self, $ranges, $argref) = @_;
	my $arg = merge_opts($argref, array => undef);

	my @ranges;

	if (ref $ranges eq 'ARRAY') {
		if (ref $ranges->[0] eq 'ARRAY') {
			@ranges = $ranges;
		}
		else {
			@ranges = ( [ $ranges->[0] ] );
		}
	}
	else {
		@ranges = map { [ $_ ] } $self->{hdb}->hAllChromNames;
	}
	foreach my $intv (@ranges) {
		print STDERR join(" ", @$intv), "\n" if $self->verbose;
		if (@$intv > 1) {
			validate_elem(@$intv);
		}
		else {
			croak "Invalid chromosome name" unless $intv->[0] =~ /^\w+$/;
		}
	}

	my $ii = 0;
  	my $sth = $self->prepare(@{$ranges->[0]}); # statement handle

  	if ($arg->{array}) {
  		return iterator {
  			my $dbrow = $sth->fetchrow_arrayref;
  			if ( !defined $dbrow ) {
  				if ($sth->err) {
  					croak "fetch error : " . $sth->errstr;
  				}
  				while ( !defined $dbrow ) {
  					if ( $ii < $#$ranges ) {
  						$sth   = $self->prepare(@{$ranges->[++$ii]});
  						$dbrow = $sth->fetchrow_arrayref;
  						if (wantarray) {
  							return @$dbrow;
  						}
  						else {
  							return $dbrow;
  						}
  					}
  					else {
  						return;
  					}
  				}
  			}
  			if (wantarray) {
  				return @$dbrow;
  			}
  			else {
  				return $dbrow;
  			}
  		};
  	}
  	else {
  		return iterator {
  			my $dbrow = $sth->fetchrow_hashref;
  			if ( !defined $dbrow ) {
  				if ($sth->err) {
  					croak "fetch error : " . $sth->errstr;
  				}
  				while ( !defined $dbrow ) {
  					if ( $ii < $#$ranges ) {
  						$sth   = $self->prepare(@{$ranges->[++$ii]});
  						$dbrow = $sth->fetchrow_hashref;
  						if (wantarray) {
  							return %$dbrow;
  						}
  						else {
  							return $dbrow;
  						}
  					}
  					else {
  						return;
  					}
  				}
  			}
  			if (wantarray) {
  				return %$dbrow;
  			}
  			else {
  				return $dbrow;
  			}
  		};
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

    perldoc Genome::UCSC::Query


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

1; # End of Genome::UCSC::Query

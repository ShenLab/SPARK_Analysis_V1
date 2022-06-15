package Utils::DB;

use strict;
use warnings;
use Carp;
use DBI;
use Iterator::Simple qw|iterator|;
use List::MoreUtils qw|all|;
use Utils::Iterator;
use Utils::Hash qw|merge_opts|;


=head1 NAME

Utils::DB - Utility interface for MySQL database.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

	use Utils::DB;

	my $dbo = Utils::DB->new({ host => 'genome-mysql.cse.ucsc.edu', db => 'hg19',
		user => 'genomep', pass => 'password' });
	my $ntrans = $dbo->getQuickScalar("SELECT count(*) FROM refGene WHERE cdsStart<txStart");
	my @dmd = $dbo->getQuickScalar("SELECT * FROM from refGene WHERE name2='DMD'");

	my $git = igrep { $_->{cdsStart}<$_->{txStart} } $dbo->iter();
	while($gdat = $git->()) {
		...
	}


=head1 SUBROUTINES/METHODS

=head2 $self->new OPTIONS

Create a new object.
Required input options: host, db, user, pass; Optional: port.

=cut

sub new {
	my ($class, $argref) = @_;
	my $arg = merge_opts($argref, host => undef, db => undef,
		user => undef, pass => undef, port => 3306);
	my %self = (%$arg);
	$self{dbh} = _dbconnect(%$arg);
	return bless \%self, $class;
}

sub _dbconnect {
	my %arg = @_;
	my $dbh = DBI->connect("DBI:mysql:host=$arg{host};database=$arg{db};port=$arg{port}",
		$arg{user}, $arg{pass})	or croak("Failed connecting to the database on $arg{host}");
	return $dbh;
}

=head2 $self->getQuickScalar SQL

Return the scalar result of single row/single column query.

=cut

sub getQuickScalar {
  my ($self, $sql) = @_;
  my $res = $self->{dbh}->selectall_arrayref($sql);
  if (defined $res && @$res > 0) {
	return $res->[0][0];
  }
  else {
	return undef;
  }
}

=head2 $self->getQuickArray SQL

Return the array for single column query or just keep the
first field for each row.

=cut 

sub getQuickArray {
  my ($self, $sql) = @_;
  return map { $_->[0] } @{ $self->{dbh}->selectall_arrayref($sql) // () };
}

=head2 $self->getQuickHash SQL

Return a hash for two column query.

=cut

sub getQuickHash {
  my ($self, $sql) = @_;
  return map { $_->[0] => $_->[1] } @{ $self->{dbh}->selectall_arrayref($sql) // () };
}

=head2 $self->getQuickArrayRef SQL

Return the an array of hashref for multi-column query

=cut

sub getQuickArrayRef {
  my ($self, $sql) = @_;
  return @{ $self->{dbh}->selectall_arrayref($sql, { Slice => {} }) // () };
}


=head2 $self->listTables [PATTERN]

List table names in the database.

=cut

sub listTables {
  my ($self, $pattern) = @_;
  if (defined $pattern) {
  	return $self->getQuickArray("SHOW TABLES LIKE '$pattern'");
  }
  else {
  	return $self->getQuickArray("SHOW TABLES");
  }
}

=head2 $self->listFields TABLE

Return a list of fields in the table.

=cut

sub listFields {
  my ($self, $table) = @_;
  return $self->getQuickArray("DESC $table");
}


=head2 $self->tableExists TABLE

Test if a table exists
May be slow if too many tables exist

=cut

sub tableExists {
  my ($self, $table) = @_;
  return $self->getQuickArray(<<EOT
     SELECT DISTINCT TABLE_NAME FROM information_schema.COLUMNS
     WHERE TABLE_SCHEMA = '$self->{db}' AND TABLE_NAME = '$table'
EOT
               );
}

=head2 $self->tablesExist TABLE1,TABL2,...

Test if all tables exists

=cut 

sub tablesExist {
  my $self = shift;
  return 1 if ( all { $self->tableExists($_) > 0 } @_ );
}

=head2 $self->tableWildExists PATTERN

Return TRUE if table (which can include SQL wildcards) exists.

=cut 

sub tableWildExists {
  my ($self, $table) = @_;
  return 1 if $self->getQuickArray("SHOW TABLES LIKE '$table'") > 0;
}

=head2 $self->getColumnNum TABLE

Return the number of columns in a table.

=cut

sub getColumnNum {
  my ($self, $table) = @_;
  return scalar $self->listFields($table);
}

=head2 $self->getTableSize TABLE

Return row count if a table exists, undef if it does not

=cut

sub getTableSize {
  my ($self, $table) = @_;
  return $self->getQuickScalar("SELECT count(*) FROM $table");
}

=head2 $self->getFieldInfo TABLE

Return field information for table.

=cut

sub getFieldInfo {
  my ($self, $table) = @_;
  return $self->{dbh}->selectall_hashref("SHOW COLUMNS FROM $table","Field");
}

=head2 $self->getFieldIndex TABLE, FIELD

Returns index of field in a row from table, or undef if it
does not exist

=cut

sub getFieldIndex {
  my ($self, $table, $field) = @_;
  my @fnames = $self->listFields;
  for ( my $ii = 0; $ii < @fnames; $ii ++ ) {
	if ( $fnames[$ii] eq $field ) {
	  return $ii;
	}
  }
}


=head2 $self->getPrimaryKey TABLE

Get primary key if any for table return empty array if none

=cut

sub getPrimaryKey {
  my ($self, $table) = @_;
  return map { $_->{Field} }
	grep { $_->{Key} eq 'PRI' }
	  @{ $self->{dbh}->selectall_arrayref("DESC $table",
										  { Slice =>{} } ) // () };
}


=head2 $self->getEnumDef TABLE, COLNAME

Get the definitions of a enum column in a table, returning an
array of enum values.

=cut
 
sub getEnumDef {
  my ($self, $table, $colname) = @_;
  my $res = $self->{dbh}->selectall_hashref("DESC $table", "Field");
  if (defined $res->{$colname} && $res->{$colname}{Type} =~ /^set|enum/) {
	my $def = $res->{$colname}{Type};
	$def =~ s/^set|enum\((.+)\)$/$1/;
	$def =~ s/['\\]//g;
	return map { $_ =~ s/'//g; $_ } split(',', $def)
  }
}

=head2 $self->getRandSample TABLE, NUM, SEED

Get random sample from a database table given seed

=cut

sub getRandSample {
  my ($self, $table, $num, $seed) = @_;
  my $res = $self->{dbh}->selectall_arrayref(<<EOF
      SELECT * FROM $table ORDER BY RAND($seed) LIMIT $num
EOF
											 ,{ Slice =>{} }
											);
  return $res;
}


=head2 $self->iter SQL [, OPTIONS]

Return an L<Iterator::Simple> iterator.

See also L<Iterator::DBI> pacakge.

Options:
By default, the return value by calling iterator is hash or hashref (depending on the context),
users can ask it to return arrayref by switching on the C<array> option.

=cut

sub iter {
  my ($self, $sql, $argref) = @_;
  #return idb_rows($self->{dbh}, $sql);

  my $arg = merge_opts($argref, array => undef);

  my $sth = $self->{dbh}->prepare($sql);
  $sth->execute();

  if ($arg->{array}) {
    return iterator {
      my $dbrow = $sth->fetchrow_arrayref;
      if ( !defined $dbrow ) {
        return;
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
        return;
      }
      if (wantarray) {
        return %{$dbrow};
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

Please report any bugs or feature requests to C<bug-utils at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Utils>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Utils::DB


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

1; # End of Utils::DB

package Genome::Ensembl;

use strict;
use warnings;
use Carp;
use List::Util qw|max|;
use List::MoreUtils qw|all|;
use Data::Dumper;
use Utils::Hash qw|merge_opts|;
use Net::FTP;

use base qw|Exporter|;
our @EXPORT_OK = qw|load_ensembl curr_ver|;
our %EXPORT_TAGS = ( 'all' => \@EXPORT_OK );


=head1 NAME

Genome::Ensembl - Utility functions to Ensembl database!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

    use Genome::Ensembl;

    my $reg = load_ensembl({ version => 75 });

    my $gene_adaptor = $reg->get_adaptor('human', 'core', 'gene');
    my @genes = @{$gene_adaptor->fetch_all};

    while(my $gene = shift @genes) {
        print $fout $gene->stable_id, "\t",
              $gene->canonical_transcript->stable_id, ".", 
              $gene->canonical_transcript->version, "\n";
    }


=head1 EXPORT

=over 5

=item * load_ensembl

=item * curr_ver

=head1 SUBROUTINES/METHODS

=head2 load_ensembl APIDIR [, OPTIONS]

Add ensembl API to PERL5LIB path during run time, and return database registry.

Options:

* basedir : base directory for bioperl and ensembl modules default to ~/lib/EnsAPI

* version : the release version used, we will match the database version with API,
            default is the latest release returned from C<curr_ver>

* species : which species database will be used, default: homo_sapiens

=cut

sub load_ensembl {
	my ($argref) = @_;
	my $arg = merge_opts($argref, basedir => "$ENV{HOME}/Software/EnsAPI",
		version => undef, species => "homo_sapiens");
	my $version = $arg->{version};
	if (!defined $arg->{version}) {
		$version = curr_ver();
	}

	my @modules = qw|ensembl ensembl-compara ensembl-variation ensembl-funcgen ensembl-io|;
	croak "Must provide basedir for Ensembl API" unless defined $arg->{basedir} && -d $arg->{basedir} &&
		(all { -d "$arg->{basedir}/release$version/$_" } @modules);
	croak "Cannot find BioPerl-1.6.1 in $arg->{basedir}" unless -d "$arg->{basedir}/BioPerl-1.6.1";

	push @INC, "$arg->{basedir}/BioPerl-1.6.1";
	foreach my $subdir (@modules) {
		push @INC, "$arg->{basedir}/release$version/$subdir/modules";
	}

	require Bio::EnsEMBL::Registry;
	my $reg = "Bio::EnsEMBL::Registry";
	$reg->load_registry_from_db (
		-host => 'ensembldb.ensembl.org',  -user => 'anonymous',
 		-species  => $arg->{species}, -db_version => $version, -verbose => 1);
	return $reg;
}


=head2 curr_ver

Return the current version of ensembl database.

The version is obtained by listing the directory in the FTP.
ftp://ftp.ensembl.org/pub/


=cut

sub curr_ver {
	my $ftp_site = 'ftp.ensembl.org';
    my $ftp = Net::FTP->new($ftp_site) 
        or croak "Could not connect to $ftp_site: $!";
    $ftp->login("anonymous",'-anonymous@')
      or croak "Cannot login ", $ftp->message;
    $ftp->cwd("/pub")
      or croak "Cannot change working directory ", $ftp->message;

    my @remote_files = $ftp->ls();
    my @rels = map { $_ =~ s/release\-//; $_ } grep { /release-(\d+)/ } @remote_files;
    return max(@rels);
}


=head2 legacy_date VERSION

Return the date for the legacy version.
In scalar context, return "jul2015" useful for construction URL.
In array context, return ("Jul", 2015).

Ensembl 91: Dec 2017
Ensembl 90: Aug 2017
Ensembl 89: May 2017
Ensembl 88: Mar 2017
Ensembl 87: Dec 2016
Ensembl 86: Oct 2016
Ensembl 85: Jul 2016
Ensembl 84: Mar 2016
Ensembl 83: Dec 2015
Ensembl 82: Sep 2015
Ensembl 81: Jul 2015
Ensembl 80: May 2015
Ensembl 79: Mar 2015
Ensembl 78: Dec 2014
Ensembl 77: Oct 2014
Ensembl 76: Aug 2014
Ensembl 75: Feb 2014
Ensembl 74: Dec 2013
Ensembl 67: May 2012
Ensembl 54: May 2009

=cut

my $DATE = <<EOF;
Ensembl 91: Dec 2017
Ensembl 90: Aug 2017
Ensembl 89: May 2017
Ensembl 88: Mar 2017
Ensembl 87: Dec 2016
Ensembl 86: Oct 2016
Ensembl 85: Jul 2016
Ensembl 84: Mar 2016
Ensembl 83: Dec 2015
Ensembl 82: Sep 2015
Ensembl 81: Jul 2015
Ensembl 80: May 2015
Ensembl 79: Mar 2015
Ensembl 78: Dec 2014
Ensembl 77: Oct 2014
Ensembl 76: Aug 2014
Ensembl 75: Feb 2014
Ensembl 74: Dec 2013
Ensembl 67: May 2012
Ensembl 54: May 2009
EOF

sub legacy_date {
	my ($version) = @_;
	my %DATE = map { my ($k, $v) = split(/:\s+/, $_);
				     my @v = split(/\s+/, $v); 
				     $k => [@v] 
				    } split(/\n/, $DATE);
	if (defined $DATE{"Ensembl $version"}) {
		if (wantarray) {
			return @{$DATE{"Ensembl $version"}};	
		}
		else {
			return join("", map { lc($_) } @{$DATE{"Ensembl $version"}});
		}
	}
	else {
		return;
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

    perldoc Genome::Ensembl


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

1; # End of Genome::Ensembl

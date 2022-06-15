package Genome::UCSC::Liftover;

use strict;
use warnings;
use Carp;
use Data::Dumper;
use File::Basename;
use Utils::File qw|open_file|;
use Utils::Hash qw|merge_opts set_default chk_default|;
use IntervalTree;

use base qw|Exporter|;

our @EXPORT_OK = qw|check_chain_chrpref|;


=head1 NAME

Genome::UCSC::Liftover - Perl implementation of liftover genome coordinate conversion.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

	use Genome::UCSC::Liftover;

	my $lift = Genome::UCSC::Liftover->new("hg19ToHg38.over.chain.gz");

	open my $fin, "gzcat CDH_MGH_VQSR.vcf.gz | grep -v '^#' | cut -f1-5 |" or croak "Cannot open pipe";
	my $fbed = IO::File->new("markers.bed", "w");
	my $fout = IO::File->new("markers.hg38.lifted", "w");
	my $fump = IO::File->new("markers.hg38.unmapped", "w");
	while(<$fin>) {
		chomp;
		my ($chrom, $pos, $name, $ref, $alt) =  split;
		print $fbed join("\t", "chr".$chrom, $pos-1, $pos, "$chrom:$pos", 100, "+"), "\n";
		my @res = $lift->query("chr".$chrom, $pos);
		if (@res == 1) {
			my ($tg_chrom, $tg_pos, $tg_strand, $tg_score) = @{$res[0]};
				print $_ ,"\t", join("\t", $tg_chrom, $tg_pos, $tg_strand), "\n";
			}
		}
		else {
			print $fump join("\t", $chrom, $pos, "+"), "\n";
		}
	}


=head1 SUBROUTINES/METHODS

=head2 Utility function

check_chain_chrprefix: check if query and target genome assembly use 'chr' prefix in chromsome names.

=cut


sub check_chain_chrpref {
	my ($chainpath) = @_;
	my ($qchr, $tchr);
	my $fch;
	if ($chainpath =~ /\.gz$/) {
		open $fch, "<:gzip", $chainpath or die "Cannot open gzipped chain";
	}
	else {
		open $fch, $chainpath or die "Cannot open chain";
	}
	my $header = <$fch>;
	my ($tName, $qName) = (split(/\s+/, $header))[2,7];
	if ($tName =~ /^chr/) {
		$tchr = 1;
	}
	else {
		$tchr = 0;
	}
	if ($qName =~ /^chr/) {
		$qchr = 1;
	}
	else {
		$qchr = 0;
	}
	return ($qchr, $tchr);
}



=head2 $class->new CHAINFILE

Create an object from chain file.

=cut

my @chain_fields = qw|score source_name source_size source_strand source_start source_end
	target_name target_size target_strand target_start target_end id|;

sub new {
	my ($class, $chain_file) = @_;

	# Load chains
	my @chains;
	my $fin = open_file($chain_file);
	my %chain;
	my ($sfrom, $tfrom);
	my $prencol ;
	while(<$fin>) {
		next if /^#/ || /^\s+$/;
		my @vals = split;
		if (@vals == 12 || @vals == 13) {
			croak "Incorrect chain line\n$_" unless $vals[0] eq 'chain';
			# chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 1
			if (defined $prencol) {
				croak "Expecting one number on the end of last alignment block" unless $prencol == 1;
				push @chains, { %chain };
			}
			shift @vals;
			push @vals, undef if @vals == 11;
			%chain = ();
			@chain{@chain_fields} = @vals;
			$chain{blocks} = [];
			($sfrom, $tfrom) = ($chain{source_start}, $chain{target_start});
		}
		elsif (@vals == 3){
			my ($size, $sgap, $tgap) = @vals;
			push @{$chain{blocks}}, [$sfrom, $sfrom+$size, $tfrom, $tfrom+$size];
			$sfrom += $size + $sgap;
			$tfrom += $size + $tgap;
		}
		elsif (@vals == 1) {
			my $size = $vals[0];
			push @{$chain{blocks}}, [$sfrom, $sfrom+$size, $tfrom, $tfrom+$size];
			if ($sfrom+$size != $chain{source_end} ||
				$tfrom+$size != $chain{target_end}) {
				croak "Alignment blocks do not match specified block sizes."
			}

		}
		$prencol = @vals;
	}
	# Create chain index from list of chains
	my %chain_index;
	my (%source_size, %target_size);
	foreach my $c (@chains) {
		chk_default(\%source_size, $c->{source_name}, $c->{source_size});
		chk_default(\%target_size, $c->{target_name}, $c->{target_size});
		set_default(\%chain_index, $c->{source_name}, IntervalTree->new());
		my $tree = $chain_index{$c->{source_name}};
		foreach my $bk (@{$c->{blocks}}) {
			$tree->insert_interval(IntervalTree::Interval->new($bk->[0], $bk->[1], 
				[ $bk->[2], $bk->[3], $c ]));
		}
	}
	return bless { file => $chain_file, index => \%chain_index}, $class;
}

=head2 $self->query CHROM, POS [,OPTIONS]

Returns a B<list> of possible conversions for a given chromosome position.
      
If chromosome is completely unknown to the LiftOver, none is returned.
        
B<NOTE> that coordinates stored in chain files are 0-based, and even at negative strand
are relative to the beginning of the genome. I.e. position 0 strand + is the first position
of the genome. Position 0 strand - is also the first position of the genome (and the last
position of reverse-complemented genome). 

The START and END provided in the argument are 1-based inclusive, compatible with
L<Genome::Ranges> module.


Options
C<unique>: when unique is set to 1 (default), only uniquely mapped position will be returned.
Otherwise, all positions will be retured (sorted by chain scores).

=cut

sub query {
	my ($self, $chrom, $position, $argref) = @_;
	my $arg = merge_opts($argref, unique => 1);

	my $index = $self->{index};
	return unless defined $index->{$chrom};
	# Convert to internal representation
	$position -= 1;

	# First query for the chains, sorted by decreasing chain scores
	my $overlaps = $index->{$chrom}->find($position, $position+1);
	return unless @$overlaps;
	
	my @results;
	foreach my $block (@$overlaps) {
		my ($source_start, $source_end, $data) =
			($block->{start}, $block->{end}, $block->{value});
		my ($target_start, $target_end, $chain) = @$data;
		my $result_pos = $target_start + ($position - $source_start);
		if ($chain->{target_strand} eq '-') {
			$result_pos = $chain->{target_size} - 1 - $result_pos;
		}
		push @results, [$chain->{target_name}, $result_pos+1, 
			$chain->{target_strand}, $chain->{score}];

	}
	if (@results > 1) {
		@results = sort {  $b->[3] <=> $a->[3] } @results;
	} 
	unless ($arg->{unique}) {
		if (wantarray) {
			return @results;
		}
		else {
			return \@results;
		}
	}
	else {
		if (@results == 1) {
			if (wantarray) {
				return @{$results[0]};
			}
			else {
				return $results[0];
			}
		}
		else {
			return;
		}
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

    perldoc Genome::UCSC::Liftover


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

1; # End of Genome::UCSC::Liftover

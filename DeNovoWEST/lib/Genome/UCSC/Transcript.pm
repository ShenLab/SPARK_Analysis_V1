package Genome::UCSC::Transcript;

use strict;
use warnings;
use Carp;
use IO::File;
use Cwd;
use Set::Scalar;
use Storable qw(dclone);
use Data::Dumper;
use List::MoreUtils qw(all notall any none uniq mesh pairwise);
use Set::IntSpan;
use Array::IntSpan;
use Utils::Seq qw|dna2peptide rev_comp|;
use Genome::Ranges qw|range_to_spec|;


=head1 NAME

Genome::UCSC::Transcript - The great new Genome::UCSC::Transcript!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use Genome::UCSC::Transcript;

    my $foo = Genome::UCSC::Transcript->new();
    ...

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS


=head2 Constructor

  Construct a transcript object from one row of UCSC gene tables or gene
  prediction format files. Each row is read as a hashref contain at least
  the following fields: name, strand, chrom, (txStart/txEnd [optional]), 
  cdsStart/cdsEnd, and exonStarts/exonEnds

  The transcript object is internally represented by a hashref, with the
  fields of value type: _name, _strand, _chrom, _txStart, _txEnd, _cdsStart,
  _cdsEnd.  And two fields of interlaced intervals of introns and exons.
  If cds exist, we add fields for coding_exon and 3 and 5 prime UTRs. These
  range objects are internally represented by Array::IntSpan objects.
  Note that our definition of coding exon do not contain UTR portion.

  Necessary work has been done to deal with subtleties of strandness and half-open
  zero-based range in the UCSC data tables. All the internal representation use
  one-based inclusive ranges.

=cut

sub new {
	my ($class, $row) = @_;
	my %gene =
	( _name   => $row->{name}     || do { croak "name field required!" },
		_strand => $row->{strand} eq '+' ?  1 :
		$row->{strand} eq '-' ? -1 :
		do { croak "unrecognized strand!" },
		_chrom  => $row->{chrom}    || do { croak "chrom field required!" },
	  # Change cdsStart/End txStart/End to their normal meanings
	  # But, it may causes some problems to determine the range of a gene
	  	_txStart  => !defined $row->{txStart} ? undef :
	  		( $row->{strand} eq '+' ? $row->{txStart}+1 : $row->{txEnd} ),
	  	_txEnd    => !defined $row->{txEnd}   ? undef :
	  		( $row->{strand} eq '+' ? $row->{txEnd} :  $row->{txStart}+1 ),
	  	_cdsStart => $row->{cdsStart} == $row->{cdsEnd} ? undef :
	  		( $row->{strand} eq '+' ? $row->{cdsStart}+1 : $row->{cdsEnd} ),
	  	_cdsEnd   => $row->{cdsStart} == $row->{cdsEnd} ? undef :
	  		( $row->{strand} eq '+' ? $row->{cdsEnd} : $row->{cdsStart}+1 ),
	);

	my @exonStarts = split(',', $row->{exonStarts});
	my @exonEnds   = split(',', $row->{exonEnds});

	if ( defined $gene{_txStart} && defined $gene{_txEnd} ) {
		unless ( $exonStarts[0] == $row->{txStart} &&
			$exonEnds[-1]  == $row->{txEnd}   ) {
			print Dumper \%gene;
			croak "Exon starts and ends should be consistent with txStart/End!\n";
		}
	}
	else {
		if ( $gene{_strand} > 0 ) {
			$gene{_txStart} = $exonStarts[0] + 1;$gene{_txEnd} = $exonEnds[-1];
		}
		else {
			$gene{_txEnd} = $exonStarts[0] + 1; $gene{_txStart} = $exonEnds[-1];
		}
	}

	if ( defined $gene{_cdsStart} && defined $gene{_cdsEnd} ) {
		unless ( $exonStarts[0] <= $row->{cdsStart} &&
			$exonEnds[-1] >= $row->{cdsEnd} ) {
			print Dumper \%gene;
			print STDERR "Exon starts and ends should be consistent with cdsStart/End\n";
		}
	}

	my $gene_spans = Set::IntSpan->new([[ $gene{_strand} > 0 ?
		($gene{_txStart},$gene{_txEnd}) : ($gene{_txEnd},  $gene{_txStart}) ]]);

	# Some predicted gene may have zero length exons or introns,
  	# we discard them, and merge spuriously segmented ones
  	my $exon_spans = Set::IntSpan->new([ pairwise { [$a+1, $b] } @exonStarts, @exonEnds ]);
  	my $intron_spans = $gene_spans - $exon_spans;
  	my @exons = $exon_spans->spans();
  	my @introns = $intron_spans->spans();

	# Add exons and introns fields, and label them from 5' end to 3' end
	# labels are integers starting from 1.
	if ($gene{_strand} > 0) {
  		$gene{_exons} = Array::IntSpan->new( 
  			map { [$exons[$_][0],$exons[$_][1],$_+1] } 0..$#exons );
  		$gene{_introns} = Array::IntSpan->new( 
  			map { [$introns[$_][0],$introns[$_][1],$_+1] } 0..$#introns );
  	}
  	else {
  		$gene{_exons} = Array::IntSpan->new( 
  			map { [$exons[$_][0],$exons[$_][1],@exons-$_] } 0..$#exons );
  		$gene{_introns} = Array::IntSpan->new( 
  			map { [$introns[$_][0],$introns[$_][1],@introns-$_] } 0..$#introns );
  	}

	# Add cds and UTRs fields if they exist
	if ( $row->{cdsStart} != $row->{cdsEnd} ) {
 		if ($gene{_strand} > 0) {
  			$gene{_coding_exons} = $gene{_exons}->get_range($gene{_cdsStart}, $gene{_cdsEnd});
  			$gene{_five_utrs} = $gene{_exons}->get_range($gene{_txStart}, $gene{_cdsStart}-1);
  			$gene{_three_utrs} = $gene{_exons}->get_range($gene{_cdsEnd}+1,  $gene{_txEnd});
	  		# Adjust the labels to reflect the orders of sub-structures
	  			$_->[2] -= ($gene{_coding_exons}[0][2]-1) foreach ( @{$gene{_coding_exons}} );
	  			$_->[2] -= ($gene{_five_utrs}[0][2]-1) foreach ( @{$gene{_five_utrs}} );
	  			$_->[2] -= ($gene{_three_utrs}[0][2]-1) foreach ( @{$gene{_three_utrs}} );
		}
		else {
			$gene{_coding_exons} = $gene{_exons}->get_range($gene{_cdsEnd}, $gene{_cdsStart});
			$gene{_five_utrs} = $gene{_exons}->get_range($gene{_cdsStart}+1,  $gene{_txStart});
			$gene{_three_utrs} = $gene{_exons}->get_range($gene{_txEnd}, $gene{_cdsEnd}-1);
			$_->[2] -= ($gene{_coding_exons}[-1][2]-1) foreach ( @{$gene{_coding_exons}} );
			$_->[2] -= ($gene{_five_utrs}[-1][2]-1) foreach ( @{$gene{_five_utrs}} );
			$_->[2] -= ($gene{_three_utrs}[-1][2]-1) foreach ( @{$gene{_three_utrs}} );
		}
	} 
	else {
		$gene{_coding_exons} = Array::IntSpan->new();
		$gene{_five_utrs} = Array::IntSpan->new();
		$gene{_three_utrs} = Array::IntSpan->new();
	}

	return bless \%gene, $class;
}

=head2 Some Query Methods

$gene->is_coding() : Check if it is a protein coding gene: 0/1.

$gene->range()     : Return the genome interval of the gene: array ref

$gene->first_coding_exon() / $gene->last_coding_exon() :
    Return [st, ed, ind] of the first/last coding exon

$gene->codon_pos()
    Return a list of condon positions.

=cut

sub is_coding {  @{ $_[0]->coding_exons } > 0 }

sub range { $_[0]->strand > 0 ? [ $_[0]->txStart, $_[0]->txEnd ] :
			                    [ $_[0]->txEnd , $_[0]->txStart] ; }

sub first_coding_exon {
	my $gene = shift;
	$gene->get_exon( $gene->exons->lookup($gene->cdsStart) );
}

sub last_coding_exon {
 	my $gene = shift;
  	$gene->get_exon( $gene->exons->lookup($gene->cdsEnd) );
}

sub condon_pos {
	my $gene = shift;
	return () unless $gene->is_coding;
	my @coding_exons = @{ $gene->coding_exons() };
	my @positions;
  	foreach my $exon (@coding_exons) {
		push @positions, ($exon->[0]..$exon->[1]);
 	}
  	if ($gene->strand < 0) {
		@positions = reverse(@positions);
  	}
  	my @codon_positions;
  	for (my $ii = 0; $ii < int(@positions/3); $ii ++) {
		push @codon_positions, [@positions[(3*$ii)..(3*$ii+2)]];
  	}
  	return @codon_positions;
}

=head2 Get methods

All queriable fields: chrom, name, strand, txStart, txEnd, cdsStart, cdsEnd,
exons, introns, cds, five_utrs, three_utrs. If the field is represented as complex data
structure like exons, a deep clone is created and returned, to avoid tampering
the original object data.

For complex structure, also support fetch by index method. For example, get_exon(2),
will fetch the interval of the third exon (unblessed). The index are one-based, to
keep consistent with the label. The strandness will be regarded. So the third exon
is always the third one counted from the 5'end.

=cut

my $accessible = Set::Scalar->new( qw|_name _chrom _strand
									  _txStart _txEnd _cdsStart _cdsEnd
									  _exons _introns _coding_exons _five_utrs _three_utrs| );
my $complex    = Set::Scalar->new( qw|_exons _introns _coding_exons _five_utrs _three_utrs| );

our $AUTOLOAD;

sub AUTOLOAD {
  my ($self) = @_;
  my ($attribute) = ($AUTOLOAD =~ /::(\w+)$/);
  return if $attribute eq 'DESTROY'; # ignore DESTORY method

  if ($attribute =~ /^get(_\w+)$/) {
	my $field = $1; $field .= 's'; # to plural form
	croak "Cannot access element from complex field $field" unless $complex->has($field);
	{
	  no strict 'refs';
	  *{$AUTOLOAD} = sub {
		my ($self, $index) = @_;
		my $num = @{$self->{$field}};
		if ( $self->strand > 0) { # AUTOLOAD is re-entrant
		  return $index >= 1 && $index <= $num ?
			dclone($self->{$field}[$index-1])    : undef;
		}
		else {
		  return $index >= 1 && $index <= $num ?
			dclone($self->{$field}[$num-$index]) : undef;
		}
	  };
	}
  }
  elsif ($attribute =~ /^num(_\w+)$/) {
	my $field = $1;
	croak "Cannot access complex field $field" unless $complex->has($field);
	{
	  no strict 'refs';
	  *{$AUTOLOAD} = sub {
		my ($self, $index) = @_;
		my $num = @{$self->{$field}};
		return $num;
	  };
	}
  }
  else {
	$attribute = '_'.$attribute;  # to private name
	croak "Cannot access $attribute attribute" unless $accessible->has($attribute);
	{
	  no strict 'refs';
	  *{$AUTOLOAD} = sub {
		my ($self) = @_;
		return ref $self->{$attribute} ? dclone($self->{$attribute}) : $self->{$attribute};
	  };
	}
  }
  goto &{$AUTOLOAD};
}


=head2 mRNA and peptide sequence

my $cds = $gene->cds($tbo);

my $peptide_sequence = $gene->peptide;

=cut

sub _seq_for_blocks {
  my ($chrom, $blocks, $tbo) = @_;
  my $seq;
  foreach my $bk (@{$blocks}) {
	$seq .= $tbo->fetch(range_to_spec($chrom, $bk->[0], $bk->[1]));
  }
  return $seq;
}

sub cds {
  my ($self, $tbo, $quite) = @_;
  my $cds = _seq_for_blocks($self->chrom, $self->coding_exons, $tbo);
  unless ($quite) {
	warn "CDS DNA length is not a multiple of three?" if length($cds) % 3;
  }
  return $self->strand > 0 ? $cds : rev_comp($cds);
}

sub peptide {
  my $self = shift;
  return dna2peptide($self->cds(@_));
}

=head2 CDS/Transcript-to-Genome Mapping

Construct subroutines to map relative positions on cds dna or mRNA
to the genome backbone, or vice versa. It can be used to map nonsyn-SNPs
to genomic positions.

To provide more flexible construction of aligned block, we can provide an
optional callback, that accept $self as the first argument, to produce a list
of blocks.

Return a pair of subroutine reference that can be called to convert betweem
genomic position and relative position on the cds.

=cut

sub align {
  my ($self, $cb) = @_;
  # or you can call Genome::Annotation::Transcript::exons to get transcript
  # to genome mapping
  $cb = \&Genome::Annotation::Transcript::coding_exons unless defined $cb;
  my $gene_blocks = $cb->($self);

  my (@ali_blocks, @local_to_genome, @genome_to_local);
  if ( $self->strand > 0 ) {
	my $pstart = 1;
	foreach my $bk (@$gene_blocks) {
	  # Query is the local sequence, and should always be on the forward strand
	  my $val =  { qName   => $self->name,
				   qStart  => $pstart,
				   qEnd    => $pstart+$bk->[1]-$bk->[0],
				   tName   => $self->chrom,
				   tStart  => $bk->[0],
				   tEnd    => $bk->[1],
				   tStrand => 1,
				 };
	  push @local_to_genome => [ $val->{qStart},$val->{qEnd},$val ];
	  push @genome_to_local => [ $val->{tStart},$val->{tEnd},$val ];
	  $pstart = $val->{qEnd}+1;
	}
  }
  else {
	my $length = 0;
	foreach my $bk (@$gene_blocks) {
	  $length += $bk->[1]-$bk->[0]+1;
	}
	foreach my $bk (@$gene_blocks) {
	  my $val =  { qName   => $self->name,
				   qStart  => $length-($bk->[1]-$bk->[0]),
				   qEnd    => $length,
				   tName   => $self->chrom,
				   tStart  => $bk->[1],
				   tEnd    => $bk->[0],
				   tStrand => -1,
				 };
	  unshift @local_to_genome => [ $val->{qStart},$val->{qEnd},$val ];
	  push @genome_to_local => [ $val->{tEnd},$val->{tStart},$val ];
	  $length = $val->{qStart}-1;
	}
  }
  my $local_to_genome = Array::IntSpan->new(@local_to_genome);
  my $genome_to_local = Array::IntSpan->new(@genome_to_local);

  # A pair of subroutines to carry out the mapping process
  my $local_to_genome_mapper = sub {
	my ($local_pos) = @_;
	if (my $val = $local_to_genome->lookup($local_pos)) {
	  return ($local_pos-$val->{qStart})*$val->{tStrand}+$val->{tStart};
	}
	return undef;
  };

  my $genome_to_local_mapper = sub {
	my ($genome_pos) = @_;
	if (my $val = $genome_to_local->lookup($genome_pos)) {
	  return ($genome_pos-$val->{tStart})*$val->{tStrand}+$val->{qStart};
	}
	return undef;
  };

  return ($local_to_genome_mapper, $genome_to_local_mapper);
}


=head1 AUTHOR

Xueya Zhou, C<< <xueyazhou at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-genome at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Genome>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Genome::UCSC::Transcript


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

1; # End of Genome::UCSC::Transcript

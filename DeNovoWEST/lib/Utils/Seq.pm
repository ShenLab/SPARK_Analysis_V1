package Utils::Seq;

use strict;
use warnings;
use Carp;
use IO::File;
use IO::Detect;
use Iterator::Simple;
use List::MoreUtils qw(mesh all);
use File::Temp qw(tempfile);

our @EXPORT_OK = qw(read_fasta iter_fasta extract_fasta print_seq print_seqs
					rev_comp toRNA fa_count is_snv is_sym is_asym is_ti is_tv
					is_aa is_dna is_rna	iub2regex is_fourD  codon2AA iub3to1 iub1to3
					dna2peptide translate_DNA
					);

use base qw(Exporter);

=head1 NAME

Utils::Seq - Simple biological sequence manipulation tool!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Quick summary of what the module does.
    ...

=head1 EXPORT

readFasta iterFasta printSeq printSeqs revComp toRNA faCount 
iub2Regex isFourD  codon2AA iub3to1 iub1to3 dna2Peptide translateDNA

=head1 SUBROUTINES/METHODS

=head2 read_fasta FILE

  Extract all sequences from fasta format.

=cut

sub read_fasta {
  my($file) = @_;
  my $fh;
  if (is_filehandle($file)) {
  	$fh = $file;
  }
  else {
  	$fh = IO::File->new($file) or croak "Cannot open FASTA file $file for reading";
  }

  # Declare and initialize variables
  my %sequences = ();
  my $seqname;
  while ( my $line =<$fh> ) {
	# discard blank line
	if ($line =~ /^\s*$/) {
	  next;
	  # discard comment line
	} elsif($line =~ /^\s*#/) {
	  next;
	  # discard fasta header line
	} elsif($line =~ /^>(.*)/) {
	  $seqname = $1;
	  next;
	  # keep line, add to sequence string
	} else {
	  $sequences{$seqname} .= $line;
	}
  }
  # remove non-sequence data (in this case, whitespace) from $sequence string
  foreach my $id (keys %sequences) {
	$sequences{$id} =~ s/\s//g;
  }
  return \%sequences;
}


=head2 iter_fasta FILE

  Provide an iterator to sequences in fasta file. Each call return an arrayref
  [SeqName, Seq].

  We assume that each sequence may be separated in different lines, but no
  blanks lines can be found within sequence blocks.

=cut

sub iter_fasta {
  my($file) = @_;
  my $fh;
  if (is_filehandle($file)) {
  	$fh = $file;
  }
  else {
  	$fh = IO::File->new($file) or croak "Cannot open FASTA file $file for reading";
  }

  my ($src);
  while(<$fh>) {
	if (/^>(.*)/) {
	  $src = $1;
	  last;
	}
  }
  return iterator {
	   my $seq;
	   while(<$fh>) {
		 if(/^>(.*)/) {
		   my @out = ($src, $seq);
		   $src = $1;
		   $seq = "";
		   return \@out;
		 }
		 else {
		   next if /^\s*$/;
		   chomp;
		   $seq .= $_;
		 }
	   }
	   if ($seq) {
		 my @out = ($src, $seq);
		 $seq = "";
		 return \@out;
	   }
	   return;
	 };
}


=head2 extract_fasta FILE, NAME

   Extract sequence of given name from fasta file

=cut

sub extract_fasta {
    my ($file, $seqname) = @_;
    open my $fin, $file or die "Cannot open $file";
    my $seq;
    while(<$fin>) {
        last if /^>$seqname\s+/;
    }
    while(<$fin>) {
        last if /^>/;
        chomp;
        $seq .= $_;
    }
    close $fin;
    return $seq;
}


=head2 rev_comp DNA

   Reverse complement the DNA sequence : pass through non-agtc chars.
   It can deal with IUPAC code appropriately.

=cut

my %rc = ('A'=>'T', 'T'=>'A', 'G'=>'C', 'C'=>'G',
		  'a'=>'t', 't'=>'a', 'g'=>'c', 'c'=>'g',
		  'W'=>'S', 'S'=>'W', 'w'=>'s', 's'=>'w',
		  'Y'=>'R', 'R'=>'Y', 'r'=>'y', 'y'=>'r',
		  'M'=>'K', 'K'=>'M', 'm'=>'k', 'k'=>'m',
	    'B'=>'V', 'V'=>'B', 'b'=>'v', 'v'=>'b',
		  'D'=>'H', 'H'=>'D', 'h'=>'d', 'd'=>'h',
      );

sub rev_comp {
  my ($seq) = @_;
  my $rcSeq = reverse $seq;
  for (my $i = 0;  $i < length($rcSeq);  $i++) {
    my $base = substr($rcSeq, $i, 1);
    my $cBase = $rc{$base} || $base;
    substr($rcSeq, $i, 1, $cBase);
  }
  return $rcSeq;
}

=head2 toRNA DNA

   Substitute T to U.

=cut

sub toRNA {
  my ($seq) = @_;
  $seq =~ s/T/U/g;
  return $seq;
}


=head2 fa_count SEQ

  Count the occurance of each base type in the sequence.

=cut

sub fa_count {
  my ($seq) = @_;
  my %count = (A => 0, G => 0, C => 0, T => 0, CG => 0);
  for (my $ii = 0; $ii < length($seq); ++$ii) {
	$count{uc(substr($seq, $ii, 1))} ++;
	if ($ii < length($seq) - 1 && uc(substr($seq, $ii, 2)) eq 'CG') {
	  $count{CG} ++;
	}
  }
  my $totbp = $count{G}+$count{C}+$count{A}+$count{T};
  if ($totbp > 0) {
	$count{GCPerc}  = ($count{G}+$count{C})/$totbp;
	$count{CpGPerc} = $count{CG}/$totbp;
  }
  else {
	$count{GCPerc} = undef;
	$count{CpGPerc} = undef;
  }
  return \%count;
}

=head2 print_seq SEQ, FH [,LEN]

  Format and print sequence data. Default output length per line is 50.

=cut

sub print_seq {
  my ($sequence, $fh, $length) = @_;
  $fh = \*STDOUT unless defined $fh;
  $length = 50 unless defined $length;
  # Print sequence in lines of $length
  for ( my $pos = 0 ; $pos < length($sequence) ; $pos += $length ) {
	print $fh substr($sequence, $pos, $length), "\n";
  }
}

=head2 printSeqs SEQREF, FH [,LEN]

  Format and print a list of sequences.

=cut

sub print_seqs {
  my ($seq_ref, $fh, $length) = @_;
  foreach my $id (keys %$seq_ref) {
	print ">$id\n";
	printSeq($seq_ref->{$id}, $fh, $length);
  }
}

# IUPAC code for ambiguous DNA sequence
my %iub2charset =
  (
   A => 'A',
   C => 'C',
   G => 'G',
   T => 'T',
   R => '[GA]',
   Y => '[CT]',
   M => '[AC]',
   K => '[GT]',
   S => '[GC]',
   W => '[AT]',
   B => '[CGT]',
   D => '[AGT]',
   H => '[ACT]',
   V => '[ACG]',
   N => '[ACGT]',
  );


=head1 is_snv

Test if two different alleles are both ACGT.

=cut

sub is_snv {
  my @a = @_;
  if ($a[0] =~ /^[ACGT]$/ && $a[1] =~ /^[ACGT]$/) {
    return 1;  
  }
  else {
    return 0;
  }
}

=head2 is_sym/is_asym ALLELE1, ALLELE2

Test if two alleles are symmetric or not.

=cut

sub is_sym {
  my (@a) = @_;
  if ($a[0] eq 'C' && $a[1] eq 'G' ||
    $a[0] eq 'G' && $a[1] eq 'C' ||
    $a[0] eq 'A' && $a[1] eq 'T' ||
    $a[0] eq 'T' && $a[1] eq 'A'
   ) {
  return 1;
  }
  else {
  return 0;
  }
}

sub is_asym {
  my (@a) = @_;
  if ($a[0] eq 'C' && $a[1] eq 'A' ||
    $a[0] eq 'C' && $a[1] eq 'T' ||
    $a[0] eq 'G' && $a[1] eq 'A' ||
    $a[0] eq 'G' && $a[1] eq 'T' ||
    $a[0] eq 'A' && $a[1] eq 'G' ||
    $a[0] eq 'A' && $a[1] eq 'C' ||
    $a[0] eq 'T' && $a[1] eq 'G' ||
    $a[0] eq 'T' && $a[1] eq 'C'
   ) {
  return 1;
  }
  else {
  return 0;
  }
}


=head2 is_ti/is_tv ALLELE1 ALLELE2

Test of two alleles are transition or transversion

=cut

sub is_ti {
  my (@a) = @_;
  if ($a[0] eq 'A' && $a[1] eq 'G' ||
    $a[0] eq 'G' && $a[1] eq 'A' ||
    $a[0] eq 'C' && $a[1] eq 'T' ||
    $a[0] eq 'T' && $a[1] eq 'C'
   ) {
  return 1;
  }
  else {
  return 0;
  }
}

sub is_tv {
  my (@a) = @_;
  if ($a[0] eq 'A' && $a[1] eq 'C' ||
    $a[0] eq 'A' && $a[1] eq 'T' ||
    $a[0] eq 'C' && $a[1] eq 'A' ||
    $a[0] eq 'C' && $a[1] eq 'G' ||
    $a[0] eq 'G' && $a[1] eq 'C' ||
    $a[0] eq 'G' && $a[1] eq 'T' ||
    $a[0] eq 'T' && $a[1] eq 'A' ||
    $a[0] eq 'T' && $a[1] eq 'G'
   ) {
  return 1;
  }
  else {
  return 0;
  }
}


=head2 is_dna/is_rna/is_aa

Test if the sequence is a valid DNA/RNA/AA sequence.

=cut

sub is_dna {
	my $DNAseq = shift @_;
	if ($DNAseq =~ /^[ACGT]+$/i) {
		return 1;
	}
	return 0;
}

sub is_rna {
	my $RNAseq = shift @_;
	if ($RNAseq =~ /^[ACGU]+$/i) {
		return 1;
	}
	return 0;
}

sub is_aa  {
  my $AAseq = shift @_;
  if ($AAseq =~ /^[ACDEFGHIKLMNPQRSTVWY]+$/) {
	return 1;
  }
  return 0;
}


=head2 iub2Regex IUB_SEQ

  Return the regex for IUPAC sequence.

  e.g. $re = qr/iub2Regex("TRAGN")/;

=cut

sub iub2regex {
  my($iub) = @_;
  my $regexp = '';
  # Translate each character in the iub sequence
  for ( my $i = 0 ; $i < length($iub) ; ++$i ) {
	$regexp .= $iub2charset{substr($iub, $i, 1)};
    }
  return $regexp;
}

## Below are subroutines for DNA/protein transformation.

# The universal genetic code table
# If you want to deal with special case like some virus, resort to bioperl!
# q(_) represent the stop codon.
my %genetic_code =
  (
   'TCA' => 'S',				# Serine
   'TCC' => 'S',				# Serine
   'TCG' => 'S',				# Serine
   'TCT' => 'S',				# Serine
   'TTC' => 'F',				# Phenylalanine
   'TTT' => 'F',				# Phenylalanine
   'TTA' => 'L',				# Leucine
   'TTG' => 'L',				# Leucine
   'TAC' => 'Y',				# Tyrosine
   'TAT' => 'Y',				# Tyrosine
   'TAA' => '_',				# Stop
   'TAG' => '_',				# Stop
   'TGC' => 'C',				# Cysteine
   'TGT' => 'C',				# Cysteine
   'TGA' => '_',				# Stop
   'TGG' => 'W',				# Tryptophan
   'CTA' => 'L',				# Leucine
   'CTC' => 'L',				# Leucine
   'CTG' => 'L',				# Leucine
   'CTT' => 'L',				# Leucine
   'CCA' => 'P',				# Proline
   'CCC' => 'P',				# Proline
   'CCG' => 'P',				# Proline
   'CCT' => 'P',				# Proline
   'CAC' => 'H',				# Histidine
   'CAT' => 'H',				# Histidine
   'CAA' => 'Q',				# Glutamine
   'CAG' => 'Q',				# Glutamine
   'CGA' => 'R',				# Arginine
   'CGC' => 'R',				# Arginine
   'CGG' => 'R',				# Arginine
   'CGT' => 'R',				# Arginine
   'ATA' => 'I',				# Isoleucine
   'ATC' => 'I',				# Isoleucine
   'ATT' => 'I',				# Isoleucine
   'ATG' => 'M',				# Methionine
   'ACA' => 'T',				# Threonine
   'ACC' => 'T',				# Threonine
   'ACG' => 'T',				# Threonine
   'ACT' => 'T',				# Threonine
   'AAC' => 'N',				# Asparagine
   'AAT' => 'N',				# Asparagine
   'AAA' => 'K',				# Lysine
   'AAG' => 'K',				# Lysine
   'AGC' => 'S',				# Serine
   'AGT' => 'S',				# Serine
   'AGA' => 'R',				# Arginine
   'AGG' => 'R',				# Arginine
   'GTA' => 'V',				# Valine
   'GTC' => 'V',				# Valine
   'GTG' => 'V',				# Valine
   'GTT' => 'V',				# Valine
   'GCA' => 'A',				# Alanine
   'GCC' => 'A',				# Alanine
   'GCG' => 'A',				# Alanine
   'GCT' => 'A',				# Alanine
   'GAC' => 'D',				# Aspartic Acid
   'GAT' => 'D',				# Aspartic Acid
   'GAA' => 'E',				# Glutamic Acid
   'GAG' => 'E',				# Glutamic Acid
   'GGA' => 'G',				# Glycine
   'GGC' => 'G',				# Glycine
   'GGG' => 'G',				# Glycine
   'GGT' => 'G',				# Glycine
  );

# Four fold degenerate codons

my %fourD_codons;
{
  #my ($i, $j);
  #foreach my $cp ( comp2 { $i . $j }
  #				   i => [qw(A C G T)],
  #				   j => [qw(A C G T)] ) {
  #	my $AA = $genetic_code{$cp.'A'};
  #	if ( all {$genetic_code{$cp.$_} eq $AA} qw(A C G T) ) {
  #	  %fourD_codons = ( %fourD_codons,
  #						( map {($cp.$_, $genetic_code{$cp.$_})} qw|A C G T| )
  #					  );
  #	}
  #}
  my @first2 = ();
  for my $i (qw|A C G T|) {
	for my $j (qw|A C G T|) {
	  push @first2, $i.$j;
	}
  }
  foreach my $cp (@first2) {
	my $AA = $genetic_code{$cp.'A'};
	if ( all {$genetic_code{$cp.$_} eq $AA} qw(A C G T) ) {
	  %fourD_codons = ( %fourD_codons,
						( map {($cp.$_, $genetic_code{$cp.$_})} qw|A C G T| )
					  );
	}
  }
}

# String of 3-letter IUB amino acid code to 1-letter code
# Note: Adding selenocysteine
my %three2one
  = (
	 'ALA' => 'A',
	 'VAL' => 'V',
	 'LEU' => 'L',
	 'ILE' => 'I',
	 'PRO' => 'P',
	 'TRP' => 'W',
	 'PHE' => 'F',
	 'MET' => 'M',
	 'GLY' => 'G',
	 'SER' => 'S',
	 'THR' => 'T',
	 'TYR' => 'Y',
	 'CYS' => 'C',
	 'ASN' => 'N',
	 'GLN' => 'Q',
	 'LYS' => 'K',
	 'ARG' => 'R',
	 'HIS' => 'H',
	 'ASP' => 'D',
	 'GLU' => 'E',
   'SEC' => 'U'
    );

my %one2three = mesh @{[ values %three2one ]}, @{[ keys %three2one ]};

=head2 DNA/Protein Code Transformation

  codon2AA : transform three letter DNA condon to one letter AA code
  isFourD   : check if the condon is four fold degenerate.
  iub3to1   : transform 3-letter AA seq into 1-letter AA seq
  iub1to3   : the vice versa

=cut


sub codon2AA {
  my ($codon) = @_;
  $codon = uc $codon;
  if ( defined $genetic_code{$codon} ) {
	return $genetic_code{$codon};
  }
  else{
	return 'X';
  }
}

sub is_fourD {
  my ($codon) = @_;
  $codon = uc $codon;
  if ( defined $fourD_codons{$codon} ) {
	return $fourD_codons{$codon};
  }
  else {
	return undef;
  }
}

sub iub3to1 {
  my($input, $strict) = @_;
  $input =~ s/\n/ /g;
  my $seq = '';
  # This use of split separates on any contiguous whitespace
  my @code3 = split(' ', $input);
  foreach my $code (@code3) {
	   # A little error checking
	   if(not defined $three2one{$code}) {
        if ($strict) {
            die "Cannot find 1-letter code for $code";
        }
        else {
            $seq .= 'X';
        }
	   }
     else {
        $seq .= $three2one{$code};
     }
  }
  return $seq;
}

sub iub1to3 {
  my($input, $strict) = @_;
  $input =~ s/\s/ /g;
  my $seq = '';
  for (my $ii = 0; $ii < length($input); $ii ++) {
	   my $code = substr($input, $ii, 1);
	   # A little error checking
	   if(not defined $one2three{$code}) {
        if ($strict) {
            die "Cannot find 3-letter code for $code";
        }
        else {
            $seq .= ($ii == 0 ? '':' ') . 'XXX'; 
        }
	   }
     else {
         $seq .= ($ii == 0 ? '':' ') . $one2three{$code};
     } 
  }
  return $seq;
}

=head2 DNA/protein translation

  $pep = dna2Peptide($dna);

  $prot = translateDNA($dna, $start, $end);

  where $start and $end are one-based postions.

=cut

sub dna2peptide {
  my($dna) = @_;
  my $protein = '';
  # Translate each three-base codon to an amino acid, and append to a protein 
  for(my $i=0; $i < (length($dna) - 2) ; $i += 3) {
	$protein .= codon2AA( substr($dna,$i,3) );
  }
  return $protein;
}

sub translate_DNA {
  my($seq, $start, $end) = @_;
  my $protein;
  # To make the subroutine easier to use, you won't need to specify
  #  the end point-it will just go to the end of the sequence
  #  by default.
  $start = 1 unless defined $start;
  $end = length($seq) unless defined $end;

  # Finally, calculate and return the translation
  return dna2peptide ( substr ( $seq, $start - 1, $end -$start + 1) );
}



=head1 AUTHOR

Xueya Zhou, C<< <xueyazhou at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-utils at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Utils>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Utils::Seq


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

1; # End of Utils::Seq

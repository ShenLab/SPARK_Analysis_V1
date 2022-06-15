package Genome::UCSC::TwoBit;

use strict;
use warnings;
use Carp;
use IO::File;
use Fcntl qw(:seek);
use Set::IntSpan;
use Data::Dumper;
use Genome::Ranges;
use Utils::Hash qw|merge_opts|;

use constant            ONE_BYTE => 1;
use constant           FOUR_BYTE => 4;
use constant       BITS_PER_BYTE => 8;
use constant BASES_PER_FOUR_BYTE => 16;

# DNA Convertion Table
my %BITMAP = ('00' => 'T', '01' => 'C', '10' => 'A', '11' => 'G');
my @BYTEMAP = map{
  join '', map $BITMAP{ $_ }, unpack '(A2)4', unpack 'B8', chr
} 0 .. 255;


=head1 NAME

Genome::UCSC::TwoBit - Read genome sequence from UCSC twobit file.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

	use Genome::UCSC::TwoBit;
	my $tb = Genome->new("/gbdb/hg19/hg19.2bit");
	my $fin = IO::File->new("Range_file.bed");
	while(<$fin>) {
	    my ($chr, $st, $ed) = split;
	    # NOTE: genomic positions are half-open, zero-based
	    my $seq = $tb->fetch("$chr:$st-$ed", {noMask => 1})
	    ...
	}


=head1 SUBROUTINES/METHODS

=head2 $self->new FILE

Construct the object for 2bit file.

=cut

sub new {
	my ($class, $file) = @_;
	open(my $fh, '<', $file) or croak "Unable to open '$file' for reading: $!";
 	# Read header
 	my $raw = '';
 	sysread($fh, $raw, FOUR_BYTE * 4);
  	# Parse header
  	my ($sig, $ver, $count, $reserved) = unpack('l4', $raw);
  	my $self = {TBF => $fh, SIG => $sig, VER => $ver, CNT => $count, RSV => $reserved};

  	# Validate (signature, reverse byte order, version)
  	# Currently only support version 0
  	unless($reserved eq '0' && $ver eq '0' && $sig eq hex("0x1A412743") ) {
  		print Dumper $self;
  		croak "Incorrect header";
  	}
  	# Find offsets
  	my (%offsets);
  	my ($size, $cnt, $name) = ('', '', '');
  	for (1..$count) {
		# Read size of record name
		sysread($fh, $raw, ONE_BYTE);
		$size = unpack('C', $raw);
		# Read name of reacord
		sysread($fh, $name, $size);
		# Read and store offset for each seq
		sysread($fh, $raw, FOUR_BYTE);
		$offsets{$name} = unpack('l', $raw);
	}
 	# Then fetch sequence attributes (sequence lengths, N and M blocks) for faster retrieval
  	# And set new offsets for where sequence blocks actually begin
  	my (%newoffsets, %seqlens);
  	foreach my $name (keys %offsets) {
  		my (@nblocks, @mblocks);
  		# print STDERR $name, "\n";
  		sysseek($fh, $offsets{$name}, SEEK_SET);
  		sysread($fh, $raw, FOUR_BYTE);
  		$seqlens{$name} = unpack('l', $raw);
		# N/M blocks
		foreach my $blocks (\@nblocks, \@mblocks) {
	  		# fetch the block count
		  	sysread($fh, $raw, FOUR_BYTE);
		  	$cnt = unpack('l', $raw);
		  	if ($cnt) {
	  			sysread($fh, $raw, FOUR_BYTE * $cnt);
	  			my @starts = unpack("l$cnt", $raw);
	  			sysread($fh, $raw, FOUR_BYTE * $cnt);
	  			my @lens = unpack("l$cnt", $raw);
	 	 		@$blocks = map { [ $starts[$_], $starts[$_]+$lens[$_]-1 ] } 0..$#starts;
	  		}
	  	}
		# Throw away reserved field
		sysread($fh, $raw, FOUR_BYTE);
		# Determine the new offset
		$newoffsets{$name} = sysseek($fh, 0, SEEK_CUR);
		# Define block span, position are inclusive
		$self->{NBLOCK}{$name} = Set::IntSpan->new(\@nblocks);
		$self->{MBLOCK}{$name} = Set::IntSpan->new(\@mblocks);
	}
	$self->{OFFSET} = \%offsets;
	$self->{NEWOFFSET} = \%newoffsets;
	$self->{SEQLEN} = \%seqlens;
	return bless $self, $class;
}


=head2 $self->exists CHROM [, POS]

Test if a chromosome exist and if so test if the position is within total length.
Do this test before calling other functions can avoid croak.

=cut

sub exists {
	my ($self, $chrom, $start, $end) = @_;
    if (defined $start) {
       $end = $start unless defined $end;
       if (exists  $self->{SEQLEN}{$chrom}) {
            if ($start > 0 && $end >= $start && $end <= $self->{SEQLEN}{$chrom}) {
                return 1;
            }
            else {
                return 0;
            }
        }
        else {
        	return 0;
        }
    }
    else {
        if (exists  $self->{SEQLEN}{$chrom}) { 
            return 1;
        }
        else {
            return 0;
        }
    }
}


=head2 $self->size CHROM [, OPTIONS]

Return the length of a chromosome.

noNs: default is 0, if switched on, the length does not count Ns.

=cut

sub size {
	my ($self, $chrom, $argref) = @_;
	my $arg = merge_opts($argref, noNs => 0);

	croak "Chromosome $chrom cannot be found" unless $self->{SEQLEN}{$chrom};
	if ($arg->{noNs}) {
		return $self->{SEQLEN}{$chrom};
	}
	else {
		my $nlen;
		if (defined $self->{NBLOCK}{$chrom}) {
			$nlen = $self->{NBLOCK}{$chrom}->size();
		}
		else {
			$nlen = 0;
		}
		return $self->{SEQLEN}{$chrom} - $nlen;
	}

}


=head2 $self->gaps

Return all gaps (regions with N's) found in the genome sequence.
The return value is a L<Genome::Ranges> object.

=cut

sub gaps {
	my ($self) = @_;
	my %nbed;
	foreach my $chrom (keys %{$self->{NBLOCK}}) {
		push @{$nbed{$chrom}}, $self->{NBLOCK}->spans();
	}
	return Genome::Ranges->new(\%nbed);
}

=head2 $self->get_slice CHROM, START, END [, OPTIONS]

Fetch sequence given a genomic range.

B<NOTE> Internally genomic positions are half-open, zero-based; but externally
we use the inhouse representation of genomic intervals that is 1-based.

=head3 Options

=over 5

=item * noMask

Convert sequence to all upper case; otherwise, repeat-masked sequence will appear in lower cases.
The default is 1.

=back

=cut

sub get_slice {
	my ($self, $chrom, $chromStart, $chromEnd, $argref) = @_;
	my $arg = merge_opts($argref, noMask => 1);

	$chromStart -= 1;
	$chromEnd = $chromStart + 1 unless defined $chromEnd;

	croak "Incorrect genomic range specification"
		unless defined $chrom && defined $chromStart && $chromStart >= 0
			&& defined $chromEnd && $chromEnd > $chromStart;

	if (!defined $self->{SEQLEN}{$chrom}) {
		croak "Cannot find chromosome $chrom";
	}
	if ($chromEnd > $self->{SEQLEN}{$chrom}) {
		croak "Interval beyond chromosome boundaries: $chrom:$chromStart-$chromEnd";
	}

  	# Seek to the seq blocks
  	sysseek($self->{TBF}, $self->{NEWOFFSET}{$chrom}, SEEK_SET);

  	my $packedStart = ($chromStart >> 2);
  	my $packedEnd = (($chromEnd + 3) >> 2);
  	my $packByteCount = $packedEnd - $packedStart;
  	sysseek($self->{TBF}, $packedStart, SEEK_CUR);

  	my ($raw, $dna, $size, $cnt) = ('', '', '', '');
  	sysread($self->{TBF}, $raw, $packByteCount);

  	#my $bytes = ((int( ($chromEnd-$chromStart+1) / BASES_PER_FOUR_BYTE)) + 1) * FOUR_BYTE;
  	#$dna = join '', map $BITMAP{$_},
  	# unpack('(A2)*', unpack("B" . $bytes * BITS_PER_BYTE , $raw));
  	$dna = join("", map { $BYTEMAP[$_] } unpack 'C*', $raw);

  	# Trim out header and tail
  	if ($chromStart % 4) {
		my $trimStart = $chromStart % 4;
		substr($dna, 0, $trimStart, '');
  	}
  	if ($chromEnd % 4) {
		my $trimEnd = 4-($chromEnd % 4);
		substr($dna, length($dna)-$trimEnd, $trimEnd, '');
  	}

  	# Find overlapping N-blocks, using run-list specification
  	my $query_rng = Set::IntSpan->new([[$chromStart, $chromEnd-1]]);
  	my $noverlap = $query_rng*$self->{NBLOCK}{$chrom};

  	if (!$noverlap->empty) {
  		print "Overlap with N\n";
  		my @spans = $noverlap->spans();
  		foreach my $span (@spans) {
  			my ($st, $ed) = @$span;
  			my $len = $ed - $st + 1;
  			my $st_rel = $st - $chromStart;
  			substr($dna, $st_rel, $len, "N"x$len);
  		}
  	}

  	unless ($arg->{noMask}) {
  		my $moverlap = $query_rng*$self->{MBLOCK}{$chrom};
  		if (!$moverlap->empty) {
  			print "Overlap with Mask\n";
  			my @spans = $moverlap->spans();
  			foreach my $span (@spans) {
  				my ($st, $ed) = @$span;
  				my $len = $ed - $st + 1;
  				my $st_rel = $st - $chromStart;
  				substr($dna, $st_rel, $len, lc(substr($dna, $st_rel, $len)));
  			}
  		}
  	}

  	return $dna;
}

=head2 $self->get_base CHROM, POS

The same as above, but only fetch sequence at one bp.

=cut

sub get_base {
	my ($self, $chrom, $pos) = @_;
	return $self->get_slice($chrom, $pos, $pos);
}

=head1 AUTHOR

Xueya Zhou, C<< <xueyazhou at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-genome at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Genome>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Genome::UCSC::TwoBit


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

1; # End of Genome::UCSC::TwoBit

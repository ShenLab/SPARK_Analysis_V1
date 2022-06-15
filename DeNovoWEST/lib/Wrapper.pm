package Wrapper;


use strict;
use warnings;
use Carp;
use List::Util qw|sum first|;
use List::MoreUtils qw|uniq all any|;
use Perl6::Slurp;
use File::Temp qw|tempdir|;
use File::Path qw|make_path|;
use File::Basename qw|basename|;
use Utils::Hash qw|merge_opts|;
use Utils::File qw|swap_files open_file|;
use Utils::File::Iter qw|iter_file|;
use Genet::File::VCF qw|slurp_vcf_header|;

use base qw|Exporter|;


our @EXPORT_OK = qw|bam_sampids bam_rgids bam_ismapped vcf_sampids vcf_rename
					plink_run file_exists extract_bam iter_sigdat iter_normdp|;

=head1 SYNOPSIS

Utility functions that wrap the output from some some executables.

=head2 bam_sampids/bam_rgids 

Extract sample ID or read group ID from BAM/CRAM file

=cut

sub bam_sampids {
 	my ($bamfile) = @_;
 	#croak "Cannot find BAM file: $bamfile" unless -f $bamfile;
	my @sampid;
	open my $fin, "htsfile -h $bamfile | " or die "Cannot open pipe";
	while(<$fin>) {
		if (/^\@RG/) {
			if (/SM:(\S+)/) {
				push @sampid, $1;
			}
			else {
				croak "Cannot find read group ID: $_";
			}
		}
	}
	my @uqid = uniq sort @sampid;
	return @uqid;
}

sub bam_rgids {
	my ($bamfile) = @_;
	#croak "Cannot find BAM file: $bamfile" unless -f $bamfile;
	my @rgid;
	open my $fin, "htsfile -h $bamfile | " or die "Cannot open pipe";
	while(<$fin>) {
		if (/^\@RG/) {
			if (/ID:(\S+)/) {
				push @rgid, $1;
			}
			else {
				croak "Cannot find read group ID: $_";
			}
		}
	}
	return @rgid;
}

=head2 bam_ismapped 

Test if a bam file is mapped.
Currently, it only looks for sequence dictionary in the header.
Return the first sequence name.

=cut

sub bam_ismapped {
	my ($bamfile) = @_;
	#croak "Cannot find BAM file: $bamfile" unless -f $bamfile;
	my $flag;
	open my $fin, "htsfile -h $bamfile | " or die "Cannot open pipe";
	while(<$fin>) {
		if (/^\@SQ/) {
			my ($chrom) = ($_ =~ /SN:(\w+)/);
			unless(defined $chrom) {
				die "Cannot find sequence name in header line begin with \@SQ";
			}
			$flag = $chrom;
			last;
		}
	}
	return $flag;
}


=head2 vcf_sampids

Sample IDs in VCF file.
Wrapper of vcfsamplenames from vcflib if bgzip'ed and index.
Otherwise, it will parse VCF using inhouse lib.

Note: we have modified function to use bcftools

=cut 

sub vcf_sampids {
	my ($vcf) = @_;
	croak "Cannot find VCF file $vcf" unless -f $vcf;
	my @iids;
	open my $fin, "bcftools query -l $vcf |" or die "Cannot open bcftools pipe";
	#if ($vcf =~ /gz$/ && -f "$vcf.tbi") {
	#	open my $fin, "vcfsamplenames $vcf |" or die "Cannot open pipe";
		while(<$fin>) {
			push @iids, (split)[0];
		}
		my @uqid = uniq sort @iids;
		croak "Duplicated sample names exist in VCF" unless @uqid == @iids;
	#}
	#else {
	#	my ($header, $sampids) = slurp_vcf_header($vcf);
	#	@iids = @$sampids;
	#}
	return @iids;
}

=head2 vcf_rename 

Rename sample IDs in VCF file in place.
Wrapper of bcftools reheader function.

=cut

sub vcf_rename {
	my ($vcf, $rename, $wrkdir) = @_;
	croak "Cannot find VCF file $vcf" unless -f $vcf;
	my @origids = vcf_sampids($vcf);
	my @rename;
	while(my ($oldiid, $newiid) = each %$rename) {
		next unless grep { $oldiid eq $_ } @origids;
		next if $oldiid eq $newiid;
		push @rename, [$oldiid, $newiid];
	}
	if (@rename) {
		unless($wrkdir) {
			$wrkdir = tempdir(CLEANUP => 1);
		}
		else {
			unless(-d $wrkdir) {
				make_path $wrkdir;
			}
		}
		open my $fout, ">$wrkdir/rename.txt" or die "Cannot write to rename";
		my $fbase = basename($vcf);
		system(qq|bcftools reheader -s $wrkdir/rename.txt -o $wrkdir/$fbase $vcf|);
		swap_files("$wrkdir/$fbase", $vcf);
		return 1;
	}
	else {
		return 0;
	}
}


=head2 plink_run CMD

Run plink command line, return undef if errors found in log file.

=cut

sub plink_run {
	my ($cmd) = @_;
	system($cmd);
	my $outfile;
	if ($cmd =~ /\-\-out\s+(\S+)/) {
		$outfile = $1;
	}
	else {
		$outfile = 'plink';
	}
	croak "Cannot find log file" unless -f "$outfile.log";
	my @errs;
	open my $fin, "$outfile.log" or croak "Cannot open $outfile.log";
	while(<$fin>) {
		if (/Error/i) {
			push @errs, $_;
		}
	}
	if (@errs) {
		return 0;
	}
	else {
		return 1;
	}
}

=head2 file_check

Test if file exist on a remote server.

=cut

sub file_exists {
	my ($server, $path) = @_;
	if ($server =~ /localhost/) {
		if (-f $path) {
			return 1;
		}
		else {
			return 0;
		}
	}
	else {
		system qq|ssh $server 'test -e $path'|;
		my $rc = $? >> 8;
		if ($rc) {
			return 0;
		}
		else {
			return 1;
		}
	}
}

=head2 extract_bam

Extract given regions from BAM/CRAM files.

Files can be stored in a remote server (accessible via SSH). 

Users need to specify the absolute path to samtools and reference sequence (for cram).

Region specification can be either "chr:st-ed" string, or a list of such string,
or a Genome::Ranges object.

=cut 

sub extract_bam {
	my ($file, $output, $argref) = @_;
	my $arg = merge_opts($argref, server => "localhost", outfmt => "bam",
		samtools => undef, fasta => undef, region => undef);
	croak "Cannot find data file $file" unless file_exists($arg->{server}, $file);
	croak "Cannot find samtools on $arg->{server}" 
		unless defined $arg->{samtools} && file_exists($arg->{server}, $arg->{samtools});
	croak "Must provide region specification or a (local) bed file" unless defined $arg->{region};

	my ($outsuffix, $outopt);
	if ($arg->{outfmt} eq 'bam') {
		$outsuffix = 'bam';
		$outopt = '-b'
	}
	elsif ($arg->{outfmt} eq 'cram') {
		$outsuffix = 'cram';
		$outopt = '-C';
	}
	else {
		croak "Cannot recognize output format: $arg->{outfmt}";
	}
	if ($file =~ /\.cram$/) {
		unless (defined $arg->{fasta} && file_exists($arg->{server}, $arg->{fasta})) {
			croak "Cannot find genome reference file on $arg->{server}";
		}
		(my $idxfile = $file) =~ s/cram$/crai/;
		unless (file_exists($arg->{server}, "$file.crai") || 
				file_exists($arg->{server}, $idxfile)) {
			croak "Cannot find data file index $file.crai or $idxfile on $arg->{server}";
		}
		if (-f $arg->{region}) {
			# If regions are provided in file, it must be readable to the server before this operation
			# here we use a pipeline to redirect 
			system(qq{cat $arg->{region} | }.
				qq|ssh $arg->{server} '$arg->{samtools} view |.
				qq|-T $arg->{fasta} $outopt -h -L - -M -o - $file' |.
				qq|>$output.$outsuffix|);
		}
		elsif ($arg->{region} =~ /^\w:\d+-\d+$/) {
			system(qq|ssh $arg->{server} '$arg->{samtools} view |.
				qq|-T $arg->{fasta} $outopt -h -o - $file $arg->{region}' |.
				qq|>$output.$outsuffix|);
		}
		else {
			croak "Cannot recognize region: $arg->{region}";
		}
	}
	elsif ($file =~ /\.bam$/) {
		(my $idxfile = $file) =~ s/bam$/bai/;
		unless (file_exists($arg->{server}, "$file.bai") || 
				file_exists($arg->{server}, $idxfile)) {
			croak "Cannot find data file index $file.bai or $idxfile on $arg->{server}";
		}
		if (-f $arg->{region}) {
			system(qq{cat $arg->{region} | }.
				qq|ssh $arg->{server} '$arg->{samtools} view |.
				qq|-T $arg->{fasta} $outopt -h -L $arg->{region} -M -o - $file'|.
				qq|>$output.$outsuffix|);
		}
		elsif ($arg->{region} =~ /^\w:\d+-\d+$/) {
			system(qq|ssh $arg->{server} '$arg->{samtools} view |.
				qq|-T $arg->{fasta} $outopt -h -o - $file $arg->{region}' |.
				qq|>$output.$outsuffix|);
		}
		else {
			croak "Cannot recognize region: $arg->{region}";
		}
	}
	else {
		croak "Data file $file has unsupported suffix";
	}
	if (-f "$output.$outsuffix") {
		system("samtools index $output.$outsuffix");
		return 1;
	}
	else {
		croak "Cannot extract reads from $file";
	}
}

=head2 iter_sigdat

Iterate over sample array signals for given chromosome or region

=cut

sub iter_sigdat {
	my ($sigdir, $iid, $chrom, $start, $end) = @_;
	my $spec;
	if (defined $start && defined $end) {
		$spec = "$chrom:$start-$end";
	}
	else {
		$spec = $chrom;
	}
	open my $fin, "tabix -h -s 2 -b 3 -e 3 $sigdir/$iid.txt.gz $spec |"
		or die "Cannot open tabix pipe for $iid $spec";
	my @fields;
	{
		my $fin = open_file("$sigdir/$iid.txt.gz");
		my $header = <$fin>; chomp($header);
		@fields = split(/\t/, $header);
	}
	my $f_lrr = first { /\.Log R Ratio$/ } @fields;
	my $f_baf = first { /\.B Allele Freq$/ } @fields;
	unless(defined $f_lrr && defined $f_baf) {
		die "Cannot find LRR and/or BAF field in $sigdir/$iid.txt.gz!";
	}
	unless ($f_lrr eq "$iid.Log R Ratio" && $f_baf eq "$iid.B Allele Freq") {
		my ($sampid) = ($f_lrr =~ /^(.+)\.Log R Ratio$/);
		warn "Sample in the signal data has been renamed $iid => $sampid";
	}
	my $it = iter_file($fin, { fsep => qr/\t/, alias => { $f_lrr => "LRR", $f_baf => "BAF" } });
	return $it;
}

=head2 iter_normdp

Iterate over CLAMMS-style normalized depth of coverage data.

=cut

sub iter_normdp {
	my ($sigdir, $iid, $chrom, $start, $end) = @_;
	my $spec;
	if (defined $start && defined $end) {
		$spec = "$chrom:$start-$end";
	}
	else {
		$spec = $chrom;
	}
	open my $fin, "tabix -p bed $sigdir/$iid.norm.cov.bed.gz $spec |"
		or die "Cannot open tabix pipe for $iid $spec";
	my $it = iter_file($fin, { fsep => qr/\t/, header => 0, 
							   fields => [qw|Chrom Start End NormDP DipMean Sigma|] });
	return $it;
}

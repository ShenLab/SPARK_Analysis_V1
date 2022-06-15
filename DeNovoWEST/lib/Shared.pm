package Shared;

use strict;
use warnings;
use Carp;
use IO::Dir;
use Tie::IxHash;
use Data::Dumper;
use Storable qw|dclone|;
use Cwd qw|abs_path|;
use FindBin qw|$Bin $Script|;
use File::Basename;
use File::Which qw|which|;
use File::Temp qw|tempdir|;
use File::Path qw|make_path|;
use Utils::Stat qw|mean|;
use List::Util qw|sum min max first|;
use List::MoreUtils qw|uniq all any|;
use Perl6::Slurp;
use Config::Std;
use List::BinarySearch qw|:all|;
use Storable qw|dclone|;
use Genet::Ped;
use String::Random;
use Graph::Directed;
use Graph::Undirected;
use Lingua::EN::Numbers::Ordinate;
use Utils::Parser qw|sql_query|;
use Utils::Hash qw|merge_opts peek_hash chk_default|;
use Utils::File::Iter qw|iter_file|;
use Utils::List qw|all_combs all_pairs|;
use Genome::UCSC qw|hg_chr|;
use Genome::Ranges qw|validate_elem range_overlap|;
use Genome::UCSC::BinKeeper;

use Wrapper qw|file_exists|;

use base qw|Exporter|;
our @EXPORT_OK = qw|check_conf read_rename read_remove read_list read_dir 
read_sex read_geneset read_twinsibs read_twindups read_sibgraph 
parse_fstr parse_tabfile merge_cols parse_filters parse_ranges parse_bins parse_calibr 
chain_renames merge_tabfiles slurp_xref expand_dat struct_dat fam_rels reltype fam_kins
read_chrlen split_chrs comb_intervals read_probes find_proberng region_bk|;


=head1 SYNPOSIS

Utility functions related to read/write data files.

=head2

Check the data read from config file with config file template 

The default template is presumed to have the same file name as the script, and have conf suffix.

=cut

sub check_conf {
	my ($conf, $template) = @_;
	unless(defined $template) {
		$template = "$Bin/$Script";
		$template =~ s/\.pl$/.conf/;
	}
	unless(-f $template) {
		croak "Cannot find config file template!";
	}
	read_config $template => my %ref;
	foreach my $section (sort keys %$conf) {
		next if $section eq 'SGE' || $section eq 'BASH';
		foreach my $key (sort keys %{$conf->{$section}}) {
			unless(defined $ref{$section} && defined $ref{$section}{$key}) {
				croak "$section.$key is not a valid config option!";
			}
		}
	}
	return 1;
}


=head2 read_list LIST, OPTIONS

The list file should contain two columns: path to the data file and sample ID
Sample ID can be different from IDs encoded in the data file (e.g., SM field in BAM file).
If second column is not provided, the basename of the file (after stripping off suffix, if provided) 
will be used as sample ID. Alternatively, provide a callback to extract sample ID from file.
Priority to assign sample ID will be given to the second column, followed by callback, then basename.

Options:
  * suffix: A list of or a single allowable sufficies. If a list is provided, file types will be 
            determined from the suffix, but ONLY ONE file type should be allowed.
  * rename: A sample rename list, to substitute the old ID with the new one.
  * multi:  To allow multiple data files for one sample. By default, one sample is required to have 
            only one associated data file.
  * remove: A sample removal list. Should use sample ID after rename.
  * ignore: regex, ignore the lines matching this pattern in the list.
  			By default, ignores /^#/ and empty lines.
  * callback: Callback function to return sample IDs in the data file.
  * remote: Hostname of the remote server that files are hosted (SSH)
  * notest: Skip testing the existance of files, faster but risky

=cut 

sub read_list {
	my ($list, $argref) = @_;
	my $arg = merge_opts($argref, suffix => undef, rename => undef, multi => undef,
		ignore => qr/^#/, remove => undef, remote => undef, notest => undef, callback => undef);

	if (-f $list) {
		print STDERR "Reading list $list\n";
	}
	elsif (-d $list) {
		carp "The provided list file is a directory";
		return read_dir($list, $argref);
	}
	else {
		croak "Cannot find file $list";
	}

	my %rename = read_rename($arg->{rename});
	my %remove = read_remove($arg->{remove});

	#my %files;
	tie my %files, 'Tie::IxHash';
	my $count = 0;
	my $filetype;
	{
		open my $fin, $list or die "Cannot open data file list";
		while(<$fin>) {
			next if /^\s*$/;
			if (defined $arg->{ignore}) {
				next if /$arg->{ignore}/;
			}
			my ($datafile, $iid) = split(/\s+/, $_);

			unless($arg->{notest}) {
				if ($arg->{remote}) {
					croak "Cannot find data file $datafile on server $arg->{remote}" 
						unless file_exists($arg->{remote}, $datafile);
				}
				else {
					croak "Cannot find data file $datafile" unless -f $datafile;
				}
			}
			
			my $basename = basename($datafile);
			my ($fbase, $ftype) = _parse_filename($basename, $arg->{suffix});

			unless (defined $fbase && defined $ftype) {
				croak "Cannot determine file type for $basename";
			}

			if (defined $filetype) {
				unless ($filetype eq $ftype) {
					croak "Only one file type is allowed in a list: $fbase --> $filetype <> $ftype";
				}
			}
			else {
				$filetype = $ftype;
			}

			unless ($iid) {
				if ($arg->{callback}) {
					my @iids = $arg->{callback}->($datafile);
					unless (@iids == 1) {
						croak "Data file $datafile contains multiple sample IDs";
					}
					else {
						$iid = $iids[0];
					}
				}
				else {
					$iid = $fbase;
				}
			}

			$iid = $rename{$iid} if defined $rename{$iid};
			next if defined $remove{$iid};

			if ($arg->{multi}) {
				push @{$files{$iid}} => $datafile;
			}
			else {
				if (defined $files{$iid}) {
					croak "Data file for $iid was already found!";
				}
				else {
					$files{$iid} = $datafile;
				}
			}
			$count ++;
			if ($count % 1000 == 0) {
				print STDERR "$count $ftype files read\n"
			}
		}
	}
	unless (keys %files) {
		croak "No files can be found in list $list";
	}
	if (wantarray) {
		return (\%files, $filetype);
	}
	else {
		return \%files;
	}	
}



=head2 read_dir DIR, OPTIONS

Read sample data from a directory. 

The function provide the same interface as read_list.

=cut

sub read_dir {
	my ($indir, $argref) = @_;
	my $arg = merge_opts($argref, suffix => undef, rename => undef, multi => undef,
		ignore => undef, remove => undef, callback => undef);

	if (-d $indir) {
		print STDERR "Reading directory $indir\n";
	}
	elsif (-f $indir) {
		croak "The provided directory is a file";
	}
	else {
		croak "Cannot find directory $indir";
	}

	my %rename = read_rename($arg->{rename});
	my %remove = read_remove($arg->{remove});

	# Read files
	#my %files;
	tie my %files, 'Tie::IxHash';
	my $filetype;
	{
		my @files = grep { !/^\./ && -f "$indir/$_" } IO::Dir->new($indir)->read();
		my $count = 0;
		foreach my $basename (@files) {
			my $datafile = "$indir/$basename";
			my ($fbase, $ftype) = _parse_filename($basename, $arg->{suffix});

			next unless defined $ftype && defined $fbase;

			if (defined $filetype) {
				unless ($filetype eq $ftype) {
					croak "Only one file type is allowed in a list: $fbase --> $filetype <> $ftype";
				}
			}
			else {
				$filetype = $ftype;
			}
			my $iid;
			if ($arg->{callback}) {
				my @iids = $arg->{callback}->($datafile);
				unless (@iids == 1) {
					croak "Data file $datafile contains multiple sample IDs";
				}
				else {
					$iid = $iids[0];
				}
			}
			else {
				$iid = $fbase;
			}

			$iid = $rename{$iid} if defined $rename{$iid};
			next if defined $remove{$iid};

			if ($arg->{multi}) {
				push @{$files{$iid}} => $datafile;
			}
			else {
				if (defined $files{$iid}) {
					croak "Data file for $iid was already found!";
				}
				else {
					$files{$iid} = $datafile;
				}
			}
			$count ++;
			if ($count % 1000 == 0) {
				print STDERR "$count $ftype files read\n"
			}
		}
	}
	unless (keys %files) {
		croak "No files can be found in directory $indir";
	}
	if (wantarray) {
		return (\%files, $filetype);
	}
	else {
		return \%files;
	}
}

sub _parse_filename {
	my ($basename, $suffix) = @_;
	my ($ftype, $fbase);
	if (defined $suffix) {
		if (ref $suffix) {
			if (ref $suffix eq 'HASH') {
				my %suffix = %$suffix;
				while(my ($type, $regex) = each %suffix) {
					if ($basename =~ /^(.+)$regex/) {
						$fbase = $1;
						$ftype = $type;
						last;
					}
				}
			}
			elsif (ref $suffix eq 'ARRAY') {
				foreach my $type (@{$suffix}) {
					if ($basename =~ /^(.+)\.($type)$/) {
						$fbase = $1;
						$ftype = $2;
						last;
					}
				}
			}
			elsif (ref $suffix eq 'Regexp') {
				if ($basename =~ /^(.+)$suffix$/) {
					$fbase = $1;
					$ftype = 'Unknown';
				}
			}
			else {
				croak "Incorrect object type of suffix: should be HASH/ARRAY/Regexp";
			}	
		}
		else {
			my $type = $suffix;
			if ($basename =~ /^(.+)\.($type)$/) {
				$fbase = $1;
				$ftype = $2;
			}
		}
	}
	else {
		$fbase = $basename;
		$ftype = "NA";
	}
	#unless (defined $fbase && defined $ftype) {
	#	croak "Cannot determine the file type for $basename";
	#}
	# print $fbase, "\t", $ftype, "\n";
	return ($fbase, $ftype);
}

# Helper functions
sub read_rename {
	my ($file) = @_;
	my %rename;
	return unless defined $file;
	if (ref $file eq 'HASH') {
		%rename = %$file;
	}
	else {
		open my $fin, $file or die "Cannot open rename list: $file";
		while(<$fin>) {
			next if /^\s*$/;
			my ($oldid, $newid) = (split)[0,1];
			if (defined $rename{$oldid}) {
				croak "Alias for $oldid has already been defined";
			}
			$rename{$oldid} = $newid;
		}
	}
	if (wantarray) {
		return %rename;
	}
	else {
		return \%rename;
	}
}

sub read_remove {
	my ($file) = @_;
	return unless defined $file;
	my %remove;
	if (ref $file eq 'HASH') {
		%remove = %$file;
	}
	elsif (ref $file eq 'ARRAY') {
		@remove{@$file} = @$file;
	}
	else {
		open my $fin, $file or die "Cannot open removal list: $file";
		while(<$fin>) {
			next if /^\s*$/;
			my $iid = (split)[0];
			$remove{$iid} = 1;
		}
	}
	if (wantarray) {
		return %remove;
	}
	else {
		return \%remove;
	}
}

sub read_sex {
	my ($file) = @_;
	my %sex;
	if (ref $file eq 'HASH') {
		%sex = %$file;
		unless(all { $_ eq '1' || $_ eq '2' || $_ =~ /^m/i || $_ =~ /^f/i } values %sex) {
			croak "Incorrect code for sample sex";
		}
		foreach my $iid (keys %sex) {
			my $sexcode = $sex{$iid};
			if ($sexcode =~ /^m/i) {
				$sex{$iid} = 1;
			}
			elsif ($sexcode =~ /^f/i) {
				$sex{$iid} = 2;
			}
			else {
				croak "Cannot determine sex $sexcode for $iid";
			}
		}
	}
	else {
		open my $fin, $file or die "Cannot open rename list: $file";
		while(<$fin>) {
			next if /^\s*$/;
			my ($iid, $sex) = (split)[0,1];
			if (defined $sex{$iid}) {
				croak "Sex $iid has already been defined";
			}
			if ($sex eq '1' || $sex =~ /^m/i) {
				$sex{$iid} = 1;
			}
			elsif ($sex eq '2' || $sex =~ /^f/i) {
				$sex{$iid} = 2;
			}
			else {
				croak "Cannot determine sex $sex for $iid";
			}
			#$sex{$iid} = $sex;
		}
	}
	if (wantarray) {
		return %sex;
	}
	else {
		return \%sex;
	}
}


=head2 read_geneset File, Fields

Gene sets can be read from files in the following five ways:

1. A list of gene in a given column of a table file
	File=FilePath(:SetName)    <-- FilePath cannot contain [:]; SetName is optional, default to file basename without suffix
	Fields=ColName or ColNum   <-- ColName cannot contain any of [,"']
2. Genes categorized in a table file
	File=FilePath
	Fields=GeneID,GeneClass    <-- GeneClass must be a categorical variable
3. Subset of genes defined by SQL like filtering expression in a gene table
This method is more flexible than 2 to define and select gene sets
	File=FilePath
	Fields=GeneID,"Filter Expression 1":SetName1,'Filter Expression 2':SetName2
    <-- Filter expression are quoted and followed by an optional gene set name.
4. Gene sets defined by GMT file (must have .gmt extension to the file name)
	File=filepath   
	Fields=SetName1:Alias1,SetName2:Alias2   <--Can select one or more sets here, may be left empty
5. Gene sets defined by GMX file (must be an excel file)
	File=FilePath(:1)    <--Should be sheet number after file name
	Fields=SetName1,SetName2:Alias2 <- Should set names defined in the first row

In all except Method5, input file should be tab-delmited text file. A gene set can only be 
defined in one file, not distributed across many files.

Under scalar context it will return a hash of gene sets found $geneset{Set} = [G1, G2, ...].
Under array context it will also return a hash of hash storing gene set membership for each
gene $gsmember{GeneID}{Set} = 1. %gsmember can store the number of times each gene appear in 
the set and can be used as weights.

=cut

sub read_geneset {
	my ($file, $fields) = @_;
	my (@gsfiles, @gsfields);
	if (ref $file eq 'ARRAY') {
		@gsfiles = @$file;
		@gsfields = @$fields;
	}
	else {
		@gsfiles = ($file);
		@gsfields = ($fields);
	}

	tie my %geneset, "Tie::IxHash";
	for(my $ii = 0; $ii < @gsfiles; $ii ++) {
		my ($filepath, $filesupp) = split(':', $gsfiles[$ii]);

		my (%gs, %gskeep);
		if ($filepath =~ /\.gmt$/) {
			print STDERR "Reading gene sets from GMT file: $filepath\n";
			if ($gsfields[$ii]) {
				%gskeep = parse_fstr($gsfields[$ii]);
			}
			open my $fin, $filepath or die "Cannot open $filepath";
			while(<$fin>) {
				my @dat = split;
				unless(@dat >= 3) {
					die "Gene set $dat[0] contains no gene?"
				}
				if (defined $gs{$dat[0]}) {
					die "Gene set $dat[0] has already been defined in GMT file";
				}
				$gs{$dat[0]} = [ @dat[2..$#dat] ];
			}
		}
		elsif ($filepath =~ /\.(xlsx|xls|ods)$/) {
			print STDERR "Reading gene sets from GMX spreadsheet: $filepath\n";
			unless(defined $filesupp) {
				die "Must provide sheet num or name for spreadsheet";
			}
			if ($gsfields[$ii]) {
				%gskeep = parse_fstr($gsfields[$ii]);
			}
			my $book = Spreadsheet::Read->new($filepath, strip => 3);
			my $sheet = $book->sheet($filesupp);
			for my $col (1..$sheet->maxcol) {
    			my $setname = $sheet->cell($col, 1);
    			if (defined $gs{$setname}) {
    				die "Gene set $setname has already been defined in GMX spreadsheet";
    			}
			    for(my $ii = 2; $ii <= $sheet->maxrow; $ii ++) {
					my $gid = $sheet->cell($col, $ii);
					last unless $gid;
					push @{$gs{$setname}} => $gid;
    			}
			}
		}
		else {
			my @filefields = split(',', $gsfields[$ii]);
			if (@filefields == 1) {
				print STDERR "Reading gene set from gene list $filepath\n";
				my $set;
				if ($filesupp) {
					$set = $filesupp;
				}
				else {
					$set = (split(q|\.|, basename($filepath)))[0];
				}
				my ($it, $fnames, $keyfields) = parse_tabfile($filepath, $gsfields[$ii], 1, 1);
				while(my $dat = $it->()) {
					my ($gid) = @{$dat}{@$keyfields};
					push @{$gs{$set}}, $gid;
				}
			}
			else {
				if (@filefields == 2 && $filefields[1] !~ /^['"]/) {
					print STDERR "Reading gene set from two columns of gene table $filepath\n";
					if ($filesupp) {
						warn "Filesupp $filesupp will be ignored";
					}
					my ($it, $fnames, $keyfields) = parse_tabfile($filepath, $gsfields[$ii], 2, 2);
					while(my $dat = $it->()) {
						my ($gid, $set) = @{$dat}{@$keyfields};
						push @{$gs{$set}}, $gid;
					}
				}
				elsif (@filefields >= 2 &&
					   (all { $filefields[$_] =~ /^['"]/ } 1..$#filefields)) { 
					print STDERR "Reading gene set from table file with filter expression\n";
					if ($filesupp) {
						warn "Filesupp $filesupp will be ignored";
					}
					my ($it, $fnames) = iter_file($filepath, { fsep => qr/\t/ });
					# First field must be GeneID
					my $f_gid = shift @filefields;
					my (@filters, @labels);
					foreach my $filefield (@filefields) {
						my ($filter, $label) = split(':', $filefield);
						print $filter, "\t", $label, "\n";
						if (defined $label) {
							push @labels, $label;
						}
						else {
							push @labels, $filter;
						}
						$filter =~ s/^["']//; $filter =~ s/["']$//;
						my ($cb, $tokens) = sql_query($filter, 1);
						foreach my $tok (@$tokens) {
							if ($tok->[0] eq 'FIELD') {
								unless(grep { $_ eq  $tok->[1] } @$fnames) {
									die "Cannot find field $tok->[1] in input file: $filepath";
								}
							}
						}
						unless(grep { $_ eq $f_gid } @$fnames) {
							die "Cannot find gene ID field $f_gid from input file: $filepath";
						}
						push @filters, $cb;
					}
					while (my $dat = $it->()) {
						for (my $jj = 0; $jj < @filters; $jj ++) {
							if ($filters[$jj]->($dat)) {
								push @{$gs{$labels[$jj]}}, $dat->{$f_gid}; 
							}
						}
					}
				}
				else {
					die "Cannot parse gene set file: $gsfiles[$ii]\nFields: $gsfields[$ii]";
				}
			}
		}
		
		if (%gskeep) {
			foreach my $set (keys %gskeep) {
				unless(defined $gs{$set}) {
					die "Cannot find gene set $set from $filepath";
				}
				if (defined $geneset{$gskeep{$set}}) {
					die "Gene set $gskeep{$set} has already been defined in previous files";
				}
				$geneset{$gskeep{$set}} = $gs{$set};
			}
		}
		else {
			foreach my $set (keys %gs) {
				if (defined $geneset{$set}) {
					die "Gene set $set has already been defined before";
				}
				$geneset{$set} = $gs{$set};
			}
		}
	}
	# Then store gene set membership for genes if needed, it stores the occurence of genes
	if (wantarray) {
		my %gsmember;
		while(my ($set, $genes) = each %geneset) {
			foreach my $gid (@$genes) {
				$gsmember{$gid}{$set} ++;
			}
		}
		return (\%geneset, \%gsmember);
	}
	else {
		return \%geneset;
	}
}


=head2 read_twinsibs

Read twin and sibs.

Twin or sib list is a two column file, each row is a pair of samples.
We will store those pairwise relations in a undirected graph and find out
connected components. So for each sample, we can get all his twin or sibs
from the connected graph component.

=cut

sub read_sibgraph {
	my $g = Graph::Undirected->new;
	foreach my $list (grep { defined $_ } @_) {
		next unless defined $list;
		if (-f $list) {
			print STDERR "Reading twin/sib pairs from file: $list\n";
			open my $fin, $list or die "Cannot open $list";
			while(<$fin>) {
				my @iids = split;
				unless(@iids == 2) {
					die "Incorrect number of columns of twin or sib pair list!";
				}
				$g->add_edge($iids[0], $iids[1]);
			}	
		}
		elsif (ref $list eq 'ARRAY') {
			print STDERR "Reading twin/sib pairs from array\n";
			foreach my $pair (@$list) {
				unless(@$pair == 2) {
					die "Incorrect number of elements for the pair!";
				}
				$g->add_edge($pair->[0], $pair->[1]);
			}
		}
		elsif (ref $list eq 'HASH') {
			print STDERR "Reading twin/sib pair from 2D-Hash\n";
			foreach my $pair (sort keys %$list) {
				my @iids = split($;, $pair);
				unless(@iids == 2) {
					die "Incorrect number of elements for the pair: $pair!";
				}
				$g->add_edge($iids[0], $iids[1]);
			}
		}
		else {
			die "Cannot recognize the format of twin/sib pair list";
		}
	}
	return $g;
}

# read_twinsibs store all other twin/sibs for each sample
sub read_twinsibs {
	my $g = read_sibgraph(@_);
	my %twinsibs;
	my @ccs = $g->connected_components();
	print STDERR "Number of connnected components=", scalar(@ccs), "\n";
	foreach my $cc (@ccs) {
		unless(@$cc >= 2) {
			die "Incorrect number of samples in a connected twin-sib cluster!";
		}
		foreach my $iid (sort @$cc) {
			$twinsibs{$iid} = [grep { $_ ne $iid } sort @$cc];
		}
	}
	if (wantarray) {
		return %twinsibs;
	}
	else {
		return \%twinsibs;
	}
}


# read_twindups will store the primary ID for all other members in the twin/dup cluster
sub read_twindups {
	my $g = read_sibgraph(@_);
	my %twindups;
	my @ccs = $g->connected_components();
	print STDERR "Number of connnected components=", scalar(@ccs), "\n";
	foreach my $cc (@ccs) {
		unless(@$cc >= 2) {
		#	print Dumper $cc; 
			warn "More than two samples in a connected twin-dups cluster!";
		}
		#foreach my $iid (sort @$cc) {
		#	$twindups{$iid} = [grep { $_ ne $iid } sort @$cc];
		#}
		my @samps = sort @$cc;
		my $mainid = shift @samps; 
		foreach my $secid (@samps) {
			if (defined $twindups{$secid}) {
				die "Secondary ID $secid has been renamed before: $twindups{$secid}!"
			}
			$twindups{$secid} = $mainid;
		}
	}
	if (wantarray) {
		return %twindups;
	}
	else {
		return \%twindups;
	}
}


=head2 chain_renames FILES...

Chain multiple rename files to create a final mapping list from the original
ID to the final ID.

This function will use read_rename function to read the mapping list from
the individual mapping file. 

For chaining purpose, We also need to check that individual rename list
does not have 1-to-many mapping.

=cut

sub chain_renames {
	my @files = @_;
	my @maps;
	foreach my $file (@files) {
		my $idmap = read_rename($file);
		my %revidmap;
		while(my ($oldid, $newid) = each %$idmap) {
			if (defined $revidmap{$newid}) {
				die "Multiple samples mapped to $newid in $file";
			}
			$revidmap{$newid} = $oldid;
		}
		unshift @maps, \%revidmap;
	}
	my @idpairs;
	# now chaining, from the final new ID to original ID
	for(my $jj = 0; $jj < @maps; $jj ++) {
		foreach my $finalid (sort keys %{$maps[$jj]}) {
			my $previd = $finalid;
			for(my $ii = $jj; $ii < @maps; $ii ++) {
				if (defined $maps[$ii]{$previd}) {
					my $tmpid = $maps[$ii]{$previd};
					delete $maps[$ii]{$previd};
					$previd = $tmpid;
				}
			}
			push @idpairs, [$previd, $finalid];
		}
	}
	return @idpairs;
}


=head2 parse_fstr STRING [, ORDER]

Parse field name specification string.

An example: 

Fields=IID,COHORT:Cohort,CHDCLASS:CHDClass,EM,NDD,SEX:Gender

will be parsed to hash ref: { IID => IID, COHORT => Cohort, CHDCLASS => CHDClass, EM => EM, NDD => NDD, SEX => Gender }

Can be used as alias for iter_file.

=cut

sub parse_fstr {
	my ($fstr, $order) = @_;
	my %fnames;
	if ($order) {
		tie %fnames, 'Tie::IxHash';
	}
	
	if (!defined $fstr) {
		return \%fnames;
	}

  	# Parse fields
  	my @fields;
  	if (ref $fstr eq 'ARRAY') {
  		#@fields = @$fstr;
  		foreach my $string (@$fstr) {
  			push @fields, split(',', $string);
  		}
  	}
  	else {
  		@fields = split(',', $fstr);
  	}
  	
  	foreach my $field (@fields) {
  		my ($fnm, $lab) = split(':', $field);
  		$lab = $fnm unless defined $lab;
		#print $fnm, "\t", $lab // $fnm, "\n";
		chk_default(\%fnames, $fnm, $lab);
		# $fnames{$fnm} = $lab // $fnm;
	}
	my @labs = values %fnames;
	my @luq = uniq sort @labs;
	if (@labs != @luq) {
		croak "Non-unique label names!";
	}
	if (wantarray) {
		if ($order) {
			warn "Order will be ignored if return hash";
		}
		return %fnames;
	}
	else {
		return \%fnames;
	}
}

=head2 parse_tabfile INFILE, FSTR, N_MIN[, N_MAX]

Commonly used subroutine for parsing tab separated data file.

Should provide a comma separated list of required fields FSTR.

Field names after N_MIN'th are optional.

=cut

sub parse_tabfile {
	my ($infile, $fieldstr, $n_min, $n_max) = @_;
	$n_max = $n_min unless defined $n_max;
	my $input_header = 1;
	my @input_fields = split(',', $fieldstr);
	croak "Must provide $n_min~$n_max fields" 
		unless @input_fields >= $n_min && @input_fields <= $n_max;
	if (all { /^\d+$/ } @input_fields) {
		$input_header = 0;
	}
	my ($it, $fnames) = iter_file($infile eq '-' ? \*STDIN : $infile , { fsep => qr/\t/, header => $input_header });
	for(my $ii = 0; $ii < @input_fields; $ii ++) {
		my $field = $input_fields[$ii];
		if ($ii < $n_min) {
			unless(grep { $field eq $_ } @$fnames) {
				croak "Cannot find field $field in the input file $infile";
			}
		}
		elsif ($ii >= $n_max) {
			croak "More than $n_max fields are given: $fieldstr";
		}
		else {
			unless(grep { $field eq $_ } @$fnames) {
				carp "Cannot find field $field in the input file";
				pop @input_fields;
			}
		}
	}

	return ($it, $fnames, \@input_fields);
}


=head2 merge_cols

When merging contents from multipe table files, the input tables may contain different columns.
We will take the union of columns  while keep the order of original columns as much as possible 
(using the algorithms shown below).

Ordering algorithm

1 2 3 4 5 6 
A B E F G H
1 2 3 4 5 6
A B D E G H
1.5 3 4.5 6
A   B C   D

1.2 2.3 4.5  3.5 4  5  6  6
A   B   C    E   F  G  H  D

=cut

sub merge_cols {
	my ($cols, $priority) = @_;

	# Union of all column names
	my @allcols = uniq sort map { @$_ } values %$cols;

	my (%posord, %posordmean, %grprnk);
	foreach my $label (sort { $priority->{$a} <=> $priority->{$b} } keys %$cols) {
		my @fields = @{$cols->{$label}};
		for (my $ii = 0; $ii < @fields; $ii ++) {
			push @{$posord{$fields[$ii]}}, $ii * scalar(@allcols)/scalar(@fields);
			unless (defined $grprnk{$fields[$ii]}) {
				$grprnk{$fields[$ii]} = $priority->{$label};
			}
		}
	}
	foreach my $field (keys %posord) {
		$posordmean{$field} = mean($posord{$field});
	}
	my @mergecols = sort { $posordmean{$a} <=> $posordmean{$b} || $grprnk{$a} <=> $grprnk{$b} } keys %posordmean;

	if(wantarray) {
		return @mergecols;
	}
	else {
		return \@mergecols;
	}
}


=head2 parse_filters

Parse a group of filters.

=cut

sub parse_filters {
	my ($labels, $filters) = @_;
	# Find out all filter names and expression for defining sample subsets
	my (@filternames, @filterexprs);
	unless (defined $labels && defined $filters) {
		die "Must provide both filter label and expression";
	}
	if (ref $labels eq 'ARRAY') {
		push @filternames => @{$labels};
		unless(ref $filters eq 'ARRAY' &&
			@{$filters} == @{$labels}) {
			print Dumper $filters, $labels;
			die "Numbers of filter labels and expressions do not match";
		}
		push @filterexprs => @{$filters};
	}
	else {
		push @filternames => $labels;
		if (ref $filters eq 'ARRAY') {
			die "Numbers of filter labels and expressions do not match";
		}
		push @filterexprs => $filters;
	}
	# Parse filter expression and store relevant fields
	my (@callbacks, %fields);
	for(my $ii = 0; $ii < @filterexprs; $ii ++) {
		$filterexprs[$ii] =~ s/^["']//; $filterexprs[$ii] =~ s/["']$//;
		print STDERR $filternames[$ii], ": ", $filterexprs[$ii], "\n";
		my ($cb, $tokens) = sql_query($filterexprs[$ii], 1);
		push @callbacks, $cb;
		foreach my $tok (@$tokens) {
			$fields{$tok->[1]} = 1 if $tok->[0] eq 'FIELD';
		}
	}
	return (\@filternames, \@callbacks, \%fields);
}

=head2 parse_ranges

Valid range format: 1~5, =5, >6...

=cut 

sub parse_ranges {
	my ($rngstr) = @_;
	tie my %counts, 'Tie::IxHash';
	foreach my $range (split(',', $rngstr)) {
		if (defined $counts{$range}) {
			die "Gene counts range for $range has been defined!";
		}
		if ($range =~ /^\d+$/) {
			$counts{"=".$range} = [$range];
		}
		elsif ($range =~ /^(\d+)[\-\~](\d+)$/) {
			my ($min, $max) = ($1, $2);
			unless($min < $max) {
				die "Invalid range: $range";
			}
			$counts{$range} = [$min, $max];
		}
		elsif ($range =~ /^([>=<]+)(\d+)$/) {
			my ($mod, $num) = ($1, $2);
			if ($mod eq '>') {
				$counts{$range} = [$num+1, 3e10];
			}
			elsif ($mod eq '>=') {
				$counts{$range} = [$num, 3e10];
			}
			elsif ($mod eq '<') {
				$counts{$range} = [0, $num-1];
			}
			elsif ($mod eq '<=') {
				$counts{$range} = [0, $num];
			}
			elsif ($mod eq '=' || $mod eq '==') {
				$counts{$range} = [$num];
			}
			else {
				die "Invalid range: $range";
			}
		}
		else {
			die "Cannot recognize range: $range";
		}
	}
	if (wantarray) {
		return %counts;
	}
	else {
		return \%counts;
	}
}

=head2 parse_bins

Valid bin format: 1:21:5,30 (Start:End:Step)

Missing data should be represented as "." in the bin, and can only appear at the beginning.

The score specified by the bin must be ordered!

=cut 

sub parse_bins {
	my ($binstr) = @_;
	my @sbins;
	foreach my $score (split(',', $binstr)) {
		if ($score =~ /:/) {
			my ($start, $end, $step) = split(':', $score);
			unless(defined $start && defined $end && defined $step) {
				die "Cannot parse Start:End:Step from $score";
			}
			unless($start =~ /^-?\d[\.0-9]*$/ && $end =~ /^-?\d[\.0-9]*$/ && $step =~ /^\d[\.0-9]*$/) {
				die "Incorrect Start:End:Step: $score";
			}
			my $nsteps = int(($end-$start)/$step);
			for(my $ii = 0; $ii <= $nsteps; $ii ++) {
				push @sbins, $start+$ii*$step;
			}
		}
		else {
			unless($score eq '.' || $score =~ /^-?\d[\.0-9]*$/) {
				die "Score is not a number: $score";
			}
			if ($score eq "." && @sbins > 0)  {
				die "Missing data can only appear at the beginning!";
			}
			push @sbins, $score;
		}
	}
	unless(scalar(@sbins) == scalar(uniq sort @sbins)) {
		die "Scores in the bin are not unique!";
	}
	my @nonmis = grep { $_ ne "." } @sbins;
	unless(@nonmis > 0) {
		die "No nonmissing score can be found in the bin!";
	}
	for(my $ii = 0; $ii < @sbins-1; $ii ++) {
		next if $sbins[$ii] eq ".";
		unless($sbins[$ii+1] > $sbins[$ii]) {
			die "Scores in the bin are not ordered!";
		}
	}
	if (wantarray) {
		return @sbins;
	}
	else {
		return \@sbins;
	}
}

=head2 parse_calibr

Parse calibration parameter grid.

=cut

sub parse_calibr {
	my ($config, $multi) = @_;
	#my $protofilt = $conf{Calibrate}{Filter};
	my $protofilt = $config->{Filter};
	unless(defined $protofilt) {
		die "Cannot find filter prototype!";
	}
	$protofilt =~ s/^['"]//; $protofilt =~ s/['"]$//;
	my %params;
	# Parse and check parameters
	foreach my $varname (grep { $_ ne 'Filter' } keys %$config) {
		# Check if the variable name can be found in the filter prototype
		# We require the parameters to be substituted appear only once in the proto-filter
		unless ($varname =~ /^[a-zA-Z]\w+$/) {
			die "Incorrect variable name $varname";
		}
		my @match = ($protofilt =~ /($varname)(?=\W|$)/g);
		unless(@match) {
			die "Cannot find variable $varname in filter: $protofilt";
		}
		else {
			if (@match > 1) {
				warn "Variable $varname appear multiple times in the filter: $protofilt";
			}
		}
		if ($config->{$varname} =~ /,/) {
			$params{$varname} = [split(',', $config->{$varname})];
		}
		elsif ($config->{$varname} =~ /:/) {
			my ($start, $end, $step) = split(':', $config->{$varname});
			unless(defined $start && defined $end && defined $step) {
				die "Must provide start:end:step for parameter $varname";
			}
			my $nsteps = int(($end-$start)/$step);
			for(my $ii = 0; $ii <= $nsteps; $ii ++) {
				push @{$params{$varname}}, $start+$ii*$step;
			}
		}
		else {
			$params{$varname} = [$config->{$varname}];
			#die "Cannot parse values from $varname: $config->{$varname}";
		}
		
	}
	# Enumerate all parameter combinations to form param grid
	my @paramgrid = all_combs(%params);
	if (wantarray) {
		return @paramgrid;
	}
	else {
		return \@paramgrid;
	}
}



=head2 slurp_xref

Slurp data from external reference
Multiple external references are allowed so is the field specification
Key will always be the first specified column in the fstr.
Return a data structure and field name list.

=cut

sub slurp_xref {
	my ($xref, $fstr) = @_;
	unless(defined $xref && defined $fstr) {
		die "Must provide external reference file and fields specification";
	}
	# Return results:
	# %xdat: will contain data, with key field alias $fkey (used to look for key field from input)
	# the data will be hashref with keys in @xfields
	my (%xdat, @xfields);

	my (@xrefs, @fstrs);
	if (ref $xref eq 'ARRAY') {
		unless(ref $fstr eq 'ARRAY' && @$xref == @$fstr) {
			print Dumper $xref, $fstr;
			die "Unequal length of xref files and field specs";
		}
		@xrefs = @$xref;
		@fstrs = @$fstr;
	}
	else {
		if (!defined $fstr) {
			die "Cannot find field spec for $xref";
		}
		elsif (ref $fstr eq 'ARRAY') {
			die "More than one field specs for one xref";
		}
		@xrefs = ($xref);
		@fstrs = ($fstr);
	}

	for(my $ii = 0; $ii < @xrefs; $ii ++) {
		$xref = $xrefs[$ii];
		$fstr = $fstrs[$ii];
	
		my $fields = parse_fstr($fstr, 1);

		my @fvals = values %$fields;
		my $fkey = shift @fvals;
		foreach my $fval (@fvals) {
			if (grep { $fval eq $_ } @xfields) {
				die "Field name $fval has been found previously!";
			}
		}
		push @xfields, @fvals;

		# Determine if we have header
		my $header = 1;
		if (all { /^\d+$/ } keys %$fields) {
			$header = 0;
		}
		# Iterate over the data file
		# TechNote: when both alias and suset are provided to iter_file, alias will be used to rename columns first before subseting
		# For files with many columns, it is likely to have name confliect.
		# So we will subset the data file using original ID, then rename them to alias
		print STDERR "Slurping $xref\n";
		my $iter;
		if ($xref =~ /\.csv$/) {
			$iter = iter_file($xref, { header => $header, select => $fields });
		}
		else {
			$iter = iter_file($xref, { fsep => qr/\t/, header => $header, select => $fields });
		}

		while(my $dat = $iter->()) {
			my $key =  $dat->{$fkey};
			foreach my $fval (@fvals) {
				# In case that multiple values found, only first one will be used.
				if (defined $xdat{$key}{$fval}) {
					print STDERR "slurp_xref $xref: Field $fval for $key has already been defined\n";
					next;
				}
				$xdat{$key}{$fval} = $dat->{$fval};
			}
		}
	}
	return (\%xdat, \@xfields);
}


=head2 expand_dat DATA, ARGS

Expand the data struct with multiple values packed in one cell.

Example: { Name => ID1, GQ => "99,50", GT => "Het,Het" }
=> { Name => ID1, GQ => 99, GT => "Het" }, { Name => ID1, GQ => 50, GT => "Het" }

Differece between fields and optional:
All columns in the "fields" should all expand or not expand, whereas fields specified
by "optional" may or may not expand. 

Under default mode, we require all expanded fields should have the same number of elements,
and croak err if field length differ. Under strict=>0 mode, we will fill in the last element
to the fields with smaller length (different from R's recycling rule).

=cut


sub expand_dat {
	my ($dat, $argref) = @_;
	my $arg = merge_opts($argref, fields => undef, optional => undef, sep => ',', strict => 1);

	unless ( (defined $arg->{fields} || defined $arg->{optional}) && defined $arg->{sep}) {
		die "Must provide fields and separator char";
	}

	my $fsep = qr/$arg->{sep}/;
	my @fields = defined $arg->{fields} ? @{$arg->{fields}} : ();
   	# Also support fields to expand optionally
   	# they can contain single value or multiple values
   	my @options = defined $arg->{optional} ? @{$arg->{optional}} : ();

    # Scan the provided fields to see if multiple values are packed
    my $dupflag;
    my @dupfields;
    if (any { $dat->{$_} =~ /$fsep/ } @fields) {
    	$dupflag = 1;
    	unless(all { $dat->{$_} =~ /$fsep/ } @fields) {
    		print Dumper $dat;
    		die "Not all fields that require expansion have multiple values";
    	}
    	push @dupfields, @fields;
    }
    if (any { $dat->{$_} =~ /$fsep/ } @options) {
    	$dupflag = 1;
    	push @dupfields, grep { $dat->{$_} =~ /$fsep/ } @options;
    }

    my @data;
    if ($dupflag) {
    	# Remove redundancies in dupfields
    	@dupfields = uniq sort @dupfields;

    	my $var = dclone($dat);
    	foreach my $fd (@dupfields) {
    		$var->{$fd} = [split($fsep, $var->{$fd}, -1)];
    	}
    	my $ndup = max (map { scalar(@{$var->{$_}}) } @dupfields); 
    	if ($arg->{strict}) {
    		unless(all { scalar(@{$var->{$_}}) == $ndup } @dupfields) {
    			my @errfields = grep { scalar(@{$var->{$_}}) != $ndup } @dupfields;
    			my %errdat;
    			@errdat{@errfields} = @{$var}{@errfields};
    			print Dumper \@dupfields, \@errfields, \%errdat;
    			print Dumper $var;
    			die "Not all expanded fields have the same length $ndup: ", join(",", @errfields);
    		}
    	}
    	else {
    		foreach my $fd (@dupfields) {
    			if (@{$var->{$fd}} != $ndup) {
    				print STDERR "Fill in field $fd\n";
    				for(my $ii = @{$var->{$fd}}; $ii < $ndup; $ii ++) {
    					push @{$var->{$fd}}, $var->{$fd}[-1];
    				}
    			}
    		}
    	}
    	for(my $ii = 0; $ii < $ndup; $ii ++) {
    		my %info;
    		foreach my $fd (sort keys %$var) {
    			if (grep { $fd eq $_ } @dupfields) {
    				$info{$fd} = $var->{$fd}[$ii];
    			}
    			else {
    				$info{$fd} = $var->{$fd};
    			}
    		}
    		push @data, \%info;
    	}
    }
    else {
    	push @data, $dat;
    }
    return @data;
}

=head2 struct_dat, ARGS

Restructure the data with provided key and values columns.

Example : { Chrom=>1, Pos=>2134, Ref=>'A', Alt => 'T', Gene=>'G1;G2',
            Trans => 'T1,T2;T3,T4', TransEff => 'Mis,Mis;Syn,Intron', 
            AAChg => 'p.G4T,pG4T;.,.' }

expand_dat => { Chrom=>1, Pos=>2134, Ref=>'A', Alt => 'T', Gene => 'G1'.
                Trans => "T1,T2", TransEff => 'Mis,Mis', AAChg => 'p.G4T,pG4T' },
              { Chrom=>1, Pos=>2134, Ref=>'A', Alt => 'T', Gene => 'G2'.
                Trans => "T3,T4", TransEff => 'Syn,Intron', AAChg => '.,.' }, 

struct_dat => { Chrom=>1, Pos=>2134, Ref=>'A', Alt => 'T', Gene => 'G1'.
                Trans => ["T1","T2"], TransEff => {T1 => 'Mis', T2 => 'Mis', 
                AAChg => {T1 => 'p.G4T', T2 => 'p.G4T' },
              { Chrom=>1, Pos=>2134, Ref=>'A', Alt => 'T', Gene => 'G2'.
                Trans => {["T3","T4"], TransEff => {T1 => 'Syn', T2 => 'Intron'}, 
                AAChg => {T1 => '.', T2 => '.'} }, 

There is an option for clone the entire data structure instead of overriding
the existing data.

Similar to expand_dat, we also provide a non-strict mode. Under this mode, when
the length of values does not equal to keys, the last element of value will be repeated

=cut


sub struct_dat {
	my ($dat, $argref) = @_;
	my $arg = merge_opts($argref, key => undef, value => undef, sep => ',', clone => undef, strict => 1);

	unless (defined $arg->{key} && defined $arg->{value}) {
		die "Must provide key and value fields";
	}
	unless (defined $arg->{sep}) {
		die "Must provide separator string";
	}
	my $fsep = qr/$arg->{sep}/;

	my $fkey = $arg->{key};
	my @fvals = ref $arg->{value} eq 'ARRAY' ? @{$arg->{value}} : ($arg->{value});

	unless(all { defined $dat->{$_} } ($fkey, @fvals)) {
		print Dumper $dat;
		die "Cannot find all key and value columns";
	}

	my $info;
	if ($arg->{clone}) {
		$info = dclone $dat;	
	}
	else {
		$info = $dat;
	}

	my @keys = split($arg->{sep}, $info->{$fkey});
	$info->{$fkey} = \@keys;
	foreach my $fval (@fvals) {
		my %dict;
		my @vals = split($arg->{sep}, $info->{$fval});
		unless(scalar(@keys) == scalar(@vals)) {
			if ($arg->{strict}) {
				print STDERR join($arg->{sep}, @keys), "\n";
				print STDERR join($arg->{sep}, @vals), "\n";
				die "Value field $fval does not contain the same number of elements as the key field";
			}
			else {
				if (scalar(@keys) > scalar(@vals)) {
					my $lastval = $vals[-1];
					push @vals => ($lastval) x (scalar(@keys)-scalar(@vals));
				}
				else {
					@vals = @vals[0..$#keys];
				}
			}
		}	
		@dict{@keys} = @vals;
		$info->{$fval} = \%dict;
	}
	return $info;
}



=head2 merge_tabfiles MAIN,F_MAIN, XREF,F_XREF, OUTFILE

Merge two table fileds, both files shoud be tab separated or csv files

File will be joined based on the values in shared fields (keys). Keys should be unique in both files.

Note: this will be "left outer join", missing values will be '.'

In the output file, fields in the main file will appear first, followed by the remaining
fields in xref file.

=cut

sub merge_tabfiles {
	my ($main, $fmain_str, $xref, $fxref_str, $outfile) = @_;
	
	my $fmain;
	if (ref $fmain_str eq 'HASH') {
		$fmain = $fmain_str;
	}
	else {
		$fmain = parse_fstr($fmain_str, 1);
	}
	my $fxref;
	if (ref $fxref_str eq 'HASH') {
		$fxref = $fxref_str;
	}
	else {
		$fxref = parse_fstr($fxref_str, 1);
	}

	my @mainfields = values %$fmain;
	# Find key and value fields for xref file
	my (@fkeys, @fvals);
	foreach my $field (values %$fxref) {
		if (grep { $field eq $_ } @mainfields) {
			push @fkeys, $field;
		}
		else {
			push @fvals, $field;
		}
	}

	unless(@fkeys > 0) {
		croak "Cannot find key fields shared by two files!";
	}
	unless(@fvals > 0) {
		croak "Empty fields for values to be fetched from xref!";
	}

	# Now parse xref file to slurp the data fields.
	my %xrefdat;
	my $header = 1;
	if (all { /^\d+$/ } keys %$fxref) {
		$header = 0;
	}
	my $it = iter_file($xref, { fsep => qr/\t/, header => $header, select => $fxref });
	while(my $dat = $it->()) {
		my $key = join("\t", @{$dat}{@fkeys});
		if (defined $xrefdat{$key}) {
			croak "Data for [$key] in the external ref has already been found!";
		}
		$xrefdat{$key} = join("\t", @{$dat}{@fvals});
	}

	# Join with the first file
	my $fout = IO::File->new($outfile, "w");
	print $fout join("\t", @mainfields, @fvals), "\n";
	$it = iter_file($main, { fsep => qr/\t/, header => $header, select => $fmain });
	while(my $dat = $it->()) {
		my $key = join("\t", @{$dat}{@fkeys});
		if (defined $xrefdat{$key}) {
			print $fout join("\t", @{$dat}{@mainfields}, $xrefdat{$key}), "\n";
		}
		else {
			print $fout join("\t", @{$dat}{@mainfields}, ('.') x scalar(@fvals)), "\n";
		}
	}
	# Return fields used to merge two files.
	return join(',', @fkeys);
}

=head2 fam_rels PEDFILE [,TWINS]

Find all relatives for each individual from a PED file, and give them human readable 
relationship description.

For first- or second- degree relatives, the actual relationship will be used.
All other related samples are connected by relationship.

Return a 2-d hash "fam_rels" { IID1 => {Samp1 => Rel, Samp2 => Rel}, IID2... }
For each individual, this hash stores familial relationship with all his/her family members
who appear in the 2nd colum of the PED file.

In array context, also return "fam_samp", { FID1 => [SAMP1, SAMP2 ...], FID2 ... }
fam_samp will include ONLY samples that appear in the second column of PED file.
Parents specified in the third or fourth columns may not all appeared in the second column.

Options:

* twins: a list of known twin pairs, should be an array ref. Twin pairs should must
share the same pair of parents.
* strict: under strict mode, both parents for each nonfounder must be specified. 
	Note: previously we require parents also appear in the PED file, now this is no longer required.
	Note: full vs half siblings will be inferred purely based on non-missing parents in the pedigree file.
* ignore: regex pattern for duplicate sample.
	Note: duplcate sample will be included as type of relation, but they will not 
	appear in the final familial sample list.
* select: a subset of selected families.
* shorten: shorten some (mostly 2nd degree) relationships.
* fullsib: add Full prefix to full sibs to distinguish those unknown Full or Half sibs.
* verbose: print out warning messages.


=cut

sub fam_rels {
	my ($pedfile, $argref) = @_;
	my $arg = merge_opts($argref, twins => undef, shorten => 0, strict => 0, fullsib => 0, 
		ignore => qr/_Re(\d*)$/, select => undef, verbose => 1);
	
	my $pattern;
	if (defined $arg->{ignore}) {
		if (ref $arg->{ignore} eq 'Regexp') {
			$pattern = $arg->{ignore};
		}
		else {
			my $str = $arg->{ignore};
			$str =~ s/^['"]//; $str =~ s/['"]$//;
			$pattern = qr/$arg->{ignore}/;
		}
	}

	my %famid = map { (split)[1,0] } slurp $pedfile;
	my %sex = map { (split)[1,4] } slurp $pedfile;

	my %keep;
	if (defined $arg->{select}) {
		unless(ref $arg->{select} eq 'ARRAY') {
			die "Selected family IDs must be stored in an array ref";
		}
		foreach my $fid (@{$arg->{select}}) {
			$keep{$fid} = 1;
		}
	}

	# Family IDs for all sample IDs appear in the 2nd,3rd,4th columns of the PED file
	my %knownid;
	open my $fin, $pedfile or die "Cannot open PED file";
	while(<$fin>) {
		my @samps = (split)[0..3];
		my $fid = shift @samps;
		foreach my $samp (grep { $_ ne '0' } @samps) {
			chk_default(\%knownid, $samp, $fid);
		}
	}

	my (%famsamp, %graph, %dad, %mom, %relct, %rels);
	# Note:
	# %famsamp: FamID => [list of samples in Family FamID]
	# %graph: FamID => relationship graph for Family FamID 
	# %dad, %mom will be ordinary hash: IID => IID's father/mother
	# %rels will be bivariate hash
	# $rels{ID1,ID2} = ID2's relation to ID1
	# $rels{ID1}{Type of relatinoship} = count
	open $fin, $pedfile or die "Cannot open PED file";
	while(<$fin>) {
		my ($fid, $iid, $dad, $mom, $sex) = (split)[0,1,2,3,4];
		if (defined $arg->{select}) {
			next unless defined $keep{$fid};
		}

		if ($arg->{strict}) {
			if ($dad ne '0' && $mom eq '0' || $dad eq '0' && $mom ne '0') {
				die "Non-founders must have both parents specified under strict mode!";
			}
		}

		push @{$famsamp{$fid}}, $iid;
		unless (defined $graph{$fid}) {
			$graph{$fid} = Graph::Undirected->new();
		}
		unless($sex eq '1' || $sex eq '2') {
			warn "Unknown sex for $iid" if $arg->{verbose};
		}
		
		if (defined $pattern && $iid =~ /$pattern/) {
			(my $origid = $iid) =~ s/$pattern//;
			if (defined $famid{$origid}) {
				croak "Sample $origid and duplicate $iid do not have the same FamID and sex"
					unless $famid{$origid} eq $fid && $sex{$origid} eq $sex;
				$relct{$origid}{Rep} ++;
				if ($relct{$origid}{Rep} > 1) {
					$rels{$origid,$iid} = ordinate($relct{$origid}{Rep}).'Rep';
				}
				else {
					$rels{$origid,$iid} = 'Rep';
				}
				$graph{$fid}->add_edge($origid, $iid);
				# We will skip the remaining part of this loop
				# This is to avoid counting duplicates as siblings since they share both parents
				next;
			}
			else {
				if ($arg->{verbose}) {
					warn "The original sample $origid for duplicate $iid does not exist in PED file";
				}
			}
		}

		if ($dad ne '0') {
			#if (!$arg->{strict} || defined $famid{$dad}) {
				if (defined $famid{$dad}) {
					croak "Incorrect sex for dad $dad" unless $sex{$dad} eq '1';
					croak "Dad $dad is not in the same family $fid as child $iid" unless $famid{$dad} eq $fid;
				}	
				else {
					warn "Sample ${iid}'s father $dad does not exist in PED file" if $arg->{verbose};
				}
				croak "Sample $iid's father has already been included in the graph" if defined $dad{$iid};
				$graph{$fid}->add_edge($iid, $dad);
				$dad{$iid} = $dad;
				$rels{$iid,$dad} = "Father";
				if ($sex eq '1') {
					$relct{$dad}{Son} ++;
					if ($relct{$dad}{Son}>1) {
						$rels{$dad,$iid} = ordinate($relct{$dad}{Son}).'Son';
					}
					else {	
						$rels{$dad,$iid} = 'Son';
					}
				}
				elsif ($sex eq '2') {
					$relct{$dad}{Daughter} ++;
					if ($relct{$dad}{Daughter}>1) {
						$rels{$dad,$iid} = ordinate($relct{$dad}{Daughter}).'Daughter';
					}
					else {
						$rels{$dad,$iid} = 'Daughter';
					}
				}
				else {
					$relct{$dad}{Offspring} ++;
					if ($relct{$dad}{Offspring}>1) {
						$rels{$dad,$iid} = ordinate($relct{$dad}{Offspring}).'Offspring';
					}
					else {
						$rels{$dad,$iid} = 'Offspring';
					}
				}
				if ($mom ne '0') {
					# $graph{$fid}->add_edge($dad, $mom);
					$graph{$fid}->add_weighted_edge($dad, $mom, 1.5);
					$relct{$dad}{Wife}{$mom} = 1; # use hash because it may appear multiple times
					if (keys %{$relct{$dad}{Wife}} > 1) {
						$rels{$dad,$mom} = ordinate(scalar(keys %{$relct{$dad}{Wife}})).'Wife';
					}
					else {
						$rels{$dad,$mom} = 'Wife';
					}
				}
			#}
			#else {
			#	warn "Sample ${iid}'s father $dad does not exist in PED file" if $arg->{verbose};;
			#}
		}
		#else {
		#	warn "Sample ${iid}'s father is unknown" if $arg->{verbose};
		#}
		if ($mom ne '0') {
			#if (!$arg->{strict} || defined $famid{$mom}) {
				if (defined $famid{$mom}) {
					croak "Incorrect sex for Mom $mom" unless $sex{$mom} eq '2';
					croak "Dad $mom is not in the same family $fid as child $iid" unless $famid{$mom} eq $fid;
				}
				else {
					warn "Sample ${iid}'s mother $mom cannot be found in the PED file" if $arg->{verbose};
				}
				croak "Sample $iid's mother has already been included in the graph" if defined $mom{$iid};
				$graph{$fid}->add_edge($iid, $mom);
				$mom{$iid} = $mom;
				$rels{$iid,$mom} = "Mother";
				if ($sex eq '1') {
					$relct{$mom}{Son} ++;
					if ($relct{$mom}{Son}>1) {
						$rels{$mom,$iid} = ordinate($relct{$mom}{Son}).'Son';
					}
					else {
						$rels{$mom,$iid} = 'Son';
					}
				}
				elsif ($sex eq '2') {
					$relct{$mom}{Daughter} ++;
					if ($relct{$mom}{Daughter}>1) {
						$rels{$mom,$iid} = ordinate($relct{$mom}{Daughter}).'Daughter';
					}
					else {
						$rels{$mom,$iid} = 'Daughter';
					}
				}
				else {
					$relct{$mom}{Offspring} ++;
					if ($relct{$mom}{Offspring}>1) {
						$rels{$mom,$iid} = ordinate($relct{$mom}{Offspring}).'Offspring';
					}
					else {
						$rels{$mom,$iid} = 'Offspring';
					}
				}
				if ($dad ne '0') {
					# $graph{$fid}->add_edge($mom, $dad);
					$graph{$fid}->add_weighted_edge($mom, $dad, 1.5);
					$relct{$mom}{Husband}{$dad} = 1; # use hash because it may appear multiple times
					if (keys %{$relct{$mom}{Husband}} > 1) {
						$rels{$mom,$dad} = ordinate(scalar(keys %{$relct{$mom}{Husband}})).'Husband';
					}
					else {
						$rels{$mom,$dad} = 'Husband';
					}
				}
			#}
			#else {
			#	warn "Sample ${iid}'s mother $mom cannot be found in the PED file" if $arg->{verbose};
			#}
		}
		#else {
		#	warn "Sample ${iid}'s mother is unknown" if $arg->{verbose};
		#}
	}

	my %twins;
	if (defined $arg->{twins}) {
		my @twinpairs;
		if (ref $arg->{twins} eq 'ARRAY') {
			@twinpairs = @{$arg->{twins}};
		}
		elsif (-f $arg->{twins}) {
			@twinpairs = map { [(split)] } slurp $arg->{twins};
		}
		else {
			die "Cannot recognize twins data type: $arg->{twins}";
		}
		foreach my $pair (@twinpairs) {
			if (@$pair == 2) {
				unless(defined $famid{$pair->[0]} && defined $famid{$pair->[1]}) {
					#die "Cannot find family ID for the twin pair $pair->[0] and $pair->[1]"
					warn "Not all twins can be found in the pedigree" if $arg->{verbose};
					next;
				}
				else {
					unless($famid{$pair->[0]} eq $famid{$pair->[1]}) {
						die "The twin pair $pair->[0] and $pair->[1] are not in the same family";
					}
					unless($sex{$pair->[0]} eq $sex{$pair->[1]}) {
						die "The twin pair $pair->[0] and $pair->[1] do not have the same sex!";
					}
				}
			}
			else {
				print Dumper $pair;
				die "Twin pair does not have the correct number of samples";
			}

			if (defined $dad{$pair->[0]} && defined $dad{$pair->[1]} && $dad{$pair->[0]} ne $dad{$pair->[1]}) {
				die "The twin pair $pair->[0] and $pair->[1] do not share the same father";
			}
			else {
				if (!defined $dad{$pair->[0]} &&  defined $dad{$pair->[1]} || 
					 defined $dad{$pair->[0]} && !defined $dad{$pair->[1]} ) {
					die "Father for one of twin pair $pair->[0] or $pair->[1] is missing";
				}
			}

			if (defined $mom{$pair->[0]} && defined $mom{$pair->[1]} && $mom{$pair->[0]} ne $mom{$pair->[1]}) {
				die "The twin pair $pair->[0] and $pair->[1] do not share the same mother";
			}
			else {
				if (!defined $mom{$pair->[0]} &&  defined $mom{$pair->[1]} ||
					 defined $mom{$pair->[0]} && !defined $mom{$pair->[1]}) {
					die "Mother for one of twin pair $pair->[0] or $pair->[1] is missing";
				}
			}
			$twins{$pair->[0],$pair->[1]} = 1;
			$twins{$pair->[1],$pair->[0]} = 1;
		}
	}

	# Find all sibs for each sample
	# For sibling that share one parent and the other parent unknown,
	# they will not be considered as half brother/sister
	# It is recommended to fill in missing persons in the pedigree to disambiguate such cases
	while(my ($iid, $fid) = each %famid) {
		if (defined $arg->{select}) {
			next unless defined $keep{$fid};
		}
		foreach my $samp (@{$famsamp{$fid}}) {
			next if $iid eq $samp;
			# Prefix for sib type
			my $prefix;
			if (defined $mom{$iid} && defined $mom{$samp} && $mom{$iid} eq $mom{$samp}) {
				if (defined $dad{$iid} && defined $dad{$samp}) {
					if ($dad{$iid} eq $dad{$samp}) {
						if (defined $twins{$iid,$samp}) {
							$prefix = 'Twin';
						}
						else {
							if ($arg->{fullsib}) {
								$prefix = 'Full';
							}
							else {
								$prefix = '';
							}
						}	
					}
					else {
						$prefix = 'Half';
					}
				}
				elsif (!defined $dad{$iid} && !defined $dad{$samp}) {
					$prefix = '';
				}
				else {
					$prefix = 'Half';
				}
			}
			elsif (defined $dad{$iid} && defined $dad{$samp} && $dad{$iid} eq $dad{$samp}) {
				if (!defined $mom{$iid} && !defined $mom{$samp}) {
					$prefix = '';
				}
				else {
					$prefix = 'Half';
				}
			}
			# If prefix is not defined, then it is not a sib pair
			if (defined $prefix) {
				unless($graph{$fid}->has_edge($iid, $samp)) {
					if ($prefix eq 'Half') {
						$graph{$fid}->add_weighted_edge($iid, $samp, 1.5);
					}
					else {
						$graph{$fid}->add_edge($iid, $samp);	
					}
				}	
				if ($sex{$samp} eq '1') {
					$relct{$iid}{$prefix."Brother"} ++;
					if ($relct{$iid}{$prefix."Brother"}>1) {
						$rels{$iid,$samp} = ordinate($relct{$iid}{$prefix."Brother"}).$prefix."Brother";
					}
					else {
						$rels{$iid,$samp} = $prefix."Brother";
					}
				}
				elsif ($sex{$samp} eq '2') {
					$relct{$iid}{$prefix."Sister"} ++;
					if ($relct{$iid}{$prefix."Sister"}>1) {
						$rels{$iid,$samp} = ordinate($relct{$iid}{$prefix."Sister"}).$prefix."Sister";
					}
					else {
						$rels{$iid,$samp} = $prefix."Sister";
					}
				}
				else {
					$relct{$iid}{$prefix."Sib"} ++;
					if ($relct{$iid}{$prefix."Sib"}>1) {
						$rels{$iid,$samp} = ordinate($relct{$iid}{$prefix."Sib"}).$prefix."Sib";
					}
					else {
						$rels{$iid,$samp} = $prefix."Sib";
					}
				}
			}
		}
	}
	#print Dumper \%rels;
	#print Dumper \%relct;

	if ($arg->{shorten}) {
		foreach my $fid (keys %famsamp) {
			if (defined $arg->{select}) {
				next unless defined $keep{$fid};
			}

			my @samps;
			foreach my $iid (@{$famsamp{$fid}}) {
				if (defined $pattern && $iid =~ /$pattern/) {
					(my $origid = $iid) =~ s/$pattern//;
					next if defined $famid{$origid};
				}
				push @samps, $iid;
			}
			my $gfam = dclone($graph{$fid});
			foreach my $iid (@samps) {
				foreach my $samp (@samps) {
					next if $iid eq $samp;
					next if defined $rels{$iid,$samp};
					# Finding shortest path between two samples using Dijkstra algorithm
					# Note: when multiple shortest paths exist, one will be randaomly chosen.
					# How to break the tie??
					my @path = $graph{$fid}->SP_Dijkstra($iid, $samp);
					if (@path) {
						my @roles;
						for(my $ii = 1; $ii < @path; $ii ++) {
							croak "Cannot find 1st-deg rel between $path[$ii-1] and $path[$ii]"
							unless defined $rels{$path[$ii-1],$path[$ii]} ;
							push @roles, $rels{$path[$ii-1],$path[$ii]};
						}
						my $relation = join("'s", @roles);
						my $shortrel = _shorten_relation($relation);
						if ($shortrel ne $relation) {
							$gfam->add_weighted_edge($iid, $samp, 1.5);
							$relct{$iid}{$shortrel} ++;
							if ($relct{$iid}{$shortrel}>1) {
								$rels{$iid,$samp} = ordinate($relct{$iid}{$shortrel}).$shortrel;
							}
							else {
								$rels{$iid,$samp} = $shortrel;
							}
						}
					}
				}
			}
			$graph{$fid} = $gfam;
		}
	}
	

	# Now label all related samples in the family
	my %famrels;
	while(my ($iid, $fid) = each %famid) {
		if (defined $arg->{select}) {
			next unless defined $keep{$fid};
		}

		# Ignore samples whose ID matches duplicated ID and original sample exist
		if (defined $pattern && $iid =~ /$pattern/) {
			(my $origid = $iid) =~ s/$pattern//;
			next if defined $famid{$origid};
		}

		foreach my $samp (@{$famsamp{$fid}}) {
			next if $iid eq $samp;
			if (defined $rels{$iid,$samp}) {
				# Direct relationship exists
				$famrels{$iid}{$samp} = $rels{$iid,$samp};
			}
			else {
				# If a pair does not have directly labeled relationship
				# we will look for connected relations from graph
				my @path = $graph{$fid}->SP_Dijkstra($iid, $samp);
				if (@path) {
					my @relation;
					for(my $ii = 1; $ii < @path; $ii ++) {
						unless (defined $rels{$path[$ii-1],$path[$ii]}) {
							print STDERR join("-", @path), "\n";
							print STDERR $rels{$path[$ii],$path[$ii-1]}, "\n";
							print STDERR $rels{$path[$ii-1],$path[$ii]}, "\n";
							croak "Cannot find relationship between $path[$ii-1] and $path[$ii]";
						}
						push @relation, $rels{$path[$ii-1],$path[$ii]};
					}
					#$famrels{$iid}{$samp} = join("'s", $iid, @relation);
					$famrels{$iid}{$samp} = join("'s", @relation);
				}
				else {
					# If a pair is not connected in the relationship graph
					# their relationship will be set to unknown
					if ($arg->{verbose}) {
						warn "Cannot determine the relationship between $iid and $samp";
					}
					$famrels{$iid}{$samp} = 'Unknown';
				}
			}
		}
	}
	if (wantarray) {
		return (\%famsamp, \%famrels);
	}
	else {
		return \%famrels;
	}
}

# Ordinate will be omitted after shorning the relationship
# Ordinate will be dealt with by the second pass to graph paths
sub _shorten_relation {
	my ($role) = @_;
	my $Ord = qr/\d*1st|\d*2nd|\d*3rd|\d*[^12]th/;
	my %Pref = (Father => 'Pat', Mother => 'Mat');
	my %Opps = (Father => 'Mother', Mother => 'Father', Husband => "Son", Wife => "Daughter");
	my %Avul = (Brother => 'Uncle', FullBrother => 'Uncle', TwinBrother => 'Uncle', HalfBrother => 'HalfUncle', 
				Sister => 'Aunt', FullSister => 'Aunt', TwinSister => 'Aunt', HalfSister => 'HalfAunt',  
				Son => 'Nephew', Daughter => 'Niece');
	my %Half = (HalfBrother => 'Half', HalfSister => 'Half', Brother => "", Sister => "", 
				FullBrother => "", FullSister => "", TwinBrother => "", TwinSister => "");

	$role =~ s/^($Ord)?(Husband|Wife)'s(Father|Mother)$/$3InLaw/;
	$role =~ s/^($Ord)?(Son|Daughter)'s($Ord)?(Husband|Wife)$/$Opps{$4}InLaw/;

	$role =~ s/^(Father|Mother)'s(Father|Mother)$/$Pref{$1}Grand$2/;
	$role =~ s/^($Ord)?(Son|Daughter)'s($Ord)?(Son|Daughter)/Grand$4/;

	$role =~ s/^(Father|Mother)'s($Ord)?((Twin|Full|Half)?(Brother|Sister))$/$Pref{$1}$Avul{$3}/;
	#$role =~ s/^($Ord)?Half(Brother|Sister)'s($Ord)?(Son|Daughter)$/Half$Avul{$4}/;
	$role =~ s/^($Ord)?((Twin|Full|Half)?(Brother|Sister))'s($Ord)?(Son|Daughter)$/$Half{$2}$Avul{$6}/;

	return $role;
}


# Shorten the relationship type to Dup/PO/FS/2nd/Other.
# Siblings with unknown FS or HS status, will be classified into a separate class. 
# To get more accurate types of relationships, use fam_kins below to estimate kinship coeff
sub reltype {
	my ($relation) = @_;
	my $Ord = qr/\d*1st|\d*2nd|\d*3rd|\d*[^12]th/;
	my $Rep = qr/'s($Ord)?Rep/;
	if ($relation =~ /^($Ord)?Rep|Twin(Brother|Sister)?($Rep)?$/) {
		return "Dup/MZ";
	}
	elsif ($relation =~ /^($Ord)?(Father|Mother|Son|Daughter)($Rep)?$/) {
		return "PO";
	}
	elsif ($relation =~ /^($Ord)?Full(Sister|Brother)($Rep)?$/) {
		return "FS";
	}
	elsif ($relation =~ /^($Ord)?Half(Sister|Brother)($Rep)?$/ ||
		   $relation =~ /^($Ord)?(Mat|Pat)(Aunt|Uncle|GrandFather|GrandMother)($Rep)?$/ ||
		   $relation =~ /^($Ord)?(Nephew|Niece|GrandSon|GrandDaughter)($Rep)?$/ ) {
		return "2nd";
	}
	elsif ($relation =~ /^($Ord)?(Sister|Brother)($Rep)?$/) {
		return "Sib";
	}
	else {
		return "Other";
	}
}


=head2 fam_kins

Calculate pairwise kinship coefficient for family members.
This is a wrapper of idcoefs, and also caluclate condensed coefficient.

Return a 2D hash, %idcoef: $idcoef{ID1}{ID2} = [kinship, ProbIBD2, ProbIBD1]
ID2 will be any one of ID1's family members that appear in the 2nd column of the pedigree file.

In array context, also return fam_samps same as fam_rels function.

Options:

* twins: a list of known MZ twin pairs, can be a list, array ref, or hashref. 
	It will be passed to read_twindups.
* strict: under strict mode, both parents must appear for nonfounders.
	Otherwise, a random ID will fill in the position for missing parent.
* ignore: regex pattern for duplicate sample.
	Duplicate will be treated as the same as MZ twin for calculating kinship.
* select: a subset of selected families for calculation
* full: output all condensed coeffs, default is to output kinship coeff and delta_7 and delta_8 
	(i.e. IBD1 and IBD2)
* wrkdir: working directory to store temp files (default will use tmpir).
* verbose: print out warning messages.

=cut

sub fam_kins {
	my ($pedfile, $argref) = @_;
	my $arg = merge_opts($argref, twins => undef, strict => undef, ignore => qr/_Re(\d*)$/, 
						select => undef, full => undef, wrkdir => undef, verbose => 1);

	my $idcoefs = which("idcoefs");
	unless($idcoefs) {
		die "Cannot find binary executable idcoefs!"
	}

	my %twindups; # store the sentinel ID for the twin group
	if (defined $arg->{twins}) {
		%twindups = read_twindups($arg->{twins});
	}

	my $pattern;
	if (defined $arg->{ignore}) {
		my $str = $arg->{ignore};
		$str =~ s/^['"]//; $str =~ s/['"]$//;
		$pattern = qr/$arg->{ignore}/;
	}

	my %famid = map { (split)[1,0] } slurp $pedfile;
	my %sex = map { (split)[1,4] } slurp $pedfile;

	my %keep;
	if (defined $arg->{select}) {
		unless(ref $arg->{select} eq 'ARRAY') {
			die "Selected family IDs must be stored in an array ref";
		}
		foreach my $fid (@{$arg->{select}}) {
			$keep{$fid} = 1;
		}
	}

	my $wrkdir;
	if (defined $arg->{wrkdir}) {
		$wrkdir = $arg->{wrkdir};
		make_path $wrkdir unless -d $wrkdir;
	}
	else {
		$wrkdir = tempdir(CLEANUP => 1);
	}

	my (%knownid, $idgen);
	open my $fin, $pedfile or die "Cannot open PED file";
	while(<$fin>) {
		my @a = split;
		for(my $ii = 1; $ii <= 3; $ii ++) {
			if ($a[$ii] ne '0') {
				#$knownid{$a[$ii]} = 1;
				chk_default(\%knownid, $a[$ii], $a[0]);
			}
		}
	}
	unless($arg->{strict}) {
		$idgen = String::Random->new(); 
	}

	# Create ped graphs, for twin/dups, only one node will be added to the graph
	my (%pedgraph, %randid, %pars, %renamed, %famsamp);
	open $fin, $pedfile or die "Cannot open PED file";
	while(<$fin>) {
		my ($fid, $iid, $dad, $mom, $sex) = (split)[0,1,2,3,4];
		if (defined $arg->{select}) {
			next unless defined $keep{$fid};
		}

		push @{$famsamp{$fid}} => $iid;
		unless(defined $pedgraph{$fid}) {
			# we will use directed graph, PO relation will be P->O
			$pedgraph{$fid} = Graph::Directed->new();
		}

		# Note in fam_rels function, we will also check that twins shared the same parents
		# here we only check sex and family ID
		if (defined $pattern && $iid =~ /$pattern/) {
			my $oldid = $iid;
			$iid =~ s/$pattern//;
			if (defined $famid{$oldid}) {
				croak "Sample $iid and duplicate $oldid do not have same FamID and sex!"
					unless $famid{$iid} eq $fid && $sex{$iid} eq $sex;
			}
			else {
				if ($arg->{verbose}) {
					warn "The original sample $iid for duplicate $oldid does not exist in PED file";
				}
			}
			# To account for more than one individual in the twin/dup group
			$renamed{$iid}{$oldid} = 1;
		}
		if (defined $twindups{$iid}) {
			my $oldid = $iid;
			$iid = $twindups{$iid};
			if (defined $famid{$oldid}) {
				croak "Sample $iid and his/her twin pair $oldid do not have same FamID and sex!"
					unless $famid{$iid} eq $fid && $sex{$iid} eq $sex;
			}
			else {
				if ($arg->{verbose}) {
					warn "The sentinel sample $iid for one of twin pairs $oldid does not exist in PED file";
				}
			}
			$renamed{$iid}{$oldid} = 1;
		}

		if ($dad ne '0') {
			if (defined $famid{$dad}) {
				die "Incorrect sex for dad $dad" unless $sex{$dad} eq '1';
				die "Dad $dad is not in the same family $fid as child $iid" unless $famid{$dad} eq $fid;
			}
		}
		if ($mom ne '0') {
			if (defined $famid{$mom}) {
				die "Incorrect sex for dad $mom" unless $sex{$mom} eq '2';
				die "Dad $mom is not in the same family $fid as child $iid" unless $famid{$mom} eq $fid;
			}
		}

		if ($dad eq '0' && $mom eq '0') {
			$pedgraph{$fid}->add_vertex($iid);
		}
		elsif ($dad ne '0' && $mom ne '0') {
			unless($pedgraph{$fid}->has_edge($dad, $iid)) {
				$pedgraph{$fid}->add_edge($dad, $iid);
				$pars{$iid}{DAD} = $dad;
			}
			unless($pedgraph{$fid}->has_edge($mom, $iid)) {
				$pedgraph{$fid}->add_edge($mom, $iid);
				$pars{$iid}{MOM} = $mom;
			}
		}
		else {
			if ($arg->{strict}) {
				die "Non-founders must have both parents specified under strict mode!";
			}
			else {
				if ($arg->{verbose}) {
					warn "One of parents for $iid is missing";	
				}
				# generate random IDs
				if ($dad eq '0') {
					$dad = $idgen->randregex('[A-Z]{10}');
					while(defined $knownid{$dad} || defined $randid{$dad}) {
						$dad = $idgen->randregex('[A-Z]{10}');
					}
					$randid{$dad} = 1;
				}
				unless($pedgraph{$fid}->has_edge($dad, $iid)) {
					$pedgraph{$fid}->add_edge($dad, $iid);
					$pars{$iid}{DAD} = $dad;
				}
				if ($mom eq '0') {
					$mom = $idgen->randregex('[A-Z]{10}');
					while(defined $knownid{$mom} || defined $randid{$mom}) {
						$dad = $idgen->randregex('[A-Z]{10}');
					}
					$randid{$mom} = 1;
				}
				unless($pedgraph{$fid}->has_edge($mom, $iid)) {
					$pedgraph{$fid}->add_edge($mom, $iid);
					$pars{$iid}{MOM} = $mom;
				}
			}
		}
	}

	# Condensed coefficients, kinship coeff will be calculated from 
	# delat1 + 0.5(delat3+delat5+delat7) + 0.25delta8
	my %coeffs;
	# Run idcoefs to calculate kinship coeffs
	foreach my $fid (sort keys %pedgraph) {
		my @sorted = $pedgraph{$fid}->topological_sort();
		# Restrict to families with at least 2 unique samples
		# If a family only contains multiple twin/dup pairs, they may be excluded at this step
		# So we should add them back here
		if (@sorted == 1 && defined $renamed{$sorted[0]} && defined $famid{$sorted[0]}) {
			# First get unique ID
			my $iid = $sorted[0]; 
			foreach my $alias (grep { defined $famid{$_} } keys %{$renamed{$iid}}) {
				$coeffs{$iid}{$alias} = [0.5, 1, 0];
				$coeffs{$alias}{$iid} = [0.5, 1, 0];
			}
		}
		next unless @sorted >= 2;
		# Map ped graph IDs to numbers for Idcoef calculation
		my (%id2num, %num2id);
		for(my $ii = 1; $ii <= @sorted; $ii ++) {
			$id2num{$sorted[$ii-1]} = $ii;
			$num2id{$ii} = $sorted[$ii-1];
		}
		# Write out idcoefs input files
		open my $fped, ">$wrkdir/$fid.pedigree" or die "Cannot write to $fid.pedigree";
		foreach my $iid (@sorted) {
			my $dad = $pars{$iid}{DAD};
			my $mom = $pars{$iid}{MOM};
			if (!defined $dad && !defined $mom) {
				print $fped join("\t", $id2num{$iid}, 0, 0), "\n";
			}
			elsif (defined $dad && defined $mom) {
				unless(defined $id2num{$dad} && defined $id2num{$mom}) {
					die "Cannot find $dad or $mom in family $fid";
				}
				print $fped join("\t", $id2num{$iid}, $id2num{$dad}, $id2num{$mom}), "\n";
			}
			else {
				die "Incomplete parent for nonfounder $iid!"
			}
		}
		close $fped;
		open my $fsamp, ">$wrkdir/$fid.samps" or die "Cannot write to $wrkdir/$fid.samps";
		foreach my $iid (@sorted) {
			if (defined $famid{$iid} || defined $renamed{$iid}) {
				print $fsamp $id2num{$iid}, "\n";
			}
		}
		close $fsamp;
		system(qq|$idcoefs -p $wrkdir/$fid.pedigree -s $wrkdir/$fid.samps -o $wrkdir/$fid.output > $wrkdir/$fid.log|);
		# Then collect results
		open my $fin, "$wrkdir/$fid.output" or die "Cannot open $wrkdir/$fid.output";
		while(<$fin>) {
			my @delta = split;
			my $num1 = shift @delta; my $iid1 = $num2id{$num1} // do { die "Cannot find IID from $num1"} ;
			my $num2 = shift @delta; my $iid2 = $num2id{$num2} // do { die "Cannot find IID from $num2" };
			# next if $num1 == $num2;
			if ($num1 == $num2 && defined $renamed{$iid1}) {
				my @iids = ($iid1);
				push @iids, grep { defined $famid{$_} } sort keys %{$renamed{$iid1}};
				#foreach my $alias (grep { defined $famid{$_} } keys) {
				#	$coeffs{$iid1}{$alias} = [0.5, 1, 0];
				#	$coeffs{$alias}{$iid1} = [0.5, 1, 0];
				#}
				foreach my $pair (all_pairs(@iids)) {
					$coeffs{$pair->[0]}{$pair->[1]} = [0.5, 1, 0];
					$coeffs{$pair->[1]}{$pair->[0]} = [0.5, 1, 0];
				}
				next;
			}
			next if $num1 == $num2;

			unshift @delta, undef;
			$delta[0] = $delta[1] + 0.5*($delta[3]+$delta[5]+$delta[7]) + 0.25*$delta[8];

			my $kins;
			if ($arg->{full}) {
				$kins = @delta;
			}
			else {
				$kins = [$delta[0], $delta[7], $delta[8]];
			}

			my @iid1 = ($iid1);
			my @iid2 = ($iid2);
			if (defined $renamed{$iid1}) {
				foreach my $alias1 (grep { defined $famid{$_} } keys %{$renamed{$iid1}}) {
					push @iid1, $alias1;
					if (defined $famid{$iid1}) {
						$coeffs{$iid1}{$alias1} = [0.5, 1, 0] unless defined $coeffs{$iid1}{$alias1};
						$coeffs{$alias1}{$iid1} = [0.5, 1, 0] unless defined $coeffs{$alias1}{$iid1};	
					}
					
				}
			}
			if (defined $renamed{$iid2}) {
				foreach my $alias2 (grep { defined $famid{$_} } keys %{$renamed{$iid2}}) {
					push @iid2, $alias2;
					if (defined $famid{$iid2}) {
						$coeffs{$iid2}{$alias2} = [0.5, 1, 0] unless defined $coeffs{$iid2}{$alias2};
						$coeffs{$alias2}{$iid2} = [0.5, 1, 0] unless defined $coeffs{$alias2}{$iid2};
					}
				}
			}

			for (my $ii = 0; $ii < @iid1; $ii ++) {
				for (my $jj = 0; $jj < @iid2; $jj ++) {
					if (defined $famid{$iid1[$ii]} && defined $famid{$iid2[$jj]}) {
						unless(defined $coeffs{$iid1[$ii]}{$iid2[$jj]}) {
							$coeffs{$iid1[$ii]}{$iid2[$jj]} = $kins;
						}
						else {
							die "Kinship coeffs for pair $iid1[$ii] and $iid2[$jj] are already defined!";
						}
						unless(defined $coeffs{$iid2[$jj]}{$iid1[$ii]}) {
							$coeffs{$iid2[$jj]}{$iid1[$ii]} = $kins;
						}
						else {
							die "Kinship coeffs for pair $iid2[$ii] and $iid1[$jj] are already defined!";
						}
					}

				} 
			}
		}
	}

	if (wantarray) {
		return (\%famsamp, \%coeffs);
	}
	else {
		return \%coeffs;
	}
}


=head2 read_chrlen

Read chromsome length from seqdict file.

=cut

sub read_chrlen {
	my ($seqdict) = @_;
	tie my %chrlen, 'Tie::IxHash';
	open my $fin, $seqdict or die "Cannot open sequence dictionary";
	while(<$fin>) {
		my ($chr, $len) = (/SN:(\S+)\tLN:(\d+)/);
		next unless defined $chr && defined $len;
		$chrlen{$chr} = $len;
	}
	if (wantarray) {
		#warn "Chromosome order will not be maintained if a hash is returned!";
		return %chrlen;
	}
	else {
		return \%chrlen;
	}
}


=head2 split_chrs SEQDICT

Split ordered chromosomes into equal length groups.

=cut

sub split_chrs {
	my ($seqdict) = @_;
	my $chrlen = read_chrlen($seqdict);
	my @chrord = keys %$chrlen;
	my @chrgrp = ([$chrord[0]]);
	my $pos = 0;
	my $maxlen = max(values %$chrlen);
	my $tmplen = $chrlen->{$chrord[0]};
	foreach my $chr (@chrord[1..$#chrord]) {
		if ($tmplen + $chrlen->{$chr} <= $maxlen) {
			$tmplen += $chrlen->{$chr};
		}
		else {
			$tmplen = $chrlen->{$chr};
			$pos ++;		
		}
		push @{$chrgrp[$pos]} => $chr;
	}
	return @chrgrp;
}


=head2 comb_intervals INPUT, [MERGECOUNT]

Combine Intervals given the expected number of intervals in the merged set
Input intervals must be sorted, and chromosomes should be numeric order.

The combined interval will not cross chromosome boundary.

If MERGECOUNT is not provided, intervals will be combined per chromosome.

=cut

sub comb_intervals {
	my ($input, $merge_count) = @_;
	_check_intervals($input);

	my $fin = IO::File->new($input) or croak "Cannot open input bed file: $input";

	my $count = 0;
	my @intvs;
	my ($m_chr, $m_start, $m_end) = ("", 0, 0);
	
	my (@merged, @original);
	while(<$fin>) {
		my ($chr, $start, $end) = (split)[0,1,2];

		if ($count == 0) {
			($m_chr, $m_start, $m_end) = ($chr, $start, $end);
			push @intvs, [$chr, $start, $end];
			$count ++;
		}
		else {
			# Extend the interval
			if ($chr eq $m_chr) {
				$m_end = $end;
				push @intvs, [$chr, $start, $end];
				$count ++;
			}
			else {
				push @merged, [$m_chr, $m_start, $m_end];
				push @original, [@intvs];
				($m_chr, $m_start, $m_end) = ($chr, $start, $end);
				@intvs = [$chr, $start, $end];
				$count = 1;
			}
		}

		if (defined $merge_count && $count == $merge_count) {
			push @merged, [$m_chr, $m_start, $m_end];
			push @original, [@intvs];
			($m_chr, $m_start, $m_end) = ("", 0, 0);
			$count = 0;
			@intvs = ();
		}
	}
	if (@intvs) {
		push @merged, [$m_chr, $m_start, $m_end];
		push @original, [@intvs];
	}

	return (\@merged, \@original);
}


sub _check_intervals {
	my ($input) = @_;
	my $fin = IO::File->new($input) or croak "Cannot open input bed file: $input";
	my ($prev_chr, $prev_start, $prev_end) = ("", 0, 0);
	while(<$fin>) {
		my ($chr, $start, $end) = (split)[0,1,2];
		validate_elem($chr, $start, $end);
		if ($prev_chr eq $chr) {
			croak "Interval start: $start <= $prev_end" unless $start > $prev_end;
		}
		if ($prev_chr ne $chr) {
			my $prev_chr_num = hg_chr($prev_chr);
			my $chr_num = hg_chr($chr);
			croak "Chrom: $chr <= $prev_chr" unless $chr_num > $prev_chr_num;
		}
		($prev_chr, $prev_start, $prev_end) = ($chr, $start, $end);
	}
	return 1;
}



=head2 region_bk dbfile, fields { callback, probefile }

Read regions from a tab delimited file to binkeeper.
Fields should be provided as a hashref of column number or name => alias. Fields whose 
aliases are Chrom,Start,End will be used to determine genomic coordinates. All fields 
will also be packed into the additional value field. 

If an optional probe file is provided, then genomic coordinates for binkeeper will be
based on the order of probes. Probe file is a three column file containing 
non-overlap genomic intervals used as probes to backbone the intervals in dbfile. 
We will use fast binary search to look for overlapping probes with certain percetage
overlap (default: 50%) for the intervals in the dbfile. 

An optional callback subroutine can be provided to populate value fields.

=cut

sub read_probes {
	my ($probefile, $fields, $argref) = @_;
	my $arg = merge_opts($argref, bed => undef);

	my $bedflag = $arg->{bed};
	if ($probefile =~  /\.bed$/ || $probefile =~ /\.bed\.gz$/) {
		$bedflag = 1 unless defined $bedflag;
		$fields = "1:Chrom,2:Start,3:End" unless defined $fields;
	}
	elsif ($probefile =~ /\.pfb$/ || $probefile =~ /\.pfb\.gz$/) {
		$fields = "Chr:Chrom,Position:Start" unless defined $fields;
	}

	unless (ref $fields eq 'HASH') {
		$fields = parse_fstr($fields);
	}
	foreach my $alias (qw|Chrom Start|) {
		unless(grep { $_ eq $alias } values %$fields) {
			die "Cannot find standard column alias $alias from $probefile";
		}
	}
	my $header = 1;
	if (all { /^\d+$/ } keys %$fields) {
		$header = 0;
	}

	my %probes;
	my $it = iter_file($probefile, { fsep => qr/\t/, header => $header, select => $fields });
	while(my $dat = $it->()) {
		my ($chrom, $start, $end) = @{$dat}{qw|Chrom Start End|};
		$start += 1 if $bedflag;
		$end = $start unless defined $end;
		validate_elem($chrom, $start, $end);

		push @{$probes{$chrom}{Starts}}, $start;
		push @{$probes{$chrom}{Ends}}, $end;
	}

	# Then we check if intervals are ordered by position
	foreach my $chrom (sort keys %probes) {
		my $starts = $probes{$chrom}{Starts};
		my $ends = $probes{$chrom}{Ends};
		my $flag;
		for(my $ii = 1; $ii < @$starts; $ii ++) {
			if ($starts->[$ii] < $starts->[$ii-1] ||
				$ends->[$ii] < $ends->[$ii-1]) {
				$flag = 1;
				last;
			}
		}
		if ($flag) {
			print STDERR "Sorting intervals on chromosome $chrom\n";
			# Sort intervals based on their mid-point
			my @midpts = map { 0.5*($starts->[$_]+$ends->[$_]) } 0..(@$starts-1);
			my @order = sort { $midpts[$a] <=> $midpts[$b] } 0..(@$starts-1);
			# Then re-order starts and ends
			$probes{$chrom}{Starts} = [ @{$starts}[@order] ];
			$probes{$chrom}{Ends} = [ @{$ends}[@order] ];
			# Then check if there are overlaps in intervals
			# While overlapping between intervals may be a not fatal err, 
			# they invalidate the assumption that start and end are both ordered
			# after sorting on mid-points. It will die if probes are overlapping.
			$starts = $probes{$chrom}{Starts};
			$ends = $probes{$chrom}{Ends};
			for(my $ii = 1; $ii < @$starts; $ii ++) {
				unless($starts->[$ii] >= $ends->[$ii-1]) {
					die "Probe intervals on $chrom are overlapping after sorting!";
				}
			}
		}
	}
	return \%probes;
}

# Return a half-open 1-based interval of ranges
# Default overlapping percentage is 0.5 requiring 50% of the probe interval
# to be included in the overlapping region.
sub find_proberng {
	my ($probes, $chrom, $start, $end, $argref) = @_;
	my $arg = merge_opts($argref, overlap => 0.5);

	unless(defined $probes->{$chrom}) {
		return ();
	}

	unless(@{$probes->{$chrom}{Starts}} == @{$probes->{$chrom}{Ends}}) {
		die "Unequal length of starts and ends for $chrom";
	}

	# Binary search on both starts and end points
	my $low_pos = binsearch_pos { $a <=> $b } $start, @{$probes->{$chrom}{Starts}};
	my $hi_pos  = binsearch_pos { $a <=> $b } $end,   @{$probes->{$chrom}{Starts}};

	if ($hi_pos == 0) {
		return ();
	}
	else {
		$low_pos -- if $low_pos > 0;
		$hi_pos -- if $hi_pos >= scalar(@{$probes->{$chrom}{Starts}});
		#print "$chrom\t$low_pos\t$hi_pos\t$probes->{$chrom}{Ends}[$low_pos]\t$probes->{$chrom}{Ends}[$hi_pos]\n";

		# Then check each intervals for overlap at both ends
		my $overlap = range_overlap($start, $end, 
									$probes->{$chrom}{Starts}[$low_pos], $probes->{$chrom}{Ends}[$low_pos]);
		if ($overlap < $arg->{overlap}*($probes->{$chrom}{Ends}[$low_pos]-$probes->{$chrom}{Starts}[$low_pos]+1)) {
			$low_pos ++;
		}
		if ($low_pos > $hi_pos) {
			return ();
		}
		$overlap = range_overlap($start, $end, 
								 $probes->{$chrom}{Starts}[$hi_pos], $probes->{$chrom}{Ends}[$hi_pos]);
		if ($overlap < $arg->{overlap}*($probes->{$chrom}{Ends}[$hi_pos]-$probes->{$chrom}{Starts}[$hi_pos]+1)) {
			$hi_pos --;
		}
		if ($low_pos > $hi_pos) {
			return ();
		}
		else {
			# Change 0-based index to 1-based index to keep compatible with BinKeeper
			return ($low_pos+1, $hi_pos+1);
		}
	}
}

sub region_bk {
	my ($dbfile, $fields, $argref) = @_;
	my $arg = merge_opts($argref, poverlap => 0.5, callback => undef, bed => undef, probes => undef, verbose => undef );

	# If a probe file is provided
	# First store probefiles into sorted array per chromosome
	# then look for overlapping intervals from dbfile
	my $probes = $arg->{probes};
	my ($probechr, $dbchr);
	if ($arg->{probes}) {
		die "Probes must be a hashref of arrays (use read_probes)" unless ref $probes eq 'HASH';
		my $firstpbchr = peek_hash($probes);
		if ($firstpbchr =~ /^chr/) {
			$probechr = 1;
		}
		else {
			$probechr = 0;
		}
	}

	unless (ref $fields eq 'HASH') {
		$fields = parse_fstr($fields);
	}

	foreach my $alias (qw|Chrom Start End|) {
		unless(grep { $_ eq $alias } values %$fields) {
			die "Cannot find standard column alias $alias from $dbfile";
		}
	}
	my $header = 1;
	if (all { /^\d+$/ } keys %$fields) {
		$header = 0;
	}

	my $bedflag = $arg->{bed};
	if ($dbfile =~ /\.bed$/ || $dbfile =~ /\.bed\.gz$/ || 
		$dbfile =~ /\.bedgraph$/ || $dbfile =~ /\.bedgraph\.gz$/) {
		$bedflag = 1;
	}

	my $it = iter_file($dbfile, { fsep => qr/\t/, header => $header, select => $fields });
	my $bk = Genome::UCSC::BinKeeper->new();
	while(my $dat = $it->()) {
		validate_elem(@{$dat}{qw|Chrom Start End|});
		unless(defined $dbchr) {
			if ($dat->{Chrom} =~ /^chr/) {
				$dbchr = 1;
			}
			else {
				$dbchr = 0;
			}
			if(defined $probechr && $probechr != $dbchr) {
				warn "Chromsome nomenclatures in probe file and database are different" if $arg->{verbose};
			}
		}

		my $chrom = $dat->{Chrom};
		my ($start, $end);
		if ($arg->{probes}) {
			if ($probechr == 1 && $dbchr == 0) {
				$chrom = hg_chrom($chrom);
			}
			elsif ($probechr == 0 && $dbchr == 1) {
				$chrom =~ s/^chr//; $chrom = 'MT' if $chrom eq 'M';
			}
			unless (defined $probes->{$chrom}) {
				warn "Cannot find probes on chromosome $dat->{Chrom}" if $arg->{verbose};
				next;
			}
			else {
				my ($low_pos, $hi_pos);
				if ($bedflag) {
					($low_pos, $hi_pos) = find_proberng($probes, $chrom, $dat->{Start}+1, $dat->{End}, { overlap => $arg->{poverlap} });
				}
				else {
					($low_pos, $hi_pos) = find_proberng($probes, $chrom, @{$dat}{qw|Start End|}, { overlap => $arg->{poverlap} });
				}
				if (defined $low_pos) {
					($start, $end) = ($low_pos, $hi_pos);
				}
				else {
					warn "Interval $chrom:$dat->{Start}-$dat->{End} has no probe coverage" if $arg->{verbose};
					next;
				}
			}
			#print join("\t", $dat->{Start}, $dat->{End}, $low_pos, $hi_pos), "\n";
		}
		else {
			$start = $bedflag ? $dat->{Start}+1 : $dat->{Start};
			$end = $dat->{End};
		}

		if (defined $arg->{callback}) {
			$bk->add($dat->{Chrom}, $start, $end, $arg->{callback}->($dat));
		}
		else {
			$bk->add($dat->{Chrom}, $start, $end, $dat);
		}
		#print join("\t", @{$dat}{qw|Chrom Start End|}, $start, $end), "\n";
	}
	return $bk;
}


1;
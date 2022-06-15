package Utils::File::Iter;

use strict;
use warnings;
use Carp;
use Hash::Util qw|lock_hash|;
use List::Util qw|max|;
use List::MoreUtils qw|any uniq|;
use IO::Detect;
use IO::Uncompress::Gunzip;
use PerlIO::gzip;
use Text::CSV;
use Spreadsheet::Read;
use Data::Table;
use Iterator::Simple qw|iterator|;
use Data::Dumper;
use Utils::Hash qw|merge_opts|;
use Utils::File qw|xls_label2col|;

use base qw|Exporter|;
our @EXPORT = qw|iter_file slurp_file data_tab|;


=head1 NAME

Utils::File::Iter - Iterator access to data files.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

    use Utils::File::Iter;

    my $iter = iter_file($data_file); # with default parameters
    # In array context, also get field names and specify tab as separator
    my ($iter, $fnames) = iter_file($data_file, { fsep => qr/\t/ }); 

    while(<$iter>) {
        my $dat = $_; # I advise not to use $_ directly, as $_ is subject to change
        ... do things with $dat ...
    }
    while(my $dat = $iter->()) {
        ... do things with $dat ...
    }

=head1 DESCRIPTION

The main function of this module is to create an iterator (an L<Iterator::Simple> object) 
for going through each line of delimiter separated data file.

Supported file format including: plain or gzipped text, csv, tsv file and spreadsheets
including ods, xls, xlsx formats.

For small files, can also slurp the entire content, see L<Data::Table> module.

=head1 SUBROUTINES/METHODS

=head2 iter_file FILE [, OPTIONS]

C<iter_file> expects a file name (or handle) and a hash reference for options.
By default, file type is inferred from file name suffix. It can be overriden
by the options, which is useful when file handle is used for input stream.

By default, gzipped files will be parsed by L<IO::Uncompress::Gunzip> module
with MultiStream option turned on to support bgzip format. But this would
incur performance issue, L<see here|https://goo.gl/uDHn2h>. To speed up gzip file
parsing, L<PerlIO::gzip> can be used for regular gzip files. This can be achieved
by either provide file handle returned by C<PerlIO::gzip>, or provide a gzip faster
option (see below).

=head3 Options

=over 5

=item C<type>

File type (string). Use this to override the default guess based on file name.
When file handle is used, default is to assume txt.
Supported file types: F<< txt, tsv, csv, xlsx, ods, xls >>

=item C<encode>

(Only for unzipped txt/csv files) Specify the input file encoding (e.g. UTF-8).

=item C<fastgz>

Indicator to use faster gzip file parser (0/1). Default is 0.
Turn on this option only for ordinary gzip files, because it does not
supprot multistream as required for parsing bgzip. 

=item C<fsep>

(Only for txt/csv file) Field separator (regex or string). 
Default is qr/\s+/ for txt file, ',' for csv file, qr/\t/ for tsv file.
Not used when praseing spreadsheets.

=item C<fquote>

(Csv file only) The character to quote fields.

=item C<header>

Indicator for the existance of header line (0/1). Default is 1, set to 0 if no header.

=item C<fields>

Provide field names when C<header=0> (array ref). 
If not provided, field names for files with header=0 is 1,2,..ncol (1-based).
Where the number of columns (ncol) is determined based on the first line.
Note: if C<header=0> and file argument is non-seekable file handle, then field
names must be provided.

=item C<alias>

Field name alias (hashref).

=item C<subset>

Subset of the fields by names. When alias is provided, it should use alias.
TechNote: subsetting is applied after aliasing. So when table is large, it's more likely to have
name conflict when alias is applied first. To circumvent this issue, use *select*

=item C<exclude>

A list of fields to be excluded. When alias is provided, it should use alias
as the field name. When subset is provided, this option will be ignored.

=item C<select>

Provide a hashref of field names and (optionally) their aliases. The iteration of the file
will then be restricted to the fields specified by this option. alias, subset, exclude
will be ignored if select is provided.

=item C<skip>

The number of lines (rows) to skip from the beginning (counts).

=item C<ignore>

(Txt/tsv file only) Ignore the lines with the pattern (regex).
By default, empty lines are also automatically ignored. 

=item C<chomp>

(Txt only) chomp off trailing white space before paring a line (0-3). Default: 0.

   chomp  leading  trailing
   -----  -------  ---------
     0      n/a      n/a
     1     chomp     n/a
     2      n/a     chomp
     3     chomp    chomp

=item C<strip>

Remove leading and/or trailing white space from every field (0-3). Default: 0.
For txt file parsing with default C<fsep>, white spaces should be removed automatically. 

   strip  leading  strailing
   -----  -------  ---------
     0      n/a      n/a
     1     strip     n/a
     2      n/a     strip
     3     strip    strip

=item C<strict>

Exit (1) or skip (0) upon error lines. Default is 1 (recommended).
Can be turned off to allow some mismatches between the number of values and fields.

=item C<sheet>

(Spreadsheet only) Specify the data sheet by name or by number.

=item C<raw>

(Spreadsheet only) Indicator for processing raw values rather than formatted values,
default: 1.

=item C<maxcol> and C<maxrow>

(Spreadsheet only) The max range of col and row, override the default.
Customized value cannot exceed default max.

=back

=cut

my @OPTIONS = qw(type encode fastgz fsep fquote header fields select alias
				subset exclude skip ignore chomp strip strict sheet raw maxcol maxrow);

my %TYPE = (txt => 'text', tsv => 'text', csv => 'csv',
			xls => 'spreadsheet', xlsx => 'spreadsheet', ods => 'spreadsheet');

my %CALL = (text => \&iter_file_text, csv => \&iter_file_csv,
			spreadsheet => \&iter_file_spreadsheet);

# determine file type
sub _type_byname {
	my ($file) = @_;
	my $gzip = $file =~ /\.gz$/ ? 1 : 0;
	$file =~ s/\.gz$//;
	my $type;
	if ($file =~ /\.(xlsx|xlx|ods)$/) {
		$type = $1;
	}
	elsif ($file =~ /\.(csv|tsv)$/) {
		$type = $1;
	}
	else {
		$type = "txt";
	}
	return ($gzip, $type);
}

sub iter_file {
	my ($file, $argref) = @_;
	if (defined $argref) {
		foreach my $opt (keys %$argref) {
			unless(grep { $_ eq $opt } @OPTIONS) {
				croak "iter_file does not support option $opt";
			}
		}
	}
	# can only specify file type if input is filehandle
	if (is_filehandle $file) {
		if (defined $argref->{type}) {
			croak "Unsupported file type: $argref->{type}"
				unless defined $TYPE{$argref->{type}};
			$CALL{$TYPE{$argref->{type}}}->($file, $argref);
		}
		else {
			iter_file_text($file, $argref);
		}
	}
	else {
		croak "File $file does not exist!" unless -f $file;
		my ($gzip, $type) = _type_byname($file);
		# override the default type inferred from suffix???
		#if (defined $argref->{type}) {
		#	if (defined $TYPE{$argref->{type}} ) {
		#		$type = $argref->{type};
		#	}
		#	else {
		#		croak "Unsupported file type: $argref->{type}";
		#	}
		#}
		if (defined $argref->{type}) {
			unless (defined $TYPE{$argref->{type}}) {
				croak "Unsupported file type: $argref->{type}";
			}
			else {
				unless ($type eq $argref->{type}) {
					carp "Incorrect file type: $argref->{type} ne $type, force using $argref->{type}";
					$type = $argref->{type};
				}
			}
		}

		my $fh;
		if ($gzip) {
			if (defined $argref->{fastgz} && $argref->{fastgz}) {
				open $fh, "<:gzip", $file or croak "Cannot open $file"; 
			}
			else {
				$fh = IO::Uncompress::Gunzip->new($file, MultiStream => 1)
					or croak "Cannot open $file";
			}
			$CALL{$TYPE{$type}}->($fh, $argref);
		}
		else {
			if ($type eq 'txt' || $type eq 'csv' || $type eq 'tsv') {
				if (defined $argref->{encode}) {
					open $fh, "<:encoding($argref->{encode})"
						or croak "Cannot open $file";
				}
				else {
					open $fh, $file or croak "Cannot open $file";
				}
				$CALL{$TYPE{$type}}->($fh, $argref);
			}
			else {
				iter_file_spreadsheet($file, $argref);
			}
		}		
	}
}

=head2 data_tab FILE [, OPTIONS]

Create a Data::Table object from file. OPTIONS can be all C<iter_file> options.

=cut

sub data_tab {
	my ($file, $argref) = @_;
	my ($iter, $fnames) = iter_file($file, $argref);
	my $tab = Data::Table->new([], $fnames);
	while(my $dat = $iter->()) {
		$tab->addRow($dat);
	}
	return $tab;
}

=head2 slurp_file FILE [, OPTIONS] 

Slurp all the content of a file. It uses C<iter_file> to extract the content,
and store in an array of hashrefs.

=cut

sub slurp_file {
	my ($file, $argref) = @_;
	my $iter = iter_file($file, $argref);
	my @out;
	while(my $dat = $iter->()) {
		push @out, $dat;
	}
	if (wantarray) {
		return @out;
	}
	else {
		return \@out;
	}
}


# This helper function returnes a full list of renamed fields
# and a list of index of selected fields that will be used by _hash_fdval helper
# Note: in this function we do alias first, then subseting or excluding
sub _alias_subfd {
	my ($fields, $alias, $subset, $exclude) = @_;

	my @fields = @$fields;
	unless (scalar(@fields) == scalar(uniq sort @fields) ) {
		#print STDERR join("\n", @fields), "\n";
		my %count;
		$count{$_} ++ foreach @fields;
		print STDERR join(",", grep { $count{$_} > 1 } sort keys %count), "\n";
		croak "There are non-unique field names parsed from file header: _alias_subfd";
	}

	my @renamefd;
	if (defined $alias) {
		if (ref $alias eq 'HASH') {
			my %fd = map { $_ => 1 } @fields;
			my @notfound = grep { !defined $fd{$_} } keys %{$alias};
			if (@notfound) {
				print STDERR join(" ", @notfound), "\n";
				croak "Field names for aliasing cannot be found";
			}
			else {
				@renamefd = map { $alias->{$_} // $_ } @fields;
			}
		}
		else {
			croak "Alias must be stored in hashref";
		}
		unless (scalar(@renamefd) == scalar(uniq sort @renamefd)) {
			print STDERR join("\n", @renamefd), "\n";
			croak "There are non-unique field names after renaming!";
		}
	}
	else {
		@renamefd = @fields;
	}
	
	my @subind;
	if (defined $subset) {
		if (defined $exclude) {
			warn "Exclude list will be ignored when subset list is provided";
		}
		if (ref $subset eq 'ARRAY') {
			my %refd = map { $_ => 1 } @renamefd;
			my @notfound = grep { !defined $refd{$_} } @{$subset};
			if (@notfound) {
				print STDERR join(" ", @notfound), "\n";
				croak "Field names for subsetting cannot be found";
			}
			else {
				foreach my $ii (0..$#renamefd) {
					if (grep { $_ eq $renamefd[$ii] } @{$subset}) {
						push @subind, $ii;
					}
				}
			}
		}
		else {
			croak "Subset must be stored in arrayref";
		}
		unless (@subind == scalar(uniq @subind)) {
			print STDERR join(" ", @subind), "\n";
			croak "There are non-unique subset index";
		}
	}
	elsif (defined $exclude) {
		if (ref $exclude eq 'ARRAY') {
			my %refd = map { $_ => 1 } @renamefd;
			my @notfound = grep { !defined $refd{$_} } @{$exclude};
			if (@notfound) {
				print STDERR join(" ", @notfound), "\n";
				carp "Some field names in the exclusion list cannot be found";
			}
			else {
				foreach my $ii (0..$#renamefd) {
					unless (grep { $_ eq $renamefd[$ii] } @{$exclude}) {
						push @subind, $ii;
					}
				}
			}
		}
	}
	else {
		@subind = 0..$#renamefd;
	}
	return (\@renamefd, \@subind);
}

# This is similar to _alias_subfd, except that we do subseting first then alias to avoid name confliect 
# in table with many columns
sub _select_subfd {
	my ($fields, $select) = @_;

	my @fields = @$fields;
	unless (scalar(@fields) == scalar(uniq sort @fields))  {
		print STDERR join("\n", @fields), "\n";
		croak "There are non-unique field names parsed from file header: _select_subfd";
	}

	unless(ref $select eq 'HASH') {
		croak "Must provide a hash ref to select and rename fields";
	}
	else {
		foreach my $sfield (keys %$select) {
			unless (grep { $sfield eq $_ } @fields) {
				croak "Selected field $sfield cannot be found from file header";
			}
		}
		my @srename = values %$select;
		unless(scalar(@srename) == scalar(uniq sort @srename)) {
			print STDERR join(" ", @srename), "\n";
			croak "There are non-unique field names after selected rename";
		}
	}

	my (@renamefd, @subind);
	for(my $ii = 0; $ii < @fields; $ii ++) {
		if (defined $select->{$fields[$ii]}) {
			push @subind, $ii;
			push @renamefd, $select->{$fields[$ii]};
		}
		else {
			push @renamefd, $fields[$ii];
		}
	}
	my @subalias = @renamefd[@subind];
	unless(scalar(@subalias) == scalar(uniq sort @subalias)) {
		print STDERR join(" ", @subalias), "\n";
		croak "There are non-unique field aliases after selection";
	}

	return (\@renamefd, \@subind);
}

# Strip white spaces in fields
sub _strip_fd {
	my ($fields, $level) = @_;
	if ($level > 0) {
		if ($level == 1) {
  			foreach (@$fields) {  s/^\s+//; } 
  		}
  		elsif ($level == 2) {
  			foreach (@$fields) {  s/\s+$//;  }
  		}
  		elsif ($level == 3) {
  			foreach (@$fields) {  s/^\s+//; s/\s+$//; }
  		}
  		else {
  			print STDERR "strip=$level is not supported\n";
  		}
	}
}

# strip while spaces in header line
sub _chomp_line {
	my ($line, $level) = @_;
	if ($level == 1) {
		$line =~ s/^\s+//;
	}
	elsif ($level == 2) {
		$line =~ s/\s+$//;
	}
	elsif ($level == 3) {
		$line =~ s/^\s+//; $line =~ s/\s+$//;
	}
	else {
		if ($level > 0) {
			print STDERR "chomp=$level is not supported\n";
		}
	}
	return $line;
}

# Fill in values to hash 
sub _hash_fdval {
	my ($fields, $vals, $subind, $strip) = @_;

  	if (@$vals != @$fields) {
  		print STDERR join(",", @$vals), "\n";
  		return undef;
  	}

  	my @subfd = @$fields[@$subind];
  	my @subval = @$vals[@$subind];
  	if ($strip > 0) {
  		_strip_fd(\@subval, $strip);
  	}

  	my %F;
  	@F{@subfd} = @subval;
  	return \%F;
}

sub iter_file_text {
	my ($fh, $argref) = @_;
	my %arg = ( 
		fsep => qr/\s+/, 	# field separator, regex or string
		header => 1, 		# if 0 means no header, default: 1
		fields => undef, 	# provide field names if header=0, 
							# will create default field names as 1,2,...
		alias  => undef,	# rename field names, hashref
		subset => undef,	# subset of the fields
							# use alias when alias is specified
		exclude => undef,
		select => undef,
		ignore => undef,	# regex pattern for comments (skip).
		skip => undef,		# number of first few lines to skip		
		chomp => 0,         # chomp off leading or trailing white space before paring a line
		strict => 1, 		# exit upon error 
		strip => 0, 		# strip white spaces
		ref $argref eq 'HASH' ? %$argref : ()
		);
	lock_hash(%arg);

	if (defined $arg{ignore}) {
		unless (ref $arg{ignore} eq 'Regexp') {
			$arg{ignore} = qr/^$arg{comment}/;
		}
	}

	my $linect = 0;
  	my @fields; # original field names for all
  	my $firstline; # for files with no header, the first line will be stored.
  	{
  		if (defined $arg{skip}) {
  			croak "Skip must be an integer" unless $arg{skip} =~ /^\d+$/;
  			<$fh> for 1..$arg{skip};
  			$linect += $arg{skip};
  		}
  		
  		unless ($arg{header}) {
  			if (defined $arg{fields}) {
  				if (ref $arg{fields} eq 'ARRAY') {
  					# for text format, we do not check nfields = ncol
  					# error will be reported later when #val <> #fields 
  					@fields = @{$arg{fields}};
  				}
  				else {
  					croak "Custom field names must be stored in array ref";
  				}
  			}
  			else {
  				#$_ = <$fh>;
  				#@fields = 1..scalar(split $arg{fsep});
  				if (defined $arg{ignore}) {
  					while(<$fh>) {
  						$linect ++;
  						last unless $_ =~ $arg{ignore};
  					}
  				}
  				else {
  					$_ = <$fh>; 
  				}
	  			chomp();
  				$_ = _chomp_line($_, $arg{chomp});

  				$firstline = $_;
  				@fields = 1..scalar(split($arg{fsep}, $firstline));
  				# Note: seek method here cannot go backward for compressed file or STDIN!
  				# Currently we only error out if input file handle is IO::Uncompress::Gunzip
  				# There should be a better solution
  				#if (ref $fh eq 'IO::Uncompress::Gunzip') {
  				#	die "File handle is non-seekable, field names must be provided!";
  				#}
  				#seek($fh, -length($_), 1);	
  			}
  		}
  		else {
  			if (defined $arg{ignore}) {
  				while(<$fh>) {
  					$linect ++;
  					last unless $_ =~ $arg{ignore};
  				}
  			}
  			else {
  				$_ = <$fh>; $linect ++;
  			}
  			chomp();

  			# header line: remove leading and trailing spaces
  			$_ = _chomp_line($_, $arg{chomp});
  			#$_ =~ s/^\s+//; $_ =~ s/\s+$//;

  			@fields = split $arg{fsep};
			if($arg{strip} > 0){
  				_strip_fd(\@fields, $arg{strip});
  			}
  		}
  	}
  	my ($renamefd, $subind);
  	if (defined $arg{select}) {
  		if(any { defined $arg{$_} } qw|alias subset exclude|) {
  			warn "Options alias/subset/exclude will be ignored when select is provided"
  		}
  	 	($renamefd, $subind) = _select_subfd(\@fields, $arg{select});
  	}
  	else {
  		($renamefd, $subind) = _alias_subfd(\@fields, $arg{alias}, $arg{subset}, $arg{exclude});
  	}
	
  	my $iter = iterator {
  		if (defined $firstline) {
  			# first line should already been chomp'ed
  			my @vals = split $arg{fsep}, $firstline, -1;
  			my $F = _hash_fdval($renamefd, \@vals, $subind, $arg{strip});
  			unless(defined $F) {
  				croak "Cannot parse the first line in a file without header!";
  			}
  			$firstline = undef;
  			return $F;
  		}

  		local $_;
  		while(<$fh>) {
  			$linect ++;
  			next if /^\s*$/;

  			if (defined $arg{ignore}) {
  				next if /$arg{ignore}/;
  			}

  			chomp();

  			$_ = _chomp_line($_, $arg{chomp});
  			# remove leading space
  			# $_ =~ s/^\s+//;

  			# Note: split will remove trailing empty fields if limit is not provided
  			# However, by specifying the number of columns as the number of fields
  			# the sanity check will be automatically by-passed.
  			# So we use a negative limit to get all trailing empty fields.
  			# The value of empty field will have an empty string but not undef!
  			#my @vals = split $arg{fsep}, $_, scalar(@fields);
  			my @vals = split $arg{fsep}, $_, -1;

  			my $F = _hash_fdval($renamefd, \@vals, $subind, $arg{strip});
	  		unless ( defined $F ) {
	  			print STDERR "At line $linect: ";
	  			if ($arg{strict}) {
	  				croak "Number of values is not equal to fields";
	  			}
	  			else {
	  				print STDERR "Number of values is not equal to fields\n";
	  				next;
	  			}
	  		}
	  		else {
	  			return $F;
	  		}
  		}
  		return;
  	};

  	# Produce iterator
  	if(wantarray) {
  		return ($iter, [@{$renamefd}[@$subind]]);
  	}
  	else {
  		return $iter;
  	}
}

sub iter_file_csv {
	my ($fh, $argref) = @_;
	my %arg = ( 
		fsep => q|,|, 		# field separator, regex or string
		fquote => q|"|,		# quotation mark
		header => 1, 		# if 0 means no header, default: 1
		fields => undef, 	# provide field names if header=0, 
							# will create default field names as 1,2,...
		alias  => undef,	# rename field names, hashref
		subset => undef,	# subset of the fields
							# use alias when alias is specified
		exclude => undef,
		select => undef,
		skip => undef,		# number of first few lines to skip		
		strict => 1, 		# exit upon error 
		strip => 0, 		# strip white spaces
		ref $argref eq 'HASH' ? %$argref : ()
		);
	lock_hash(%arg);

	my $csv = Text::CSV->new ( { binary => 1, sep => $arg{fsep}, quote => $arg{fquote} } )
		or croak "Cannot use CSV: ".Text::CSV->error_diag();

	my $linect = 0;
	my @fields;
	{
		if (defined $arg{skip}) {
			croak "Skip must be an integer" unless $arg{skip} =~ /^\d+$/;
			<$fh> for 1..$arg{skip};
			$linect += $arg{skip};
		}

		my $firstline = <$fh>; 
		$linect ++;
		my $status = $csv->parse($firstline);
		croak "Cannot parse first line" unless $status;
 		my @columns = $csv->fields(); 
 		if($arg{strip} > 0){
  			_strip_fd(\@columns, $arg{strip});
  		}

 		unless ($arg{header}) {
 			if (defined $arg{fields}) {
 				if (ref $arg{fields} eq 'ARRAY') {
 					croak "Number of fields is not equal to ncol"
 						unless @{$arg{fields}} == @columns;
 					@fields = @{$arg{fields}};
 				}
 				else {
 					croak "Custom field names must be stored in array ref";
 				}
 			}
 			else {
 				@fields = 1..@columns;
 			}
 			seek($fh, -length($firstline), 1);
 			$linect --;
 		}
 		else {
 			@fields = @columns;
 		}
 	}
 	my ($renamefd, $subind);
  	if (defined $arg{select}) {
  		if(any { defined $arg{$_} } qw|alias subset exclude|) {
  			warn "Options alias/subset/exclude will be ignored when select is provided"
  		}
  	 	($renamefd, $subind) = _select_subfd(\@fields, $arg{select});
  	}
  	else {
  		($renamefd, $subind) = _alias_subfd(\@fields, $arg{alias}, $arg{subset}, $arg{exclude});
  	}

 	my $iter = iterator {
		while ( my $row = $csv->getline( $fh ) ) {
			$linect ++;
			my $F = _hash_fdval($renamefd, $row, $subind, $arg{strip});
			unless ( defined $F ) {
	  			print STDERR "At line: $linect ";
	  			if ($arg{strict}) {
	  				croak "Number of values is not equal to fields";
	  			}
	  			else {
	  				print STDERR "Number of values is not equal to fields\n";
	  				next;
	  			}
	  		}
	  		else {
	  			return $F;
	  		}
		}
		return;
	};

	if (wantarray) {
  		return($iter, [@{$renamefd}[@$subind]]);
  	}
  	else {
  		return $iter;
  	}
}

sub iter_file_spreadsheet {
	my ($fh, $argref) = @_;
	my %arg = ( 
		type => undef,		# file type, xls/ods/xlsx [required]
		sheet => undef,		# sheet label or name [required]
		raw => 1, 			# raw (1) or formatted (0) value.
		header => 1, 		# if 0 means no header, default: 1
		fields => undef, 	# provide field names if header=0, 
							# will create default field names as 1,2,...
		alias  => undef,	# rename field names, hashref
		subset => undef,	# subset of the fields
							# use alias when alias is specified
		exclude => undef,
		select => undef,
		skip => undef,		# number of first few lines to skip		
		strict => 1, 		# exit upon error 
		strip => 0, 		# strip white spaces
		maxcol => undef, maxrow => undef, 
		ref $argref eq 'HASH' ? %$argref : ()
		);
	lock_hash(%arg);

	croak "Must provide sheet num or name" unless defined $arg{sheet};

	my $book = Spreadsheet::Read->new($fh, parser => $arg{type}, strip => $arg{strip});
	
	my $sheet = $book->sheet($arg{sheet}) or croak "Cannot find sheet $arg{sheet}";

	my $maxcol = $sheet->maxcol;
	my $maxrow = $sheet->maxrow;

	if (defined $arg{maxcol}) {
		if ($arg{maxcol} =~ /^\d+$/) {
			$maxcol = $arg{maxcol};
			croak "maxcol cannot exceed $maxcol" unless $arg{maxcol} <= $maxcol;
		}
		else {
			$maxcol = xls_label2col($arg{maxcol});
			croak "maxcol cannot exceed sheet max" unless $maxcol <= $sheet->maxrow;
		}
	}
	if (defined $arg{maxrow}) {
		croak "maxrow cannot exceed $maxrow" unless $arg{maxrow} <= $maxrow;
		$maxrow = $arg{maxrow};
	}

	my $linect = 1;
	my @fields;
	{
		if (defined $arg{skip}) {
			croak "Skip must be an integer" unless $arg{skip} =~ /^\d+$/;
			$linect += $arg{skip};
			croak "After skip, reached the end (line $linect)?" if $linect >= $maxrow;
		}

 		unless ($arg{header}) {
 			if (defined $arg{fields}) {
 				if (ref $arg{fields} eq 'ARRAY') {
 					croak "Number of fields is not equal to ncol"
 						unless @{$arg{fields}} == $maxcol;
 					@fields = @{$arg{fields}};
 				}
 				else {
 					croak "Custom field names must be stored in array ref";
 				}
 			}
 			else {
 				#@fields = 1..$maxcol;
 				# For spreadsheet, when it has no header line, defaults will be A-Z
 				@fields = map { $book->col2label($_) } 1..$maxcol;
 			}
 		}
 		else {
 			my @row = $arg{raw} ? $sheet->cellrow($linect) : $sheet->row($linect);
 			$linect ++;
 			@fields = @row[0..($maxcol-1)];
 			if (any { ! defined $_ || $_=~ /^\s+$/ } @fields) {
 				print STDERR join(",", @fields), "\n"; 
 				croak "Not all field names are defined";
 			}
 		}
 	}
 	my ($renamefd, $subind);
  	if (defined $arg{select}) {
  		if(any { defined $arg{$_} } qw|alias subset exclude|) {
  			warn "Options alias/subset/exclude will be ignored when select is provided"
  		}
  	 	($renamefd, $subind) = _select_subfd(\@fields, $arg{select});
  	}
  	else {
  		($renamefd, $subind) = _alias_subfd(\@fields, $arg{alias}, $arg{subset}, $arg{exclude});
  	}

 	my $iter = iterator {
  		while($linect <= $maxrow) {
  			my @row = $arg{raw} ? $sheet->cellrow($linect) : $sheet->row($linect);
  			$linect ++;
  			my @vals = @row[0..($maxcol-1)];
  			#print join(" ", @vals), "\n";
  			my $F = _hash_fdval($renamefd, \@vals, $subind, 0);
  			return $F;
  		}
		return;
	};

	if (wantarray) {
  		return($iter, [@{$renamefd}[@$subind]]);
  	}
  	else {
  		return $iter;
  	}
}


=head1 SEE ALSO

L<Iterator::Simple::Util::CSV>

=head1 AUTHOR

Xueya Zhou, C<< <xueyazhou at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-utils at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Utils>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.



=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Utils::File::Iter


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

1; # End of Utils::File::Iter

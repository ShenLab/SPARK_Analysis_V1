package Utils::File;

use strict;
use warnings;
use Carp;
use Cwd qw|abs_path|;
use IO::Detect;
use IO::Compress::Gzip;
use IO::Uncompress::Gunzip;
use PerlIO::gzip;
use File::Copy qw|move|;
use File::Temp;
use File::Basename;
use Utils::Hash qw|merge_opts|;

use base qw|Exporter|;
our @EXPORT_OK = qw|chmod_upx count_line count_lines open_file swap_files
					xls_col2label xls_label2col find_linkpath|;
our %EXPORT_TAGS = ( 'all' => \@EXPORT_OK );


=head1 NAME

Utils::File - File processing utilities.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Short-hand functions that are commonly used in file processing tasks.

    $fin = open_file($infile);
    $numlines = count_line($fin);

    chmod_upx($script);
    system($script);


=head1 EXPORT

=over 5

=item * chmod_upx

=item * count_line

=item * open_file

=item * xls_col2label

=item * xls_label2col

=back

=head1 SUBROUTINES/METHODS

=head2 chmod_upx FILE

Change file mode to be executable. Useful for HPC job submission.
See also L<File::chmod>.

=cut

sub chmod_upx {
	foreach my $filename (@_) {
		croak "File $filename does not exists" unless -f $filename;
		my $fh;
		open $fh, "<", $filename or croak "Cannot read $filename";
		my $perm = (stat $fh)[2] & 07777;
		chmod($perm | 0700, $fh);
	}
}

=head2 count_lines FILE/FH [, OPTIONS]

Count the number of lines of file. It is implemented to count the number of newlines ("\n").
It supports either plain or gzipped text file or file handle.

=head3 Options

Currently only one option: C<fastgz>, default is undefined.

When providing gzipped text file, by default we make use C<IO::Uncompress::Gunzip>,
with C<MultiStream> option on. This is to be compatible with bgzip format.
When oridnary gzip format is used, fastgz option can be turuned on to increase
performance. 

=cut

sub count_lines {
	my ($files, $argref) = @_;
	my $count;
	if (ref $files eq 'ARRAY') {
		foreach my $file (@$files) {
			$count += count_line($file, $argref);
		}
	}
	else {
		$count = count_line($files, $argref);
	}
	return $count;
}

sub count_line {
	my ($file, $argref) = @_;
	my $arg = merge_opts($argref, fastgz => undef);
	my $fh;
	if (is_filehandle $file) {
		$fh = $file;	
	}
	else {
		if ($file =~ /\.gz$/) {
			if (defined $arg->{fastgz} && $arg->{fastgz}) {
				open $fh, "<:gzip", $file or croak "Cannot open $file"; 
			}
			else {
				$fh = IO::Uncompress::Gunzip->new($file, MultiStream => 1)
					or croak "Cannot open $file";
			}
		}
		else {
			open $fh, $file or croak "Cannot open $file";
		}
	}
	my $count = 0;
	$count += tr/\n/\n/ while sysread($fh, $_, 2 ** 16);
	return $count;
}


=head2 open_file FILE [, OPTION]

Open file for reading, return a filehandle.
Supproted file types are text and gzip (based on suffix).

=head3 Options

C<fastgz>: use fast gz reading, should switch off for bgzip file.

C<encode>: use non-standard encoding, only for text file IO.

C<write> : open file for writing. The default is reading.
  NOTE: we only support basic option for writing. To append to existing
  file or pipe the input/output, users should use general function C<open>.

=cut

sub open_file {
	my ($file, $argref) = @_;
	my $arg = merge_opts($argref, fastgz => undef, encode => undef, write => undef);
	
	my $fh;
	if ($arg->{write}) {
		if ($file =~ /\.gz$/) {
			if (defined $arg->{fastgz} && $arg->{fastgz}) {
				open $fh, ">:gzip", $file or croak "Cannot open $file for write"; 
			}
			else {
				$fh = IO::Compress::Gzip->new($file) or croak "Cannot open $file for write";
			}
		}
		else {
			if (defined $arg->{encode}) {
				open $fh, ">:encoding($arg->{encode})", $file or croak "Cannot open $file for write";
			}
			else {
				#open $fh, ">", $file or croak "Cannot open $file";
				$fh = IO::File->new($file, "w") or croak "Cannot open $file for write";
			}
		}
		return $fh;
	}
	else {
		if ($file =~ /\.gz$/) {
			if (defined $arg->{fastgz} && $arg->{fastgz}) {
				open $fh, "<:gzip", $file or croak "Cannot open $file"; 
			}
			else {
				$fh = IO::Uncompress::Gunzip->new($file, MultiStream => 1)
					or croak "Cannot open $file";
			}
		}
		else {
			if (defined $arg->{encode}) {
				open $fh, "<:encoding($arg->{encode})", $file or croak "Cannot open $file";
			}
			else {
				#open $fh, $file or croak "Cannot open $file";
				$fh = IO::File->new($file) or croak "Cannot open $file";
			}
		}
		return $fh;
	}
}

=head2 swap_files

Swap two files.

=cut

sub swap_files {
	croak "swap_file: 2 arguments needed" unless @_ == 2;
	my ($file_1, $file_2) = @_;
	unless (-f $file_1 && -f $file_2) {
		croak "Not all files exist!";
	}
	my $tempfile = File::Temp::tempnam( "/tmp", "temp" );
	move($file_1, $tempfile);
	move($file_2, $file_1);
	move($tempfile, $file_2);
	return 1;
}


=head2 xls_col2label and xls_label2col

Convert column number and label.

See also: C<Spreadsheet::Read::col2label>

=cut

sub xls_label2col {
  my $label = shift;
  my $num = 0;
  croak "Not an Excel column label" unless $label =~ /^[A-Z]+$/;
  for ( my $ii = 0; $ii < length($label); $ii ++ ) {
	my $char = substr($label, $ii, 1);
	$num = $num * 26 + (ord($char) - 64);
  }
  return $num;
}

sub xls_col2label {
  my $num = shift;
  my $label = "";
  while($num > 0) {
	$num --;
	my $jj = $num % 26;
	$label = chr(65+$jj) . $label;
	$num = ($num-$jj)/26;
  }
  return $label;
}


=head2 find_linkpath

Find the real file that a symlink points to.
It can optionlly trace the link recursively to the final file.
Note: the function does not check if the final file exists.

=cut

sub find_linkpath {
	my ($file, $recurse) = @_;
	unless(-l $file) {
		warn "Input $file is not a symlink!";
		return $file;
	}
	my $indir = dirname($file); 
	my $path = readlink($file);
	if ($path =~ /^\//) {
		if (-l $path && $recurse) {
			return find_linkpath($path, $recurse);
		}
		else {
			return $path;
		}
	}
	else {
		if (-l $path && $recurse) {
			return find_linkpath($path, $recurse);
		}
		else {
			return abs_path("$indir/$path");			
		}
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

    perldoc Utils::File


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

1; # End of Utils::File

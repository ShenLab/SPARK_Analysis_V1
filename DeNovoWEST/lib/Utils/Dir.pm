package Utils::Dir;

use strict;
use warnings;
use Carp;
use File::Path qw|make_path|;
use File::Spec::Functions qw|abs2rel|;
use Utils::Hash qw|merge_opts|;

use base qw|Exporter|;
our @EXPORT_OK = qw|is_empty list_files make_tree dir_walk dir_walk_curried|;
our %EXPORT_TAGS = ( 'all' => \@EXPORT_OK );

=head1 NAME

Utils::Dir - Directory accessing utilities.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Useful functions commonly used when accessing directories.

	if ( ! -d $targetdir || is_empty($targetdir) ) {
		make_tree($sourcedir, $targetdir);
	}
	else {
		@files = list_files($targetdir, {include => /\.pl$/});
	}

=head1 EXPORT

=over 5

=item * is_empty

=item * list_files

=item * make_tree

=item * dir_walk

=item * dir_walk_curried

=back

=head1 SUBROUTINES/METHODS

=head2 is_empty DIR

Test if directory is empty.

Will croak if DIR does not exist or is not dir.

Return 1 if DIR is empty directory, 0 if DIR is not empty.

=cut

sub is_empty {
	my $indir = shift @_;
    croak "Directory $indir does not exist" if not -e $indir; 
    croak "$indir is not a directory" if not -d $indir;  
    opendir my $dir, $indir or    # likely a permissions issue
    	croak "Can't opendir $indir, because: $!\n";
    readdir $dir;
    readdir $dir;
    return 0 if( readdir $dir ); # 3rd times a charm
    return 1;
}

=head2 list_files DIR [, OPTIONS]

Recursively list all files in the directory.

See also L<File::Util::list_dir>, and L<File::Find>.

=head3 Options

=over 5

=item C<type>

Output file types, default: C<f> (file only). 
Use C<d> to output dir only; use C<f,d> or C<d,f> to output both file and directory.

=item C<include>

Include files whose name match a pattern (regex). When both include and exclude options
are give, include takes higher priority.

=item C<exclude>

Exclude files whose name does not match a pattern (regex).

=teim C<recur>

Recursively call the function to list all files in sub-directories, default: 1.

=item C<verbose>

Verbose mode, useful for debug, default: 0.

=back

=cut

sub _include_flag {
	my ($filename, $include, $exclude) = @_;
	if (defined $include) {
		return 0 unless $filename =~ /$include/;
	}
	if (defined $exclude) {
		return 0 if $filename =~ /$exclude/;
	}
	return 1;
}

sub list_files {
	my ($directory, $argref) = @_;
	my $arg = merge_opts($argref,
		type => 'f', include => undef, exclude => undef, recur => 1, verbose => 0);

	croak "Incorrect file type: $arg->{type}"
		unless $arg->{type} =~ /f/i || $arg->{type} =~ /d/i;

	my @files = ( );
	my $DIR;
	# Open the directory
 	unless(opendir($DIR, $directory)) {
  		croak "Cannot open directory $directory!\n";
	}
  	# Read the directory, ignoring special entries "." and ".."
  	# If file, print its name
  	# If directory, recursively print its contents
 	# Notice that we need to prepend the directory name!
	foreach my $file ( map { "$directory/$_" }
	  	grep (!/^\.\.?$/, readdir( $DIR )) ) {
		# If the directory entry is a regular file
		if (-f "$file") {
			if ($arg->{type} =~ /f/i) {
				print "$file\n" if $arg->{verbose};
				if ( _include_flag($file, $arg->{include}, $arg->{exclude}) ) {
					push @files, "$file";
				}
			}
		}
		# If the directory entry is a subdirectory
		elsif ( -d "$file") {
		 	# Here is the recursive call to this subroutine
		 	if ($arg->{type} =~ /d/i) {
		 		print "$file\n" if $arg->{verbose};
		 		if ( _include_flag($file, $arg->{include}, $arg->{exclude}) ) {
					push @files, "$file";
				}
		 	}
		 	if ($arg->{recur}) {
		 		push @files, list_files("$file", $argref);
		 	}
		}
	}
	return @files;
}

=head2 make_tree INDIR, OUTDIR

Create the directory tree as input directory.

See also: L<File::Path>.

=cut
 
sub make_tree {
  my ($indir, $outdir) = @_;
  my @filedirs = list_files($indir, { type => 'd' }); 
  foreach my $file_dir (@filedirs) {
	my $file_reldir = abs2rel($file_dir, $indir);
	my $out_filedir = "$outdir/$file_reldir";
	make_path($out_filedir) unless -d $out_filedir;
  }
}

=head2 dir_walk TOPDIR, FILEFUNC, DIRFUNC

Recurively go through all files in a directory and take actions.
FILEFUNC is callback method that receive the file name as input
DIRFUNC receive both the filename, and an array of results.

=head3 Examples

=over 5

=item * List file names in the directory

	sub print_filename { print $_[0], "\n" };
	dir_walk('.', \&print_filename, \&print_filename);

=item * Calculate sizes in each sub dirs

	sub file_size { -s $_[0] }
	sub dir_size {
	  my $dir = shift;
	  my $total = -s $dir;
	  for my $n (@_) { $total += $n };
	  print "%6d %s\n", $total, $dir;
	  return $total;
	}
	my $total_size = dir_walk('.', \&file_size, \&dir_size)

=item * Detecting dangling symbolic links

	sub dangles {
  	my $file = shift;
  	print "$file\n" if -l $file && ! -e $file;
	}
	dir_walk('.', \&dangles);

=item * Fetch all plain files

	@all_plain_files = dir_walk('.', sub {$_[0] } sub {shift; return @_ });

=back

=cut

sub dir_walk {
	my ($top, $filefunc, $dirfunc) = @_;
	my $DIR;

	if (-d $top) {
		my $file;
		unless (opendir $DIR, $top) {
			warn "Couldn't open directory top: $!; skipping.\n";
			return;
		}

		my @results;
		while ($file = readdir $DIR) {
			next if $file eq '.' || $file eq '..';
			push @results, dir_walk("$top/$file", $filefunc, $dirfunc);
		}
		return $dirfunc ? $dirfunc->($top, @results) : () ;
	} 
	else {
		return $filefunc ? $filefunc->($top): () ;
	}
}

=head2 dir_walk_curried

Curried version of dir_walk above.

=cut

sub dir_walk_curried {
	unshift @_, undef if @_ < 3;
	my ($top, $filefunc, $dirfunc) = @_;

	my $r;
	$r = sub {
		my $DIR;
		my $top = shift;
		if (-d $top) {
			my $file;
			unless (opendir $DIR, $top) {
				warn "Couldn't open directory $top: $!; skipping.\n";
				return;
			}

			my @results;
			while ($file = readdir $DIR) {
				next if $file eq '.' || $file eq '..';
				push @results, $r->("$top/$file");
			}
			return $dirfunc->($top, @results);
			} else {
				return $filefunc->($top);
			}
		};
	defined($top) ? $r->($top) : $r;
}


=head1 AUTHOR

Xueya Zhou, C<< <xueyazhou at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-utils at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Utils>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Utils::Dir


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

1; # End of Utils::Dir

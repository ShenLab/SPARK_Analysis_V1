package Utils::Number;

use strict;
use warnings;
use Utils::Hash;

use base qw|Exporter|;
our @EXPORT_OK = qw|commafy approx_eq pretty_bytes pretty_bps|;
our %EXPORT_TAGS = ( 'all' => \@EXPORT_OK );


=head1 NAME

Utils::Number - Simple string manipulation utilities!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Functions for simple commonly used string manipulations.

	$count = 129876491;
	print commafy($number); # 129,876,491
	$A = 9.123; B = 9.124;
	approx_eq($A, $B, 2) and say "A approx eq B";

=head1 EXPORT

=over 5

=item * commafy

=item * approx_eq

=back

=head1 SUBROUTINES/METHODS

=head2 commafy

Commafy a large number: 123456789.345 => 123,456,789.345

=cut

sub commafy {
	(my $num = shift) =~ s/\G(\d{1,3})(?=(?:\d\d\d)+(?:\.|$))/$1,/g;
  	return $num;
}

=head2 approx_eq NUM1, NUM2, DP

Compare two float point numbers to certain number of decimal places.
Making use of sprintf.

=cut

sub approx_eq {
	my ($A, $B, $dp) = @_;
    return sprintf("%.${dp}g", $A) eq sprintf("%.${dp}g", $B);
}

=head2 pretty_bytes NUMBER

Obtain the size of a file in a human readable format.

See also: L<Number::Bytes::Human>

=cut

sub pretty_bytes {
    my ( $size, $n ) =( shift, 0 );
    ++$n and $size /= 1024 until $size < 1024;
    if ($n > 0) {
    	return sprintf "%.2f%s", $size, ( qw[ B KB MB GB ] )[ $n ];
    }
    else {
    	return sprintf "%dB", $size;
    }   
}


=heead2 pretty_bps NUMBER

Obtain the length of DNA sequence in a human readable format

=cut

sub pretty_bps {
    my ( $len, $n ) =( shift, 0 );
    ++$n and $len /= 1000 until $len < 1000;
    if ($n > 0) {
        return sprintf "%.2f%s", $len, ( qw[ bp kb Mb Gb ] )[ $n ];
    }
    else {
        return sprintf "%dbp", $len;
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

    perldoc Utils::Number


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

1; # End of Utils::Number

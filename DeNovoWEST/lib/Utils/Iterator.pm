package Utils::Iterator;

use strict;
use warnings;
use Iterator;

=head1 NAME

Utils::Iterator - Make iterator behave like filehandle

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Tieable iterator should be able to be called with NEXTVAL($iterator).
I have make this module to tie standard Iterator object.

	tie *IT, 'Utils::Iterator', $iterator;
	$some_value = <IT>;
	# ... or ...
	while ($nextval = <IT>) {
	# do something with $nextval
	}
	$it = Iterator::simple $iterator;
	while($dat = $it->()) {
	# to something with $dat
	}


=head1 DESCRIPTION

I want to make the standard Iterator more like L<Iterator::Simple>, so I add 
the following code blocks at the compile time to jury-rug the standard module.

=cut
use Iterator::Simple;

BEGIN {
	package Iterator;
	sub NEXTVAL {
		my $self = shift;
		return if $self->is_exhausted();
		return $self->value();
	}
	sub simple {
		my ($self) = @_;
		return sub {
			return if $self->is_exhausted();
			return $self->value();
		};
	}
}

sub TIEHANDLE {
	my ($package, $iterator) = @_;
	my $self = { IT => $iterator };
	bless $self => $package;
}

sub READLINE {
	my $self = shift;
	$_=$self->{IT}->NEXTVAL();
	return $_;
}



=head1 AUTHOR

Xueya Zhou, C<< <xueyazhou at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-utils at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Utils>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Utils::Iterator


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

1; # End of Utils::Iterator

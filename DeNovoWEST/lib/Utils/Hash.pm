package Utils::Hash;

use strict;
use warnings;
use Carp;
use Tie::IxHash;
use IO::Prompt;
use Config::Std;
use Tie::IxHash;
use List::MoreUtils qw|natatime|;
use Hash::Util qw|lock_hash lock_hashref|;
use Scalar::Util qw(looks_like_number);

use base qw|Exporter|;
#our @EXPORT = qw|merge_opts|; 
our @EXPORT_OK = qw|merge_opts merge_conf append_conf
                    set_default chk_default peek_hash
                    str2hash array2hash array2hist 
                    hash_walk prompt_params|;
our %EXPORT_TAGS = ( 'all' => \@EXPORT_OK );

=head1 NAME

Utils::Hash - Simple hash utilities.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Implement some simple useful functions for hash.

	$args = merge_opts($argref, sort => 1, rev => 0, para => undef);

=head1 EXPORT

=over 5

=item * merge_opts

=back

=head1 SUBROUTINES/METHODS

=head2 merge_opts ARGS, OPTS

Merge user provided options with defaults.

ARGS is user provided options, either undef or a hash ref.
OPTS is default "key => value" pairs, it should include all allowed
option keys, no matter it is defined or not.

After calling this function, the return value is a B<reference> to a readonly hash.
The input $argref will also become readonly. This is aim to provide a safe way to
specify optional arguments to subroutines. It only allows options at given set
of keys to be set, and prevents accidental modification to the option values in 
the future.

=cut

sub merge_opts {
	my ($argref, %options) = @_;
	if (ref $argref eq 'HASH') {
		lock_hashref $argref;
		while(my ($k, $v) = each %$argref) {
			croak "Option $k is not allowed" unless exists $options{$k}; 
			$options{$k} = $v;
		}
	}
	# Note: lock_hash return a reference
	return lock_hash %options;
}

=head merge_conf FILE, PARAMS

Utility function: merge config and parameters.

This function can be used to combine parameters specified in the config file
and user provided key-value pairs. The function is similar to merge_opts,
but it deal with the data structure after parsing sectioned config files.

The limitation of this function is: if one slot is hash reference,
then additional parameter will override existing parameters. So we forbid
overriding in such cases.

=cut

sub merge_conf {
    my $file = shift @_;
    read_config $file => my %conf;
    my %known;
    my @params = @_;
    my $it = natatime 2, @params;
    while (my ($label, $value) = $it->()) {
        my ($sec, $key) = split(q|\.|, $label);
        croak "Must provide section.key as a label" unless defined $sec && defined $key;
        croak "$sec.$key cannot be found in the config" unless defined $conf{$sec}{$key};
        if (ref $conf{$sec}{$key} eq 'ARRAY') {
            croak "Cannot override $sec.$key, because it is a hash ref";
        }
        else {
            $conf{$sec}{$key} = $value;
        }
        $known{$label} = 1;
    }
    if (wantarray) {
        return %conf;
    }
    else {
        return \%conf;
    }
}


=head2 set_default HASH, SLOT, DEFAULT

If the HASH->{SLOT} has not been defined, set the value to DEFAULT.

=cut

sub set_default {
    my ($hash, $slot, $default) = @_;
    $hash->{$slot} = $default unless defined $hash && defined $hash->{$slot};
    return $hash;
}

=head2 chk_default HASH, SLOt, DEFAULT

The same as set_default, except that if the value HASH->{SLOT} has been defined,
then check that its value is equal to default.

=cut

sub chk_default {
    my ($hash, $slot, $default, $overwrite) = @_;
    $hash->{$slot} = $default unless defined $hash && defined $hash->{$slot};
    unless($hash->{$slot} eq $default) {
        if ($overwrite) {
            carp "$slot value is updated: $hash->{$slot} => $default";
            $hash->{$slot} = $default;
        }
        else {
            croak "$slot has value $hash->{$slot} different from $default";
        }
    }
    return $hash;
}

=head2 peek_hash HASH

Peek the first element of hash, return a key-value pair.

=cut

sub peek_hash {
    my ($hash) = @_;
    my $firstkey = (sort keys %$hash)[0];
    if (wantarray) {
        return ($firstkey, $hash->{$firstkey});
    }
    else {
        return $firstkey;
    }
}


=head2 str2hash

Convert between string represetation and key-value pairs. 

The default separator of pairs is ';' (psep), separator between key and value is '=' (kvsep).

This function is useful for parsing VCF and GFF files.

Example: ID=ENST00000456328.2;Parent=ENSG00000223972.5;gene_id=ENSG00000223972.5_2;transcript_id=ENST00000456328.2_1;gene_type=transcribed_unprocessed_pseudogene;gene_name=DDX11L1;transcript_type=processed_transcript;transcript_name=DDX11L1-202;level=2

Options: 
1. Use psep to define separators between pairs and kvsep to define separator between key and value
2. To preverse the order of keys, set order => 1
=cut

sub str2hash {
    my ($string, $argref) = @_;
    my $arg = merge_opts($argref, psep => ';', kvsep => '=', order => undef);

    my %dat;
    if ($arg->{order}) {
        tie %dat, 'Tie::IxHash';
    }
    #%dat = map { my @pair = (split('=', $_))[0,1] } split(';', $string);
    foreach my $pair (split($arg->{psep}, $string)) {
        my @pair = split($arg->{kvsep}, $pair, 2);
        $dat{$pair[0]} = $pair[1] // $pair[0];
    }

    if ($arg->{order}) {
        if (wantarray) {
            die "str2hash: Must return reference to keep order!"
        }
        return \%dat;
    }
    else {
        if (wantarray) {
            return %dat;
        }
        else {
            return \%dat;
        }
    }
}


=head2 array2hash

Create a hash from an array, that uses all array elements as keys.

=cut

sub array2hash {
    my %dat;
    @dat{@_} = @_;
    if (wantarray) {
        return %dat;
    }
    else {
        return \%dat;
    }
}

=head2 array2hist

Use a hash to track the number of appearance times for each element in a list.
The resulting hash has { element => counts } pairs, and can serve as a simple histogram.

=cut

sub array2hist {
    my %hist;
    foreach my $elem (@_) {
        $hist{$elem} ++;
    }
    if (wantarray) {
        return %hist;
    }
    else {
        return \%hist;
    }
}


=head2 hash_walk HASH, [], CALLBACK

Recursively visit complex data structure.

This is a general purpose subroutine to walk a hash structure. Such a hash walker
takes a code reference and simply calls that code for each leaf node in the hash.

See also: L<Data::Traverse> 

  my %data = (
    a => {
        ab => 1,
        ac => 2,
        ad => {
            ada => 3,
            adb => 4,
            adc => {
                adca => 5,
                adcb => 6,
            },
        },
    },
    b => 7,
    c => {
        ca => 8,
        cb => {
            cba => 9,
            cbb => 10,
        },
    },
  );

  sub print_keys_and_value {
    my ($k, $v, $key_list) = @_;
    printf "k = %-8s  v = %-4s  key_list = [%s]\n", $k, $v, "@$key_list";
  }

  hash_walk(\%data, [], \&print_keys_and_value);


=cut

sub hash_walk {
	my ($hash, $key_list, $callback) = @_;
	while (my ($k, $v) = each %$hash) {
        # Keep track of the hierarchy of keys, in case
        # our callback needs it.
        push @$key_list, $k;

        if (ref($v) eq 'HASH') {
            # Recurse.
            hash_walk($v, $key_list, $callback);
        }
        else {
            # Otherwise, invoke our callback, passing it
            # the current key and value, along with the
            # full parentage of that key.
            $callback->($k, $v, $key_list);
        }

        pop @$key_list;
    }
}

=head2 prompt_params PARLIST

Provide a list of variables, its description, and default values,
prompt users to update the parameter values, and return the results
as a hash.

  my @params = (SNVTR => "SNV Tranche" => 99.9, 
    VQX => "Variant Quality" => 100,
    GQ  => "Min. GQ for Autosome" => 60,
    MISS => "Max. Missing Rate for Autosome" => 0.02,
    HWE => "P-value Threshold for HWE-test for Autosome" => 1.0e-5,
    VQX => "Variant Quality for chrX" => 100,
    GQX => "Min. GQ for chrX" => 40,
    MISSX => "Max. Missing Rate for chrX" => 0.03,
    HWEX => "P-value Threshold for HWE for chrX" => 1.0e-4);

  my %conf = accept_prompt(@params);

=cut

sub prompt_params {
    my ($params) = @_;
    croak "The length of input params list must be a multiple of 3" if scalar(@$params)%3 > 0;
  
    sub menu_items {
        my @params = @_;
        push @params, "EXIT", "Accept and Exit", "Y";
        return join("\n", map {  " $_) $params[$_*3-2] ($params[$_*3-1])" } 1..@params/3)."\n---------\n";
    }

    my $npar = int(@$params/3)+1;
    while(my $nn = prompt(menu_items(@$params)."Select Parameters to Update: (1-$npar) ", -num, -default => $npar)) {
        next unless $nn >= 1 && $nn <= $npar;
        last if $nn == $npar;
        my $val;
        if (looks_like_number($params->[$nn*3-1])) {
            $val = prompt("Update ".$params->[$nn*3-2]." ", -num, -default => $params->[$nn*3-1]);
        }
        else {
             $val = prompt("Update ".$params->[$nn*3-2]." ", -default => $params->[$nn*3-1]);
        }
        $params->[$nn*3-1] = $val->{value};
    }
    print "Done\n";

    tie my %conf => "Tie::IxHash";
    for(my $ii = 0; $ii < $npar-1; $ii ++) {
        $conf{$params->[3*$ii]} = $params->[3*$ii+2];
    }
    if (wantarray) {
        return %conf;
    }
    else {
       return \%conf; 
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

    perldoc Utils::Hash


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

1; # End of Utils::Hash

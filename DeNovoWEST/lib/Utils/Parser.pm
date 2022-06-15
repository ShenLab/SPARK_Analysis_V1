package Utils::Parser;

use strict;
use warnings;
use Carp;
use Utils::Hash qw(merge_opts);
use Data::Dumper;
use HOP::Lexer qw(string_lexer);
use HOP::Stream qw(iterator_to_stream);
use HOP::Parser qw(:all);

use base qw|Exporter|;
our @EXPORT_OK = qw|sql_query|;


=head1 NAME

Utils::Parser - Simple SQL parser.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

	use Utils::File::Iter;
	$iter = iter_file($dbfile);
	$cb = sql_query("EFFECT !~ 'NONSYN' && (EXAC_AF == 'NA' || EXAC_AF<0.0001) ");
	while( $dat = $iter->() ) {
		if ($cb->($dat)) {
			say join("\t", @{$dat}{qw|CHR POS REF ALT|}), "\n";
		}
	}

=head1 EXPORT

Currently, it only exports C<sql_query> function.

See also L<Math::Expression>

=head1 SUBROUTINES/METHODS

=head2 sql_query QUERY, DEBUG

Dynamically create a call back function taking input an hash that implement the query
given by an SQL expression.

=cut

sub sql_query {
	my ($query, $debug) = @_;

	my @input_tokens = (
					['STRING', qr/' (?: \\. | [^'] )*  ' |
								  " (?: \\. | [^"] )*  "
								 /sx, # string must be quoted
					 sub { my $s = $_[1];
						   $s =~ s/.//; $s =~ s/.$//; # strip off quotes
						   $s =~ s/\\(.)/$1/g;
						   ['STRING', $s]
						 }                                     ],
					['AND',    qr/&{1,2}|\s+and\s/i            ],
					['OR',     qr/\|{1,2}|\s+or\s/i            ],
					['FIELD',  qr/[A-Z_][A-Z0-9_\.]*/i         ],
					['OP',     qr/[!<>=]=|[!=]~|<>|[<>=]/      ],
					['LPAREN', qr/[(]/                         ],
					['RPAREN', qr/[)]/                         ],
					['NUMBER', qr/[+-]?\d+ (?:\.\d*)? | \.\d+/x  ],
					['SPACE',  qr/\s+/, sub { "" }             ],
				   );
  my %string_version = ('>'  => 'gt', '>=' => 'ge', '==' => 'eq', "=" => 'eq',
						'<'  => 'lt', '<=' => 'le', '!=' => 'ne',
						'<>' => 'ne', '=~' => '=~', '!~' => '!~',
					   );
  my %numeric_version = ('>' => '>', '>=' => '>=', '==' => '==', "=" => '==',
						 '<' => '<', '<=' => '<=', '!=' => '!=', '<>' => '!=',
						);


  my $lexer = iterator_to_stream(string_lexer($query, @input_tokens));

  my @all_tokens;
  if ( $debug ) {
	my $lex = string_lexer($query, @input_tokens);
	while(my $token = $lex->()) {
		#print '[', $token->[0], "\t, ", $token->[1], ']', "\n";
		push @all_tokens, [$token->[0], $token->[1]];
	}
  }

  my ($cquery, $squery, $term);
  my $CQuery = parser { $cquery->(@_) };
  my $SQuery = parser { $squery->(@_) };
  my $Term   = parser { $term->(@_)   };

  $cquery = operator($Term,   [lookfor('OR'),  sub { $_[0] . ' || ' . $_[1]  } ]);
  $term   = operator($SQuery, [lookfor('AND'), sub { $_[0] . ' && ' . $_[1]  } ]);
  $squery = alternate(
					  T(concatenate(lookfor('LPAREN'),
									$CQuery,
									lookfor('RPAREN'),
								   ),
						sub { return '( '.$_[1].' )' }),
					  T(concatenate(lookfor('FIELD'),
									lookfor('OP'),
									lookfor('NUMBER')),
						sub {
						  my ($field, $op, $val) = @_;
						  croak "Do not support $op for numeric comparison"
							unless defined $numeric_version{$op};
						  my $cmp_code = '$F->{"'.($field).'"} '.
							$numeric_version{$op}.' '.$val;
						  return $cmp_code;
						}),
					  T(concatenate(lookfor('FIELD'),
									lookfor('OP'),
									lookfor('STRING')),
						sub {
						  my ($field, $op, $val) = @_;
						  croak "Do not support $op for lexcial comparison"
							  unless defined $string_version{$op};
						  my $cmp_code;
						  if ( $string_version{$op} eq '=~' ||
							   $string_version{$op} eq '!~' ) {
							$cmp_code = '$F->{"'.($field).'"} '.
							  $string_version{$op}.' m{'.$val.'}';
						  } else {
							$cmp_code = '$F->{"'.($field).'"} '.
							  $string_version{$op}.' q{'.$val.'}';
						  }
						  return $cmp_code
						}),
					 );

  my ($result) = $cquery->($lexer);

  print STDERR "CALLBACK CONDITION : $result\n" if $debug;

  my $call_back = eval ('sub { my ($F)= @_; ' . $result . '}') or die $!;

  if ($debug) {
  	return ($call_back, \@all_tokens)
  }
  else {
  	return $call_back;
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

    perldoc Utils::Parser


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

1; # End of Utils::HOP

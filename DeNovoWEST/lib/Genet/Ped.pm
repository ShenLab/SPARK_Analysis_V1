package Genet::Ped;

use strict;
use warnings;
use Carp;
use IO::Detect;
use Data::Dumper;
use Storable qw|dclone|;
use List::Util qw|sum|;
use List::MoreUtils qw|all uniq|;
use Graph::Directed;
use Utils::Hash qw|set_default merge_opts array2hash|;


=head1 NAME

Genet::Ped - Simple Pedigree Tools.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 DESCRIPTION

This module represent family structures in directed graphs to facilitate querying
simple relationships.

The information from a six-column ped file is necessary to construct pedigree
structures. Sample level information like gender and affection status are stored
as attributes of graph level. We only support minimal level information obtained
from ped file.

=head1 SYNOPSIS

    use Genet::Ped

    my $ped = Genet::Ped->new($pedfile);

    my @pos = $ped->get->po_pairs;
    my @sibs = $ped->get->sib_pairs;


=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=head2 Class->new PEDFILE [, OPTIONS]

New from a ped file. First 5 or 6 columns will be used.
PED file can be text file in plink FAM or merlin PED format, or an array ref to hash refs
represent each individual in PED file, or simply the resulting arrayref with FID as key.

Options:

* C<strict> : under strict mode, all samples should have non-missing sex, all parents
must also appear in the pedigree, and pedigree graph can have and only has one connected
component.

* C<ignore> : patterns of sample names to ignore, this is typical when technical duplicates
are included. Those samples should not appear in parents.

* C<nopheno> : no phenotype column.

* C<splitfam> : under non-strict mode, when multiple connected component exist, split them
into different families ordered by family size. Suffices .1,.2 will be added to the orignal 
FID to create new FID for splitted families. 

* C<sep> : separator between FID and suffix when C<splitfam> is enabled, default is '.'.

* C<verbose> : output debug information.

=cut

sub new {
	my ($class, $pedfile, $argref) = @_;
	my $arg = merge_opts($argref, strict => undef, ignore => undef, 
		nopheno => undef, splitfam => undef, sep => '.', verbose => 1);
	my (%samps, %known);
	if (ref $pedfile eq 'HASH') {
		%samps = %$pedfile;
	}
	elsif (ref $pedfile eq 'ARRAY') {
		foreach my $indv (@$pedfile) {
			if (defined $arg->{ignore}) {
				next if $indv->{IID} =~ /$arg->{ignore}/;
			}
			if ($known{$indv->{FID},$indv->{IID}}) {
				carp "Sample $indv->{IID} in Family $indv->{FID} already exists!" if $arg->{verbose};
				if ($arg->{strict}) {
					croak "Pedigree error";	
				}
			}
			$indv->{DAD} = undef if $indv->{DAD} eq '0';
			$indv->{MOM} = undef if $indv->{MOM} eq '0';
			if ($indv->{SEX} ne '1' && $indv->{SEX} ne '2') {
				carp "Cannot determine gender for $indv->{IID} in Family $indv->{FID}" if $arg->{verbose};
				if ($arg->{strict}) {
					croak "Pedigree error";	
				}
				else {
					$indv->{SEX} = undef;
				}
			}
			if (defined $indv->{AFF} && $indv->{AFF} ne '1' && $indv->{AFF} ne '2') {
				carp "Unknown phenotype for $indv->{IID} in Family $indv->{FID}" if $arg->{verbose};
				if ($arg->{strict}) {
					croak "Pedigree error";	
				}
				else {
					$indv->{AFF} = undef;
				}			
			}
			push @{$samps{$indv->{FID}}} => $indv;
			$known{$indv->{FID},$indv->{IID}} = 1;
		}
	}
	else {
		my $nlast = $arg->{nopheno} ? 4 : 5;
		open my $fin, $pedfile, or croak "Cannot read $pedfile";
		while(<$fin>) {
			my ($fid, $iid, $dad, $mom, $sex, $aff) = (split)[0..$nlast];
			if (defined $arg->{ignore}) {
				next if $iid =~ /$arg->{ignore}/;
			}
			if ($known{$fid,$iid}) {
				carp "Sample $iid in Family $fid already exists!" if $arg->{verbose};
				if ($arg->{strict}) {
					croak "Pedigree error";	
				}
			}
			$dad = undef if $dad eq '0';
			$mom = undef if $mom eq '0';
			if ($sex ne '1' && $sex ne '2') {
				carp "Cannot determine gender for $iid in Family $fid" if $arg->{verbose};
				if ($arg->{strict}) {
					croak "Pedigree error";	
				}
				else {
					$sex = undef;
				}
			}
			if (defined $aff && $aff ne '1' && $aff ne '2') {
				carp "Unknown phenotype for $iid in Family $fid" if $arg->{verbose};
				if ($arg->{strict}) {
					croak "Pedigree error";	
				}
				else {
					$aff = undef;
				}			
			}
			push @{$samps{$fid}} => { FID => $fid, IID => $iid, 
				DAD => $dad, MOM => $mom, SEX => $sex, AFF => $aff };
			$known{$fid,$iid} = 1;
		}
	}	

	my @allfams = keys %samps;
	my %peds;
	foreach my $fid (@allfams) {
		croak "*Family $fid already exist!" if defined $peds{$fid};
		my $ped = Graph::Directed->new();
		$ped->set_graph_attribute("FID", $fid);
		my $famsamp = $samps{$fid};
		foreach my $samp (@{$famsamp}) {
			$ped->add_vertex($samp->{IID});
			$ped->set_vertex_attribute($samp->{IID}, "SEX", $samp->{SEX});
			$ped->set_vertex_attribute($samp->{IID}, "AFF", $samp->{AFF});
		}
		my %PARSEX = (DAD => '1', MOM => '2');
		foreach my $samp (@{$famsamp}) {
			foreach my $par (qw|DAD MOM|) {
				if (defined $samp->{$par}) {
					if ($ped->has_vertex($samp->{$par})) {
						$ped->add_edge($samp->{$par} => $samp->{IID});
					}
					else {
						carp "$samp->{IID}'s $par ($samp->{$par}) in Family $fid does not exist" if $arg->{verbose};
						if ($arg->{strict}) {
							croak "Pedigree error";
						}
					}
					my $parsex = $ped->get_vertex_attribute($samp->{$par}, "SEX");
					if (defined $parsex && $parsex ne $PARSEX{$par}) {
						carp "$par of $samp->{IID} in Family $fid has incorrect SEX label" if $arg->{verbose};
						if ($arg->{strict}) {
							croak "Pedigree error";
						}
						else {
							$ped->set_vertex_attribute($samp->{$par}, "SEX", $PARSEX{$par});	
						}
					}
				}
			}
		}
		# Further validate ped
		my @cc = $ped->weakly_connected_components;
		if (@cc > 1) {
			carp "Pedigree $fid has more than one connected components" if $arg->{verbose};
			if ($arg->{strict}) {
				croak "Pedigree error";
			}
		}
		foreach my $samp (@{$famsamp}) {
			my @parents = $ped->predecessors($samp->{IID});
			croak "Sample $samp->{IID} in Family $fid has more than two parents?!"
				unless @parents <= 2;
			if (@parents == 2) {
				my @sex = map { $ped->get_vertex_attribute($_, "SEX") } @parents;
				croak "Parents of sample $samp->{IID} in Family $fid have incorrect SEX label"
					unless $sex[0] eq '1' && $sex[1] eq '2' || $sex[1] eq '1' && $sex[0] eq '2';
			}
		}
		if (@cc > 1 && $arg->{splitfam}) {
			my $subpeds = _sub_peds($ped, \@cc, $arg->{sep}, 1);
			foreach my $sfid (keys %$subpeds) {
				croak "**Family $sfid already exist!" if defined $peds{$sfid};
				$peds{$sfid} = $subpeds->{$sfid};
			}
		}
		else {
			$peds{$fid} = $ped;
		}	
	}
	return bless \%peds, $class;
}


=head2 $self->sub_ped FAMID, COMPONENTS

Split large family into multiple sub-families.

COMPNENTS is an array of IIDs used for reconstruct multiple sub-families. One individual can exist in more 
than one sub-families. COMPNENTS will be sorted based on the sample size. After splitting, sufficies 1,2,... 
will be appended to family IDs. Suffix separator is '.' by default, if it is conflict with other family IDs,
then it can be specified as by the <C<sep> option.

It returns a C<Genet::Ped> object, with all other methods avaiable.

=cut

sub _sub_peds {
	my ($ped, $comps, $sep, $verbose) = @_;
	my $famid = $ped->get_graph_attribute("FID");
	my %subpeds;
	for(my $ii = 1; $ii <= @$comps; $ii ++) {
		unless(all { $ped->has_vertex($_) } @{$comps->[$ii-1]}) {
			croak "Not all individuals in component $ii can be found";
		}
		$subpeds{"$famid$sep$ii"} = $ped->subgraph($comps->[$ii-1]);
		foreach my $iid (@{$comps->[$ii-1]}) {
			foreach my $attr (qw|SEX AFF|) {
				my $val = $ped->get_vertex_attribute($iid, $attr);
				$subpeds{"$famid$sep$ii"}->set_vertex_attribute($iid, $attr, $val);
			}
		}
	}
	if ($verbose) {
		if (@$comps > 1) {
			print STDERR "Split $famid into ", scalar(@$comps), " different families, ",
				"with sample sizes: ", join(", ", map { scalar(@$_) } @$comps), "\n";
		}
		else {
			print STDERR "Create a subset of $famid with sample size", scalar(@{$comps->[0]}), "\n";
		}
	}
	return \%subpeds;
}

sub sub_peds {
	my ($self, $famid, $comps, $argref) = @_;
	my $arg = merge_opts($argref, sep => '.', verbose => undef);
	my $ped = $self->{$famid};
	croak "Cannot find Family $famid" unless defined $ped;
	if (!defined $comps) {
		my @cc = $self->get_connected($famid);
		if (@cc == 1) {
			carp "Only one connected component found for $famid";
		}
		$comps = \@cc;
	}
	my $subpeds = _sub_peds($ped, $comps, $arg->{sep}, $arg->{verbose});
	return bless $subpeds, ref $self;
}


=head2 $self->clone FAMID(s)

Clone data structure for a subset of families.

=cut

sub clone {
	my ($self, $famids) = @_;
	unless (defined $famids) {
		return dclone $self;
	}
	else {
		my %peds;
		if (ref $famids eq 'ARRAY') {
			croak "Not all famids can be found" unless (all { defined $self->{$_} } @$famids);
			foreach my $famid (@$famids) {
				$peds{$famid} = dclone $self->{$famid};
			}
		}
		else {
			croak "Cannot find $famids" unless defined $self->{$famids};
			$peds{$famids} = dclone $self->{$famids};
		}
		return bless \%peds, ref $self;
	}
}

=head2 $self->write FILE [, OPTIONS]

Write pedigrees into a file. FILE can be file name or handle.

Parents that are not exist in original PED file will be excluded.

Depending on the new option, unconnected components can be in the same or difrerent families.

Options:

* C<famids> : output a subset of families (array ref).

* C<miss> : by default mssing gender and phenotypes will coded as 0.

* C<nophe> : do not output phenotype column.

=cut

sub write_ped {
	my ($self, $file, $argref) = @_;
	my $arg = merge_opts($argref, famids => undef, miss => '0', nopheno => undef);

	my @famids;
	if (defined $arg->{famids}) {
		croak "Must provide family IDs in array ref" unless ref $arg->{famids} eq 'ARRAY';
		@famids = @{$arg->{famids}};
		croak "Not all famids can be found" unless (all { defined $self->{$_} } @famids);
	}
	else {
		@famids = $self->get_famids;
	}

	my $fout;
	if (is_filehandle($file)) {
		$fout = $file;
	}
	else {
		$fout = IO::File->new($file, "w");
	}

	foreach my $fid (@famids) {
		foreach my $iid ($self->get_members($fid)) {
			my @row = ( $fid, $iid, $self->get_father($fid, $iid) // '0',
						$self->get_mother($fid, $iid) // '0',
						$self->get_sex($fid, $iid) // $arg->{miss} );
			unless ($arg->{nopheno}) {
				push @row, $self->get_aff($fid, $iid) // $arg->{miss};
			}

			print $fout join("\t", @row), "\n";
		}
	}
}


=head2 Accessors for Family

For methods accepting FAMID, when FAMID is not provided,  will return a hash of (FID => value) pairs.

=head3 get_famids

Return a list of family IDs.

=cut

sub get_famids {
	my $self = shift @_;
	return keys %$self;
}

=head3 get_sampsize

Get total number of samples

=cut

sub get_sampsize {
	my $self = shift;
	return sum(map { $self->get_famsize($_) } $self->get_famids);
}

=head3 get_members [FAMID]

Get all family members.

=cut

sub get_members {
	my ($self, $famid) = @_;
	if (defined $famid) {
		croak "Cannot find Family $famid" unless defined $self->{$famid};
		my @members = $self->{$famid}->vertices;
		return @members;
	}
	else {
		my %members;
		foreach my $famid ($self->get_famids) {
			$members{$famid} = [ $self->{$famid}->vertices ];
		}
		if (wantarray) {
			return %members;
		}
		else {
			return \%members;
		}	
	}
}

=head3 get_famsize [FAMID]

Get family sizes.

=cut

sub get_famsize {
	my ($self, $famid) = @_;
	if (defined $famid) {
		croak "Cannot find Family $famid" unless defined $self->{$famid};
		my $size = $self->{$famid}->vertices;
		return $size;
	}
	else {
		my %sizes;
		foreach my $famid ($self->get_famids) {
			$sizes{$famid} = $self->{$famid}->vertices;
		}
		if (wantarray) {
			return %sizes;
		}
		else {
			return \%sizes;
		}	
	}
}

=head3 get_founder [FAMID]

Get all founders in the family.

=cut

sub get_founders {
	my ($self, $famid) = @_;
	if (defined $famid) {
		croak "Cannot find Family $famid" unless defined $self->{$famid};
		my @founders = $self->{$famid}->predecessorless_vertices;
		return @founders;
	}
	else {
		my %founders;
		foreach my $famid ($self->get_famids) {
			$founders{$famid} = [ $self->{$famid}->predecessorless_vertices ];
		}
		if (wantarray) {
			return %founders;
		}
		else {
			return \%founders;
		}
	}
}

=head3 get_nonfounders [FAMID]

Get all offspring in the family.

=cut

sub get_nonfounders {
	my ($self, $famid) = @_;
	if (defined $famid) {
		croak "Cannot find Family $famid" unless defined $self->{$famid};
		my @offspring = $self->{$famid}->predecessorful_vertices;
		return @offspring;
	}
	else {
		my %offspring;
		foreach my $famid ($self->get_famids) {
			$offspring{$famid} = [ $self->{$famid}->predecessorful_vertices ];
		}
		if (wantarray) {
			return %offspring;
		}
		else {
			return \%offspring;
		}
	}
}

=head3 get_numgen 

Get the number of generations.

=cut

sub get_numgen {
	my ($self, $famid) = @_;
	if (defined $famid) {
		croak "Cannot find Family $famid" unless defined $self->{$famid};
		my $plen = $self->{$famid}->longest_path;
		my $num = defined $plen ? $plen + 1 : 1;
		return $num;
	}
	else {
		my %num;
		foreach my $famid ($self->get_famids) {
			$num{$famid} = $self->{$famid}->longest_path + 1;
		}
		if (wantarray) {
			return %num;
		}
		else {
			return \%num;
		}	
	}
}

=head3 get_connected [FAMID]

Get weakly connected components (cc).

Under strict mode, one family can have at most one cc. The presence of more than on cc
indicate some extended relationship in the family, but due to absence of parents, they are
not connected.

=cut

sub get_connected {
	my ($self, $famid) = @_;
	if (defined $famid) {
		croak "Cannot find Family $famid" unless defined $self->{$famid};
		my @cc = sort { scalar(@$b) <=> scalar(@$a) } $self->{$famid}->weakly_connected_components;
		return @cc;
	}
	else {
		my %cc;
		foreach my $famid ($self->get_famids) {
			my @cc = sort { scalar(@$b) <=> scalar(@$a) } $self->{$famid}->weakly_connected_components;
			$cc{$famid} = [ @cc ];
		}
		if (wantarray) {
			return %cc;
		}
		else {
			return \%cc;
		}	
	}
}

=head3 is_bipar [FAMID]

Test if all offspring have two parents.

=cut

sub is_bipar {
	my ($self, $famid) = @_;
	if (defined $famid) {
		croak "Cannot find Family $famid" unless defined $self->{$famid};
		return all { $self->get_parents($famid, $_) == 2 } $self->get_nonfounders;
	}
	else {
		my %bipar;
		foreach my $famid ($self->get_famids) {
			$bipar{$famid} = all { $self->get_parents($famid, $_) == 2 } $self->get_offspring;
		}
		if (wantarray) {
			return %bipar;
		}
		else {
			return \%bipar;
		}	
	}
}



=head2 Accessors for Individual

All methods accept FAMID and IID. IID is optional, when omitted, it will return a hash of 
(IID => value) pairs.

=head3 get_sex/aff FAMID, IID

Get gender/affection status for FAMID-IID

=cut

for my $attr (qw|SEX AFF|) {
	no strict 'refs';
	my $sub = "get_".lc($attr);
	my $code =<<'EOF';
sub {
	my ($self, $famid, $iid) = @_;
	croak "Cannot find Family $famid" unless defined $self->{$famid};
	my $ped = $self->{$famid};
	if (defined $iid) {
		croak "Cannot find $iid in Family $famid" unless $ped->has_vertex($iid);
		my $val = $ped->get_vertex_attribute($iid, "_VAR_");
		return $val;
	}
	else {
		my %val;
		foreach my $iid ($self->get_members($famid)) {
			$val{$iid} = $ped->get_vertex_attribute($iid, "_VAR_");
		}
		if (wantarray) {
			return %val;
		}
		else {
			return \%val;
		}
	}
}
EOF
	$code =~ s/_VAR_/$attr/g;
	*$sub = eval $code;
}


=for Comments

The previous hard coded version for get sex/affection status for FAMID-IID

sub get_sex {
	my ($self, $famid, $iid) = @_;
	croak "Cannot find Family $famid" unless defined $self->{$famid};
	my $ped = $self->{$famid};
	if (defined $iid) {
		croak "Cannot find $iid in Family $famid" unless $ped->has_vertex($iid);
		my $sex = $ped->get_vertex_attribute($iid, "SEX");
		return $sex;
	}
	else {
		my %sex;
		foreach my $iid ($self->get_members($famid)) {
			$sex{$iid} = $ped->get_vertex_attribute($iid, "SEX");
		}
		if (wantarray) {
			return %sex;
		}
		else {
			return \%sex;
		}
	}
}


sub get_aff {
	my ($self, $famid, $iid) = @_;
	croak "Cannot find Family $famid" unless defined $self->{$famid};
	my $ped = $self->{$famid};
	if (defined $iid) {
		croak "Cannot find $iid in Family $famid" unless $ped->has_vertex($iid);
		my $aff = $ped->get_vertex_attribute($iid, "AFF");
		return $aff;
	}
	else {
		my %aff;
		foreach my $iid ($self->get_members($famid)) {
			$aff{$iid} = $ped->get_vertex_attribute($iid, "AFF");
		}
		if (wantarray) {
			return %aff;
		}
		else {
			return \%aff;
		}
	}
}

=cut


=head3 get_parents FAMID, IID

Get parents for FAMID-IID.

The order of parents is sorted by sex, so when both parents are present father come first.

=cut 

sub _parents {
	my ($ped, $iid) = @_;
	my @parents = sort { $ped->get_vertex_attribute($a, "SEX") <=>
				  $ped->get_vertex_attribute($b, "SEX") } $ped->predecessors($iid);	
	return @parents;
}

sub get_parents {
	my ($self, $famid, $iid) = @_;
	croak "Must provide family ID" unless defined $famid;
	croak "Cannot find Family $famid" unless defined $self->{$famid};
	if (defined $iid) {
		croak "Cannot find $iid in Family $famid" unless $self->{$famid}->has_vertex($iid);
		return _parents($self->{$famid}, $iid);
	}
	else {
		my %par;
		foreach my $iid ($self->get_members($famid)) {
			my @parents = _parents($self->{$famid}, $iid);
			if (@parents > 0) {
				$par{$iid} = [@parents];
			}
		}
		if (wantarray) {
			return %par;
		}
		else {
			return \%par;
		}
	}
}

=head3 get_offspring FAMID, IID

Get offsprings for FAMID-IID


sub _offspring {
	my ($ped, $iid) = @_;
	my @children = $ped->successors($iid);
	return @children;
}

sub get_offspring {
	my ($self, $famid, $iid) = @_;
	croak "Must provide family ID" unless defined $famid;
	croak "Cannot find Family $famid" unless defined $self->{$famid};
	if (defined $iid) {
		croak "Cannot find $iid in Family $famid" unless $self->{$famid}->has_vertex($iid);
		return _offspring($self->{$famid}, $iid);
	}
	else {
		my %offspring;
		foreach my $iid ($self->get_members($famid)) {
			my @children = _offspring($self->{$famid}, $iid);
			if (@children > 0) {
				$offspring{$iid} = [@children];
			}
		}
		if (wantarray) {
			return %offspring;
		}
		else {
			return \%offspring;
		}
	}
}

=cut

=head3 get_father/mother FAMID, IID

Get father or mother for FAM-IID.

=cut


for my $attr (qw|father mother|) {
	no strict 'refs';
	my $sub = "get_".$attr;
	my $code =<<'EOF';
sub {
	my ($self, $famid, $iid) = @_;
	my $ped = $self->{$famid};
	croak "Cannot find Family $famid" unless defined $ped;
	if (defined $iid) {
		my @par = grep { $ped->get_vertex_attribute($_, "SEX") eq "_SEX_"} $ped->predecessors($iid);
		if (@par == 1) {
			return $par[0];
		}
		else {
			return undef;
		}
	}
	else {
		my %pars;
		foreach my $iid ($self->get_members($famid)) {
			my @par = grep { $ped->get_vertex_attribute($_, "SEX") eq "_SEX_" } $ped->predecessors($iid);
			if (@par == 1) {
				$pars{$iid} = $par[0];
			}
		}
		if (wantarray) {
			return %pars;
		}
		else {
			return \%pars;
		}
	}
}
EOF
	my $sex = $attr eq 'father' ? 1 : 2;
	$code =~ s/_SEX_/$sex/g;
	*$sub = eval $code;
}

=head3 is_founder FAMID, IID

Test if a individual is a founder.

=cut

sub is_founder {
	my ($self, $famid, $iid) = @_;
	croak "Must provide both FID and IID" unless defined $famid && defined $iid;
	if ( $self->get_parents($famid, $iid) > 0 ) {
		return 0;
	}
	else {
		return 1;
	}
}


=head2 Computed Accessors for Family

=head3 get_po_pairs FAMID

Get all parent-offspring pairs within families.

=cut

sub _po_pairs {
	my ($parents) = @_;
	my @relpairs;
	while(my ($iid, $pars) = each %$parents) {
		foreach my $parid (@$pars) {
			push @relpairs, [$parid, $iid];
		}
	}
	return @relpairs;
}


sub get_po_pairs {
	my ($self, $famid) = @_;
	if (defined $famid) {
		croak "Cannot find Family $famid" unless defined $self->{$famid};
		my @relpairs;
		my %parents = $self->get_parents($famid);	
		return _po_pairs(\%parents);
	}
	else {
		my %relpairs;
		foreach my $famid ($self->get_famids) {
			my %parents = $self->get_parents($famid);	
			push @{$relpairs{$famid}}, _po_pairs(\%parents);
		}
		if (wantarray) {
			return %relpairs;
		}
		else {
			return \%relpairs;
		}
	}
}


=head3 get_sib_pairs / get_all_sibpairs FAMID

get_sib_pairs : Get full sib pairs within families, who must share both parents 
  and parents should present in the ped file. 

get_all_sibpairs: Get all sib pairs, including those who only share one parent.


=cut

sub _pairwise {
	my @pairs;
	for (my $ii = 0; $ii < @_-1; $ii ++) {
		for (my $jj = $ii + 1; $jj < @_; $jj ++) {
			push @pairs, [$_[$ii], $_[$jj]];
		}
	}
	return @pairs;
}

sub _sib_pairs {
	my ($parents, $npar) = @_;
	$npar = 2 unless defined $npar;
	my %sibs;
	while(my ($iid, $pars) = each %$parents) {
		next unless @$pars >= $npar;
		my $parkey = join($;, @$pars);
		push @{$sibs{$parkey}}, $iid;
	}
	my @sibpairs;
	foreach my $sib (values %sibs) {
		push @sibpairs, _pairwise(@$sib);
	}
	return @sibpairs;
}

sub get_sib_pairs {
	my ($self, $famid) = @_;
	if (defined $famid) {
		croak "Cannot find Family $famid" unless defined $self->{$famid};
		my %parents = $self->get_parents($famid);
		return _sib_pairs(\%parents);
	}
	else {
		my %sibpairs;
		foreach my $famid ($self->get_famids) {
			my %parents = $self->get_parents($famid);
			$sibpairs{$famid} = [ _sib_pairs(\%parents) ];	
		}
		if (wantarray) {
			return %sibpairs;
		}
		else {
			return \%sibpairs;
		}
	}
}

sub get_all_sibpairs {
	my ($self, $famid) = @_;
	if (defined $famid) {
		croak "Cannot find Family $famid" unless defined $self->{$famid};
		my %parents = $self->get_parents($famid);
		return _sib_pairs(\%parents, 1);
	}
	else {
		my %sibpairs;
		foreach my $famid ($self->get_famids) {
			my %parents = $self->get_parents($famid);
			$sibpairs{$famid} = [ _sib_pairs(\%parents, 1) ];	
		}
		if (wantarray) {
			return %sibpairs;
		}
		else {
			return \%sibpairs;
		}
	}
}

=head2 get_trios FAMID

Get offspring-father-mother trios.

=cut

sub _ofm_trios {
	my ($parents) = @_;
	my @trios;
	while(my ($iid, $pars) = each %$parents) {
		my @pars = grep { defined $_ } @$pars;
		next unless @pars == 2;
		push @trios, [$iid, @pars];
	}
	return @trios;
}

sub get_trios {
	my ($self, $famid) = @_;
	if (defined $famid) {
		croak "Cannot find Family $famid" unless defined $self->{$famid};
		my %parents = $self->get_parents($famid);
		return _ofm_trios(\%parents);
	}
	else {
		my %trios;
		foreach my $famid ($self->get_famids) {
			my %parents = $self->get_parents($famid);
			$trios{$famid} = [ _ofm_trios(\%parents) ];
		}
		if (wantarray) {
			return %trios;
		}
		else {
			return \%trios;
		}
	}
}


=head1 AUTHOR

Xueya Zhou, C<< <xueyazhou at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-genet at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Genet>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Genet::Ped


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Genet>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Genet>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Genet>

=item * Search CPAN

L<http://search.cpan.org/dist/Genet/>

=back


=head1 ACKNOWLEDGEMENTS


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

1; # End of Genet::Ped

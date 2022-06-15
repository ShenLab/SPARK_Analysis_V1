package Utils::Workflow;

use strict;
use warnings;
use Carp;
use POSIX qw|ceil|;
use IO::Detect;
use IO::Scalar;
use IO::File;
use Cwd qw|abs_path|;
use File::Spec;
use File::Basename;
use File::Which;
use File::Path qw|make_path|;
use List::MoreUtils qw|all any none uniq natatime|;
use Hash::Util qw|lock_hash|;
use Storable qw|nstore retrieve|;
use Set::IntSpan;
use Data::Dumper;
use Graph::Directed;
#use Config::Simple;
use Config::Std;
use Utils::Hash qw|merge_opts|;
use Utils::Dir qw|is_empty|;
use Utils::List qw|run_list|;
use String::ShellQuote;

#use base qw|Exporter|;
#our @EXPORT_OK = qw|merge_conf|;
#our %EXPORT_TAGS = ( 'all' => \@EXPORT_OK );

=head1 NAME

Utils::Workflow - Workflow management.

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

	my $wkf = Utils::Workflow->new($rootdir, {engine => 'BASH'});
	$wkf->add(\*GATK_DOC, { name => "GATK_DoC" slots => 200, expect => { out => docfile } });
	$wkf->add(\*PREP_INTV, { name => "Prep_Intv", expect => { wrk => interval }});
	$wkf->add(\*PICARD_DOC, { name => "Picard_DoC", slots => 200, depend => "PREP_INTV",
		expect => { out => docfile2 } });

	use Config::Simple;
	my (%path, %resource);
	Config::Simple->import_from($pathfile, \%path)
	Config::Simple->import_from($cfgfile, \%resource);
	$wkf->inst(\%path);

	# Submit to SGE
	$wkf->run({ conf => \%resource, dryrun => $dryrun });

	# Make use of inline files to store script templates.

	__GATK_DOC__
	read BAMFILE _PARDIR_/_TASK_

	IID=$(basename $BAMFILE)
	IID=${IID%.bam}

	_PATH.JAVA_ -jar _PATH.GATK_  \
	    -T DepthOfCoverage \
	    -R _FASTA_ \
	    -o _OUT_/$IID.VCRome \
	    -I $BAMFILE \
	    -L _TARGET_ \
	    -ct 10 -ct 15 -ct 20 

	__PREP_INTV__
	for CHR in X Y; do
	    awk '$1=="X"' _PATH.TARGET_ > _WRK_/chr${CHR}.bed 
	    _PATH.JAVA_ -jar _PICARD_ BedToIntervalList \
	        I=_WRK_/chrX_tg.bed \
	        O=_WRK_/chrX_tg.$IID.interval_list \
	        SD=_WRK_/$BAMFILE 
	done

	__PICARD_DOC__
read BAMFILE _PARDIR_/_TASK_

	_PATH.JAVA_ -jar _PATH.PICARD_ \
	    CollectHsMetrics \
	    I=$BAMFILE \
	    O=_OUT_/$IID.hs_metrics.txt \
	    R=_FASTA_ \
	    BAIT_INTERVALS=_WRK_/chrX_tg.interval_list \
	    TARGET_INTERVALS=_WRK_/chrY.interval_list

=head2 DESCRIPTION

Each workflow is associated with a standard layout of working directory.
The standard set include C<src, par, out, wrk, tmp,> and C<log> subdirectories.
The content in those subdirs are

C<src> => bash scripts. Each task of the workflow will have its own script file in src,
          prefixed by the task name.
C<par> => parameter files for bash scripts. Each script is associated with one or more 
          parameter files in par with the same name prefix. 
C<out> => output files
C<wrk> => intermediate files, needs further inspection
C<tmp> => temporarily files, can remove after finish.
C<log> => log files, after execution, stdout and stderr of each task will be stored in two 
          log files in log directory.

Users can customize the subdir names to be more meaningful to the application.
Users can also specify additional subdirs as key-value pairs like final => 'vcfs'

Each task can be a single job or an job array with index starting from 1.
Each task will be described by a series of bash command with special
syntax reserved for future substitution (see below). Users can specify
the dependencies between different tasks.

=head3 Task Template

Workflow is comprised of a series of tasks that are dependent of each other. 
Usually, each task is a shell script with multiple commands. In defining workflow
we used template to represent the logic steps of each task. The purpose is to separate
those logical steps from options and associated execution enviroment, which are
specified by the workflow users and can be used to instantiate task scripts later.

Task template is just like normal bash scripts, except for the following reserved syntax.

=over 5

=item 1. Use C<_VARNAME_> variables to specify replacable variables 

Those variables are typically provided by the user via in config files. 
VARNAME should be C<qr/^[A-Z][A-Z0-9\.]+$/>. The dot is used to separate section
and variables in the config file. For example, C<__PATH.GATK__> will corrsponds
to the GATK software path given in the config file as

    [PATH]
    GATK=/XX/YY/GATK.jar

*NOTE1 * We also support variables with multiple values as an array ([A-Z][A-Z0-9\.]+)\[([^\[\]]*)\].
For example, C<_ARRAY[SEP]_> specifies an array variable. When substituted into the
script, multiple values will be separated by SEP that is given in the square bracket.

*NOTE 2* We noted that _VARNAME_ may be in conflict with long options in picard. 
To resolve this issue, we add an option (strict_var) to enforce that non-reserved 
VARNAME must  contain a separation dot '.'.  However, it may still run the risk of 
conflicting with the reserved variables like _INDEX_.


We recommend using L<Config::Std> to parse config files. When no sections is provided
for a variable, the section name should be changed to "default".

=item 2. Reserved C<VARNAME>

Path names: C<_SRCDIR_, _PARDIR_, _LOGDIR_, _OUTDIR_, _WRKDIR_, _TMPDIR_>,
Plus user provided subdirs.

And C<_TASK_, _INDEX_> will be substituted dynamically to the current task name 
and job array index.

=back


=head2 CLASS::new ROOT, [,ENGINE, OPTIONS]

Create a workflow management object. If rootdir already exist
and not empty, it will croak error. The default engine is C<'SGE'>.
Currently supported engines are: c<BASH, SGE>

Optionally users can specify additional sub directories, the corresponding special
variable names will be updated. For example, to add a directory "bwa", then C<BWADIR>
will be reserved for this special directory.

=cut

# Reserved variables.
my %RESVAR = (INDEX => 1, TASK => 1);

# Supported exec engines.
my %ENGINE = (BASH => 1,  SGE => 1);

# Inst values for INDEX
my %INDSUB = (BASH => '$1', SGE => '$SGE_TASK_ID');

# List of currently supported interpretors, will increase.
my %INTERP = (Rscript => '.R', perl => '.pl', python => '.py');

sub new {
	my ($class, $rootdir, $argref) = @_;
	#croak "Must provide engine" unless defined $args->{engine};
	my $engine;
	unless(defined $argref->{engine}) {
		my $qsub = which("qsub");
		if ($qsub) {
			$argref->{engine} = 'SGE';
		}
		else {
			$argref->{engine} = 'BASH'
		}
	}
	
	my $args = merge_opts($argref, engine => undef, dir => undef, force => undef, strict_var => undef);
	$engine = $args->{engine};

	croak "Engine $engine is not supported" unless defined $ENGINE{$engine};
	
	if (-d $rootdir && !is_empty($rootdir)) {
		unless ($args->{force}) {
			croak "The rootdir $rootdir is not empty" ;
		}
		else {
			carp "The rootdir $rootdir is not empty, data will be overwritten";
			unlink "$rootdir/.perldat" unless -f "$rootdir/.perldat";
		}
		
	}

	my %subdirs = ( srcdir => 'src', pardir => 'par', logdir => 'log', 
		tmpdir => 'tmp', wrkdir => 'wrk', outdir => 'out');

	if (defined $args->{dir}) {
		if (ref $args->{dir} eq 'ARRAY') {
			foreach my $sub (@{$args->{dir}}) {
				$subdirs{"${sub}dir"} = $sub;
			}
		}
		else {
			$subdirs{$args->{dir}."dir"} = $args->{dir};
		}
	}
	croak "Subdir names must all be lower case" unless(all { lc($_) eq $_ } keys %subdirs);
	
	lock_hash %subdirs;

	# Special path variables
	foreach my $sub (values %subdirs) {
		unless (-d "$rootdir/$sub") {
			make_path("$rootdir/$sub") or croak "Cannot create $sub sub-directory";	
		}
	}

	# Setup directory structure, initialize special parameters
	my $wkf = bless { rootdir => abs_path($rootdir), engine => $engine,
		subdirs => \%subdirs, tasks => {}, depends => Graph::Directed->new,
		strict_var => $args->{strict_var} },  $class;

	return $wkf;
}

=head2 dump

Dump obj to the current rootdir. Data will be stored as a hidden file ".perldat".

=cut

sub DESTROY {
	my $self = shift;
	my $fname = "$self->{rootdir}/.perldat";
	unlink $fname if -f $fname;
	#unless (-f $fname) {
	#local $Storable::Deparse = 1;
	# Remove callbacks before serilization
	foreach my $taskname ($self->get_all_tasks) {
		if (defined $self->{tasks}{$taskname}{callback}) {
			$self->{tasks}{$taskname}{callback} = undef;
		}
	}
	nstore $self, $fname;
	#}

	#else {
	#	my $oldself = retrieve $fname;
	#	croak "There is an old copy of different '.perldat'\nPlease check you rootdir!"
	#		unless Compare($self, $oldself);
	#}
}


=head2  $self->get_subdir [LABEL]

Return the absolute path of a subdirectory given by label.

If no label is provided, then return the list of all subdirs,
the return value will be array or hashref depending on context.

=cut 

sub get_subdir {
	my ($self, $sub, $quote) = @_;
	if (defined $sub) {
		croak "Cannot find $sub-subdir" unless defined $self->{subdirs}{"${sub}dir"};
		my $path = File::Spec->catfile($self->{rootdir}, $sub);
		if ($quote) {
			$path = shell_quote($path);
		}
		return $path;
	}
	else {
		if ($quote) {
			return map { shell_quote(File::Spec->catfile($self->{rootdir}, $_)) }
				values %{$self->{subdirs}};
		}
		else {
			return map { File::Spec->catfile($self->{rootdir}, $_) }
				values %{$self->{subdirs}};
		}
	}
}

=head2 $self->add TEMPLATE [,OPTIONS]

Add task into workflow. Task must be added in logical orders.
If B depends on A, then A must be added first. This prevent
potential loops in the dependency graph.

TEMPLATE can be a string, a file name or file handle.

=head3 Options

=over 5

=item * name

Task name. If not provided, it will introspect the template variable name.

=item * depend

List of previous added tasks that the current task depends on.
It can be scalar or arrayref.

NOTE: we use C<depend> option to manually deal with dependency between tasks, 
alternatively tasks can also be connected with using expected input/output.

=item * expect

The expected output files from this task in an arrayref.
For task with multiple slots, expected file from all slots should be listed.
Length of the array should equal to the nslots/step.
Expected files are specified in path relative to rootdir.
The expected files can be an empty list so the job will be mandatory to run.

=item * callback

A callback subroutine used to validate the expected output. If not provided, it will only
check the presence or absence of the file. The function will take expected files per job 
slot as input. The output will be 1 if all files are present and correct, and 0 otherwise.

=item * interp

Script interpretor. If unspecified, we assume the task is a bash script.
But users are free to use languages like perl, R, python, etc, and give appropriate
interpretors along with custom options. When non-bash script template is used,
a second bash wrapper script will also be written to the src directory.

=item * commargs

For non-bash script: command line arguments to the script. Those arguments
will be appended to bash script wrapper.

=item * nslots

When defined, it indicates the task will be a job array, and it gives the total
number of jobs.

=item * step

Step size that each job will iterate through slots.

=item * deparray

Indicate that each slot of current task depends on each slot of depended jobs.
Require the current job and all jobs it depends on have the same number of slots.

=back

=cut


sub add {
	my ($self, $template, $argref) = @_;

	if (ref $template eq "HASH") {
		croak "Cannot find script template" unless defined $template->{script};
		$argref = $template;
		$template = $template->{script};
		delete $argref->{script};
	}

	my $arg = merge_opts($argref,
		name => undef, expect => undef, callback => undef, depend => undef, deparray => undef,
		nslots => undef, start => 1, step => 1, interp => undef, commargs => undef);

	unless (defined $arg->{name} && defined $arg->{expect}) {
		croak "Must provide task name and expected output list";
	}

	my $taskname = $arg->{name};
	unless($taskname =~ /^[a-zA-z]\w*$/) {
		croak "Incorrect taskname: $taskname";
	}
	if (defined $self->{tasks}{$taskname}) {
		croak "Task $taskname already exists!";
	}

	my $fh;
	if (is_filehandle($template)) {
		$fh = $template;
	}
	else {
		if ($template !~ /\n/ && -f $template) {
			$fh = IO::File->new($template);
		}
		else {
			$fh = IO::Scalar->new(\$template);
		}
	}
	$self->{tasks}{$taskname}{template} = join("", <$fh>);

	# Order of the task is determined by the time it is added to the workflow
	$self->{tasks}{$taskname}{order} = scalar(keys %{$self->{tasks}});
	

	# Record all expected outputs
	if (defined $arg->{nslots}) {
		croak "Number of slots must an integer" 
			unless $arg->{nslots} =~ /^\d+$/ && $arg->{nslots} > 0;
		$arg->{start} = 1 unless defined $arg->{start};
 		$arg->{step} = 1 unless defined $arg->{step};
		croak "Start index must be an integer"
			unless $arg->{start} =~ /^\d+$/ && $arg->{start} > 0;
		croak "Step size must be an integer"
			unless $arg->{step} =~ /^\d+$/ && $arg->{step} > 0;
		croak "Number of slots must be a multiple of step size" 
			if $arg->{nslots} % $arg->{step} > 0;
		unless (@{$arg->{expect}} == $arg->{nslots}/$arg->{step}) {
			printf STDERR "$arg->{name}: n_exp=%d, n_slots=%d, step=%d\n",
				scalar(@{$arg->{expect}}), $arg->{nslots}, $arg->{step};
			croak "Expected output are not listed for all jobs in a jobarray";
		}
		
			
		$self->{tasks}{$taskname}{nslots} = $arg->{nslots};
		$self->{tasks}{$taskname}{start} = $arg->{start};
		$self->{tasks}{$taskname}{step} = $arg->{step};	
	}

	$self->{tasks}{$taskname}{expect} = $arg->{expect};
	if (defined $arg->{callback}) {
		$self->{tasks}{$taskname}{callback} = $arg->{callback};
	}

	# Adding edges to depencency graph
	if ($arg->{depend}) {
		if (ref $arg->{depend} eq 'ARRAY') {
			croak "Not all dependent tasks have been defined"
				unless(all { defined $self->{tasks}{$_} } @{$arg->{depend}}); 
			if ($arg->{deparray}) {
				unless(defined $arg->{nslots} &&
					(all { defined $self->{tasks}{$_}{nslots} &&
						$self->{tasks}{$_}{nslots} == $arg->{nslots} &&
						$self->{tasks}{$_}{start} == $arg->{start} } @{$arg->{depend}})) {
					foreach my $jobname ( @{$arg->{depend}}) {
						print STDERR join("\t", $jobname, $self->{tasks}{$jobname}{start}, $self->{tasks}{$jobname}{nslots}), "\n";
					}
					croak "With deparray, all tasks should have the same number of slots as the depended";
				}
				$self->{tasks}{$taskname}{deparray} = 1;
			}
			$self->{tasks}{$taskname}{depend} = $arg->{depend};
		}
		else {
			croak "Dependent task $arg->{depend} has not been defined"
				unless defined $self->{tasks}{$arg->{depend}};
			if ($arg->{deparray}) {
				unless(defined $arg->{nslots} && defined $self->{tasks}{$arg->{depend}}{nslots} &&
					$self->{tasks}{$arg->{depend}}{nslots} == $arg->{nslots} &&
					$self->{tasks}{$arg->{depend}}{start} == $arg->{start} ) {
					my $jobname = $arg->{depend};
					print STDERR join("\t", $jobname, $self->{tasks}{$jobname}{start}, $self->{tasks}{$jobname}{nslots}), "\n";
					croak "With deparray, the task should have the same number of slots and start as the depended";
				}
				$self->{tasks}{$taskname}{deparray} = 1;
			}
			$self->{tasks}{$taskname}{depend} = [ $arg->{depend} ];
		}
		foreach my $deptask ( @{$self->{tasks}{$taskname}{depend}} ) {
			$self->{depends}->add_edge($taskname, $deptask);
		}
	}
	else {
		$self->{depends}->add_vertex($taskname);
	}

	# Then find all replacable variables for this task
	my @allvars = ($self->{tasks}{$taskname}{template} =~ /(?<!_)_([A-Z][A-Z0-9\.]+)_(?!_)/g);
	if ($self->{strict_var}) {
		my @strictvars;
		foreach my $var (@allvars) {
			if ($var =~ /^[A-Z][A-Z0-9]+\.[A-Z][A-Z0-9]+$/) {
				push @strictvars, $var;
			}
			else {
				if (defined $RESVAR{$var} || exists $self->{subdirs}{lc($var)}) {
					push @strictvars, $var;
				}
				#else {
					#if ($arg->{verbose}) {
					#	warn "Under strict_var mode, $var is ignored\n";
					#}
				#}
			}
		}
		$self->{tasks}{$taskname}{vars} = [@strictvars];
	}
	else {
		$self->{tasks}{$taskname}{vars} = [@allvars];
	}
 
	my @allarrays = ($self->{tasks}{$taskname}{template} =~ /(?<!_)_([A-Z][A-Z0-9\.]+)\[[^\[\]]*\]_(?!_)/g);
	# no conflict with reserved or special variables
	foreach my $avar (@allarrays) {
		if (defined $RESVAR{$avar} || exists $self->{subdirs}{lc($avar)}) {
			croak "Array variable $avar is conflict with reserved or special variables";
		}
	}
	$self->{tasks}{$taskname}{arrays} = [@allarrays];

	# For non-bash scripts
	if (defined $arg->{interp}) {
		my $itp = basename( (split(/\s+/, $arg->{interp}))[0] );
		croak "Unsupported script interpretor: $itp" unless defined $INTERP{$itp};
		$self->{tasks}{$taskname}{interp} = $arg->{interp};
		$self->{tasks}{$taskname}{suffix} = $INTERP{$itp};
		if (defined $arg->{commargs}) {
			if (ref $arg->{commargs}) {
				croak "Command line arguments must be a string variable";
			}
			$self->{tasks}{$taskname}{commargs} = $arg->{commargs};	
		}
	}
	

	return $self;
}

=head2 $self->get_all_tasks

Return the ordered list of current tasks names.

=cut

sub get_all_tasks {
	my ($self) = @_;
	return sort { $self->{tasks}{$a}{order} <=> $self->{tasks}{$b}{order} }
		keys %{$self->{tasks}}
}


=head2 $self->check_expected TASK [,SLOT]

Check expected output files for given task and slots.
If a callback is specified for expected output, it will be called to check the output file in
addition to its presence.

Note: slot index start from 1 regardless of specified start index. When step>1, then multiple slots
will be associated with the same substask.

It will return the list of expected files that are not correct or missing.

=cut

sub get_expected {
	my ($self, $taskname, $ii) = @_;
	my $nslots = $self->get_num_slots($taskname);
	my @expected;
	if (defined $ii) {
		croak "Slot $ii cannot be found for $taskname"
			unless defined $nslots && $ii <= $nslots;
		$self->{tasks}{$taskname}{step} = 1 unless $self->{tasks}{$taskname}{step};
		my $jj = ceil($ii/$self->{tasks}{$taskname}{step});
		my $exp = $self->{tasks}{$taskname}{expect}[$jj-1];
		if (ref $exp eq 'ARRAY') {
			@expected = @$exp;
		}
		else {
			@expected = ($exp);
		}
	}
	else {
		if (defined $nslots) {
			foreach my $exp (@{$self->{tasks}{$taskname}{expect}}) {
				if (ref $exp eq 'ARRAY') {
					push @expected, @$exp;
				}
				else {
					push @expected, $exp;
				}
			}
		}
		else {
			my $exp = $self->{tasks}{$taskname}{expect};
			if (ref $exp eq 'ARRAY') {
				@expected = @$exp;
			}
			else {
				@expected = ($exp);
			}
		}
	}
	my @expfiles = map { "$self->{rootdir}/$_" } @expected;
	return @expfiles;
}

sub check_expected {
	my ($self, $taskname, $ii) = @_;
	my $nslots = $self->get_num_slots($taskname);
	my $cb = $self->{tasks}{$taskname}{callback};
	unless (defined $cb) {
		$cb = sub { @_>0 && (all { -f "$self->{rootdir}/$_" } @_) }
	}
	
	if (defined $ii) {
		croak "Slot $ii cannot be found for $taskname"
			unless defined $nslots && $ii <= $nslots;
		$self->{tasks}{$taskname}{step} = 1 unless $self->{tasks}{$taskname}{step};
		my $jj = ceil($ii/$self->{tasks}{$taskname}{step});
		my $exp = $self->{tasks}{$taskname}{expect}[$jj-1];
		if (ref $exp eq 'ARRAY') {
			return $cb->(@$exp);
		}
		else {
			return $cb->($exp);
		}
	}
	else {
		if (defined $nslots) {
			my @checked;
			foreach my $exp (@{$self->{tasks}{$taskname}{expect}}) {
				if (ref $exp eq 'ARRAY') {
					push @checked, $cb->(@$exp);
				}
				else {
					push @checked, $cb->($exp);
				}
			}
			#croak "Must specify the slot for getting expected output from job array!";
			if (all { $_ == 1 } @checked) {
				return 1;
			}
			else {
				return 0;
			}
		}
		else {
			my $exp = $self->{tasks}{$taskname}{expect};
			if (ref $exp eq 'ARRAY') {
				return $cb->(@$exp);
			}
			else {
				return $cb->($exp);
			}
		}
	}
}


=head2 $self->get_num_slots TASK

Return the number of slots for a task. 

=cut

sub get_num_slots {
	my ($self, $taskname) = @_;
	return $self->{tasks}{$taskname}{nslots};
}

=head2 $self->get_slots_index TASK

Return the index of each subtask. 

=cut

#sub get_slots_index {
#	my ($self, $taskname) = @_;
#	my $st = $self->{tasks}{$taskname}{nslots};
#	my $len = $self->{tasks}{$taskname}{nslots};
#	return $st..($st+$len-1);
#}


=head2 get_dep_array

Test if current job is an job array and depend on other arrays.

=cut

sub get_dep_array {
	my ($self, $taskname) = @_;
	if (defined $self->{tasks}{$taskname}{nslots} && 
		defined $self->{tasks}{$taskname}{deparray}) {
		return 1;
	}
	else {
		return 0;
	}
}


=head2 $self->inst [CONF, REF]

When no parameters provided, it generates a list of all user configurable variables. 

If CONF is provided, the function validate and replace variables with provided values, then write
scripts into src directory. CONF can be a config file name and will be read by L<Config::Simple>.
It can also be a data strcutre returned by L<Config::Simple> or L<Config::Std>.

If REF is provided, a sanity check on CONF will be performed. If variables not in REF is defined
in CONF, it is most likely an error in CONF.

Note 1: C<Config::Simple> and C<Config::Std> have different syntax for specifying multiple values.

Note 2: Now support substitution of array variable, and check the variable type.

=cut

sub inst {
	my ($self, $conf, $ref) = @_;
	
	my $srcdir = $self->get_subdir("src");
	#croak "Shell scripts have already been generated" unless is_empty($srcdir);
	
	my %config;
	unless (defined $conf) {
		return $self->_all_envars();
	}
	else {
		if (ref $conf eq 'HASH' || ref $conf eq 'Config::Std::Hash') {
			my @sections = keys %$conf;
			# "Plattern" the config file data
			if (all { ref $conf->{$_} eq 'HASH' } @sections) {
				foreach my $sec (@sections) {
					foreach my $var (keys %{$conf->{$sec}}) {
						if ($sec eq '') {
							$config{"default.$var"} = $conf->{$sec}{$var};
						}
						else {
							$config{"$sec.$var"} = $conf->{$sec}{$var};
						}
					}
				}
			}
			elsif (none { ref $config{$_} } @sections) {
				%config = %$conf;
			}
			else {
				croak "Incorrect data structure of conf";
			}
		}
		else {
			croak "Cannot find config file" unless -f $conf;
			Config::Simple->import_from($conf, \%config);
		}
	}

	my @allvars = $self->_all_envars();

	# check if all variables are defined
	my @undefvars = grep { !defined $config{$_} } @allvars;
	croak "The following env variables are not defined: ".join(" ", @undefvars)
		unless @undefvars == 0;

	# reverse check if special reserved vars appear in config
	my @sprsvars = grep { exists $self->{subdirs}{lc($_)} ||
						  defined $RESVAR{$_} } keys %$conf;
	croak "The following variables should be reserved: ".join(" ", @sprsvars)
		unless @sprsvars == 0;

	# now substitute variables
	foreach my $taskname ($self->get_all_tasks) {
		# All user provided env variables
		foreach my $envar ($self->_env_vars($taskname)) {
			my $envar_cf = $envar =~ /\./ ? $envar : 'default.'.$envar;
			if (ref $config{$envar_cf}) {
				croak "Variable $envar_cf must be of scalar type!";
			}
			$self->{tasks}{$taskname}{template} =~ s/(?<!_)_${envar}_(?!_)/$config{$envar_cf}/g;
		}
		# All special variables
		foreach my $spvar ($self->_sp_vars($taskname)) {
			(my $label = $spvar) =~ s/DIR$//; 
			my $path = $self->get_subdir(lc($label), 1);
			$self->{tasks}{$taskname}{template} =~ s/(?<!_)_${spvar}_(?!_)/$path/g;
		}
		# All reserved variables
		foreach my $rsvar ($self->_res_vars($taskname)) {
			if ($rsvar eq 'TASK') {
				$self->{tasks}{$taskname}{template} =~ s/(?<!_)_TASK_(?!_)/$taskname/g;
			}
			elsif ($rsvar eq 'INDEX') {
				my $engine = $self->{engine};
				$self->{tasks}{$taskname}{template} =~ s/(?<!_)_INDEX_(?!_)/$INDSUB{$engine}/g;
			}
		}
		# Now substitute all user defined array variables
		foreach my $avar ($self->_array_vars($taskname)) {
			my $avar_cf = $avar =~ /\./ ? $avar : 'default.'.$avar;
			# check that config is an array
			if (ref $config{$avar} eq 'ARRAY') {
				# for each occurence of array variable, 
				# first identify the separator (may be different at different occurence)
				while ( $self->{tasks}{$taskname}{template} =~ /(?<!_)_${avar}\[([^\[\]]*)\]_(?!_)/g ) {
					my $sep = $1;
					my $replstr = join($sep, @{$config{$avar_cf}});
					$self->{tasks}{$taskname}{template} =~ s/(?<!_)_${avar}\[([^\[\]]*)\]_(?!_)/$replstr/;
				}
			}
			else {
				#carp "Array variable $avar is scalar in config";
				$self->{tasks}{$taskname}{template} =~ s/(?<!_)_${avar}\[([^\[\]]*)\]_(?!_)/$config{$avar_cf}/g;
			}
		}

		# now writing output
		if (defined $self->{tasks}{$taskname}{suffix}) {
			my $fout = IO::File->new("$srcdir/${taskname}".
					$self->{tasks}{$taskname}{suffix}, "w");
			print $fout $self->{tasks}{$taskname}{template};
			my $fs = IO::File->new("$srcdir/$taskname.sh", "w");
			print $fs '#!/bin/bash',"\n";
			print $fs "\nset -o pipefail\nset -e\n\n";
			print $fs $self->{tasks}{$taskname}{interp}, " ",
				"$srcdir/${taskname}".$self->{tasks}{$taskname}{suffix};
			if ($self->{tasks}{$taskname}{commargs}) {
				print $fs " ", $self->{tasks}{$taskname}{commargs}, "\n";
			}
			else {
				print $fs "\n";
			}
			$fs->close;
		}
		else {
			my $fout = IO::File->new("$srcdir/$taskname.sh", "w");
			print $fout '#!/bin/bash',"\n";
			print $fout "\nset -o pipefail\nset -e\n\n";
			print $fout $self->{tasks}{$taskname}{template};
			$fout->close;
		}
		chmod(0770, "$srcdir/$taskname.sh")
	}
	return $self;
}


# helper functions for getting variables from tasks
# also include arrays
sub _all_envars {
	my ($self) = @_;
	my @allvars = uniq sort map { $self->_env_vars($_) } $self->get_all_tasks;
	push @allvars, uniq sort map { $self->_array_vars($_) } $self->get_all_tasks;

	my @env_vars;
	foreach my $var (@allvars) {
		my @subs = split(q|\.|, $var);
		if (@subs == 1) {
			push @env_vars, [q|default|, $subs[0]];
		}
		elsif (@subs == 2) {
			push @env_vars, [@subs];
		}
		else {
			croak "More than two levels of hierachy for vars is not supported: $var!"
		}	
	}
	my @envars_jn = uniq sort map { $_->[0].".".$_->[1] } @env_vars;
	croak "Inconsistent number of env_vars" unless @allvars == @envars_jn; 

	if (wantarray) {
		return @envars_jn;
	}
	else {
		my %envars;
		foreach my $var (@env_vars) {
			$envars{$var->[0]}{$var->[1]} = 1;
		}
		return \%envars;
	}
}


# all special path variables for a task
sub _sp_vars {
	my ($self, $taskname) = @_;
	return grep { exists $self->{subdirs}{lc($_)} } @{$self->{tasks}{$taskname}{vars}};
}

# all reserved variables for a task
sub _res_vars {
	my ($self, $taskname) = @_;
	return grep { defined $RESVAR{$_} }	@{$self->{tasks}{$taskname}{vars}} ;
}

# all remaining enviroment variables for a task
sub _env_vars {
	my ($self, $taskname) = @_;
	return grep { !defined $RESVAR{$_} && !exists $self->{subdirs}{lc($_)} }
		@{$self->{tasks}{$taskname}{vars}};	
}

# all array variables
sub _array_vars {
	my ($self, $taskname) = @_;
	return @{$self->{tasks}{$taskname}{arrays}};
}

=head2 $self->run OPTIONS

The workflow manager will make use the expected outputs to determine un-executed jobs.
If expected output exist, then we assume the task has been finished successfully and 
will be skipped. This default behavior will be overriden by the task dependencies. 
For example, based on the missing output file, task A is rerun. Then all subsequent 
jobs that depends on A will also be rerun regardless of the presence or absence of
the expected output file. When execution under BASH engine, outputs will be checked
after each step. 

NOTE: the existance of expected output does not guarantee the correctness of the
executions of the entire workflow. The content of log files and output files should be 
examined for further details.

=head3 Options

=over 5

=item * conf

Use CONF to provide engine specific information. It can be a data structure
or a ini format configure file, the same as CONF accepted by
C<$self->inst [CONF]> method.

- For SGE, the config should provide information about resource allocation 
(typically given in C<qsub -l> option). The default resource allocations should 
be given under the key 'default', i.e. in the default section of config file.
Then task-specific customization should be specified in sections headed by
the tasknames. If the key starts with '-', it will be interpreted as extra
command line arguments to SGE. By default, we already have "-S /bin/bash -cwd"
and other standard options providing job name, script, resources, and IO redirection. 

 - For BASH, currently it only execute all tasks sequentially. For job array
it uses parallel feature to speedup. User can specify options for C<parallel>
for examples number of jobs executed at the same time).

The runtime config parameters should be written to PARDIR directory.

=item * tasks

Specify the tasks to be executed. When this option not defined, all jobs will be
executed/submitted according to the default rules above. To execute only a subset of 
jobs, we can provide task names as arrayref, hashref or a scalar (comma or semicolon
separated).

The default execution behavior still applies. I.e., only tasks with incomplete 
expected files will be executed. And all downstream tasks that depend on the current 
task will also be mandatory to run even if they are not specified by this option,
unless the C<nochain> option is switched on.

=item * nochain

Only execute specified tasks, but do not trigger downstream mandatory tasks to run.

=item * dryrun

Only print out commands, but do not execute.

=item * interact

BASH only. Under interactive mode, the workflow will not die if expected files are missing.

=back

=cut

sub run {
	my ($self, $argref) = @_;
	my $arg = merge_opts($argref, conf => undef, qsub => "",
		 tasks => undef, nochain => undef, dryrun => undef, interact => undef);

	# Check scripts
	my $srcdir = $self->get_subdir("src");
	foreach my $taskname ($self->get_all_tasks) {
		croak "No script is found for task $taskname" unless -f "$srcdir/$taskname.sh"; 
	}
	$srcdir = shell_quote($srcdir);

	# Read or process config parameters
	my $conf = $arg->{conf};
	my %config;
	if (!defined $conf) {
		%config = ();
	}
	elsif (ref $conf eq 'HASH' || ref $conf eq 'Config::Std::Hash') {
		my @sections = keys %$conf;
	 	if (all { ref $conf->{$_} eq 'HASH' } @sections) {
	 		%config = %{$conf};
	 	}
	 	else {
	 		foreach my $sec (@sections) {
	 			my @secs = split(q|\.|, $sec);
	 			if (@secs > 2) {
	 				croak "Incorrect key for resource config: ".join(" ", @secs), "\n";
	 			}
	 			elsif (@secs == 2) {
	 				if (defined $config{$secs[0]}{$secs[1]}) {
	 					carp "Parsing config: resource for $secs[0].$secs[1] is already defined, skip.";
	 					next;
	 				}
	 				unless(defined $self->{tasks}{$secs[0]}) {
	 					carp "Parsing config: task $secs[0] cannot be found in the workflow, ignore.";
	 					next;
	 				}
	 				$config{$secs[0]}{$secs[1]} = $conf->{$sec};
	 			}
	 			else {
	 				if (defined $config{default}{$secs[0]}) {
						carp "Resource for default.$secs[0] is already defined, skip.";
	 					next;
	 				}
	 				$config{default}{$secs[0]} = $conf->{$sec};
	 			}
	 		}
	 	}
	}
	else {
		croak "Cannot find config file" unless -f $conf;
		read_config $conf => %config;
		if (defined $config{""}) {
			if (defined $config{default}) {
				carp "Default section already found";
				next;
			}
		}
		else {
			$config{default} = $config{""};
		}
	}

	# Store run-time configs if not under dryrun mode
	if ($arg->{dryrun}) {
		print STDERR "Under the dry-run mode.\n";
	}
	my $fname = "$self->{rootdir}/.confdat";
	unlink $fname if -f $fname;
	nstore \%config, $fname;
		
	# Determine tasks to run, we will start with all tasks,
	# then focus on those specified if tasks are explicitly given
	my $runmode = $arg->{tasks}; # <--- name may be confusing
	my @all_tasks = $self->get_all_tasks;
	if (defined $runmode) {
		if (ref $runmode eq 'HASH') {
			if (keys %$runmode) {
				foreach my $taskname (keys %$runmode) {
					croak "Unknown task: $taskname" unless defined $self->{tasks}{$taskname};
				}
				@all_tasks = grep { defined $runmode->{$_} } $self->get_all_tasks;
			}
			else {
				croak "No jobs found!";
			}
		}
		elsif (ref $runmode eq 'ARRAY') {
			if (@$runmode) {
				foreach my $taskname (@$runmode) {
					croak "Unknown task: $taskname" unless defined $self->{tasks}{$taskname};
				}
				@all_tasks = sort { $self->{tasks}{$a}{order} <=> $self->{tasks}{$b}{order} } @$runmode;
			}
			else {
				croak "No job is found!";
			}
		}
		else {
			@all_tasks = split(/[,;]/, $runmode);
			if (@all_tasks) {
				foreach my $taskname (@all_tasks) {
					croak "Unknown task: $taskname" unless defined $self->{tasks}{$taskname};
				}
				@all_tasks = sort { $self->{tasks}{$a}{order} <=> $self->{tasks}{$b}{order} } @all_tasks;
			}
			else {
				croak "No jobs found!";
			}	
		}
	}

	# Determine which jobs (or job slots) among all_tasks will be rerun
	# based on the expected output 
	# Mandatory rerun jobs will be added later
	my %tasks2run;
	foreach my $taskname (@all_tasks) {
		my $num_slots = $self->get_num_slots($taskname);
		# Flag the presence of unfinished task(s)
		my $flag;
		if ($num_slots) {
			# Store the index to unfinished slots of the task 
			my @unexec;
			foreach my $ii (1..$num_slots) {
				unless ($self->check_expected($taskname, $ii)) {
					push @unexec, $ii;
				}
			}
			if (@unexec > 0) {
				$tasks2run{$taskname} = Set::IntSpan->new([@unexec]);
			}
		}
		else {
			unless ($self->check_expected($taskname) ) {
				$tasks2run{$taskname} = 1;
			}
		}
	}

	unless ($arg->{nochain}) {
		my %alljobs;
		# Store job IDs for all dependent jobs given those tasks in %tasks2run
		foreach my $job (keys %tasks2run) {
			$alljobs{$job} = 1;
		}
		my @mandatory = $self->{depends}->all_predecessors(keys %tasks2run);
		foreach my $job (@mandatory) {
			$alljobs{$job} = 1;
		}
		foreach my $taskname (grep { !defined $tasks2run{$_} } @mandatory) {
			print STDERR "Task $taskname is mandatory to run.\n";
		}
		$self->_recur_lookup(\%tasks2run, \%alljobs, @mandatory);
	}

	my @tasks2run_names = sort { $self->{tasks}{$a}{order} <=> $self->{tasks}{$b}{order} } keys %tasks2run;
	if (@tasks2run_names) {
		print STDERR "Jobs to be executed: ", join(", ", @tasks2run_names), "\n";
	}
	else {
		print STDERR "All expected outputs are found, no job will be executed\n";
	}
	
	my %jobids; # for storing job IDs for SGE
	my $counter = 1;
	my $logdir = $self->get_subdir("log", 1);
	foreach my $taskname (@tasks2run_names) {
		my $task_cmd;
		# Set::IntSpan or scalar
		my $slots = $tasks2run{$taskname};

		# Now run or submit jobs
		if ($self->{engine} eq 'BASH') {		
			my $parallel = which('parallel');
			print STDERR '=====>', $taskname, "\n";
			my $task_cmd;

			if (ref $slots eq 'Set::IntSpan') {
				#croak "Cannot find parallel on the system" unless defined $parallel;
				my $step = $self->{tasks}{$taskname}{step};
				my $offset = $self->{tasks}{$taskname}{start}-1;
				my @slots = map { $_ + $offset } $slots->elements;
				croak "The number of slots that need running is not a multiple of step"
					if scalar(@slots) % $step > 0;

				my $tmpdir = $self->get_subdir("tmp");
				my $fout = IO::File->new("$tmpdir/$taskname.unexec.lst", "w") 
					or croak "Cannot write $tmpdir/$taskname.unexec.lst";
				print $fout join("\n", map { $slots[$_] } grep { $_ % $step == 0 } 0..$#slots), "\n";
				$fout->close;			
				$tmpdir = shell_quote($tmpdir);
				unless (defined $parallel) {
					carp "Cannot find parallel on the system, fall back to while loop";
					$task_cmd = <<EOF
cat $tmpdir/$taskname.unexec.lst | while read II; do
	$srcdir/$taskname.sh \$II >$logdir/$taskname.\$II.out 2>$logdir/$taskname.\$II.err
done
EOF
				}
				else {
					my $conf = $config{$taskname} // $config{default};
					my $opt = join(" ", map { length($_) > 1 ? "--$_ $conf->{$_}" : "-$_ $conf->{$_}" }	keys %$conf);
					$task_cmd = qq{cat $tmpdir/$taskname.unexec.lst | parallel --eta $opt "$srcdir/$taskname.sh {} 1>$logdir/$taskname.{}.out 2>$logdir/$taskname.{}.err"\n};
				}
			}
			else {
				# Safeguard
				if ($self->get_num_slots($taskname)) {
					croak "$taskname should not be a job array";
				}
				$task_cmd = "$srcdir/$taskname.sh > $logdir/$taskname.out 2> $logdir/$taskname.err\n";	
			}
			print $task_cmd;
			unless ($arg->{dryrun}) {
				system($task_cmd);
				unless ($self->check_expected($taskname)) {
					my @expected = $self->get_expected($taskname);
					my @notfound = grep { ! -f $_ } @expected;
 					unless ($arg->{interact}) {
						if (@expected) {
							croak "Not all expected outputs for $taskname can be found";
							print STDERR join("\n", @notfound), "\n";
						}
					}
					else {
						if (@expected) {
							carp "Not all expected output for $taskname can be found";
							print STDERR join("\n", @notfound), "\n";
						}
					}					
				}
			}
		}
		elsif ($self->{engine} eq 'SGE') {
			croak "Cannot find qsub!" unless which('qsub');

			# submit jobs and capture job IDs
			print STDERR "=====> ", $taskname, "\n";

			my $conf = $config{$taskname} // $config{default};	
			my $qsubopt = "";
			if ($conf) {
				my @resources = grep { ! /^\-/ } keys %$conf;
				my @options = grep { /^\-/ } keys %$conf;
			 	$qsubopt =  join(" ", map { $conf->{$_} ? "$_ $conf->{$_}" : $_ } @options) .
			 		' -l ' . join(',', map { $conf->{$_} ? $_.'='.$conf->{$_} : $_ } @resources)
			}
			my $qsub = "qsub -S /bin/bash -cwd -N $taskname $qsubopt";

			# Array dependency with jobarray will be specified later
			# to account for the situiation when array is split
			if (defined $self->{tasks}{$taskname}{depend} &&
				(any { defined $tasks2run{$_} } @{$self->{tasks}{$taskname}{depend}})){
				unless (defined $self->{tasks}{$taskname}{deparray}) {
					my @depjobs = grep { defined $jobids{$_} } @{$self->{tasks}{$taskname}{depend}};
					if (all { defined $jobids{$_} } @depjobs) {
						my @depids = map { @{$jobids{$_}} } @depjobs;
						$qsub .= " -hold_jid " . join(',', @depids);
					}
					else {
						croak "Cannot find IDs for jobs that $taskname depends on";
					}
				}
			}


			my $task_cmd;
			if (ref $slots eq 'Set::IntSpan') {
				$qsub .= " -o $logdir/$taskname.\\\$TASK_ID.out -e $logdir/$taskname.\\\$TASK_ID.err";
				
				my $step = $self->{tasks}{$taskname}{step};
				my $offset = $self->{tasks}{$taskname}{start}-1;
				my @ranges = run_list (map { $_ + $offset } $slots->elements);	

				foreach my $rng (@ranges) {
					my $depopt = "";
					if (defined $self->{tasks}{$taskname}{depend} &&
						defined $self->{tasks}{$taskname}{deparray} &&
						(any { defined $tasks2run{$_} } @{$self->{tasks}{$taskname}{depend}})) {
						my @depjobs = grep { defined $tasks2run{$_} } @{$self->{tasks}{$taskname}{depend}};
						if (all { defined $jobids{$_,$rng} } @depjobs) {
							my @depids = map { @{$jobids{$_,$rng}} } @depjobs;
							$depopt = "-hold_jid_ad " . join(',', @depids);
						}
						elsif (all { defined $jobids{$_} } @depjobs) {
							my @depids = map { @{$jobids{$_}} } @depjobs;
							$depopt = "-hold_jid " . join(',', @depids);
						}
						else {
							croak "Cannot find IDs for jobs that $taskname\[$rng\] depends on"
						}
					}
					if ($step > 1) {
						$task_cmd = "$qsub $depopt -t $rng:$step $srcdir/$taskname.sh\n";
					}
					else {
						$task_cmd = "$qsub $depopt -t $rng $srcdir/$taskname.sh\n";
					}
					print $task_cmd;
					unless($arg->{dryrun}) {
						my $stdout = `$task_cmd`;
						if ($stdout) {
							my $jid = _parse_jobid($stdout);
							push @{$jobids{$taskname}} => $jid;
							push @{$jobids{$taskname,$rng}} => $jid;
						}
						else {
							croak "Cannot submit job array $taskname\[$rng\]";
						}
					}
					else {
						my $jid = $counter ++;
						push @{$jobids{$taskname}} => $jid;
						push @{$jobids{$taskname,$rng}} => $jid;
					}
				}
			}
			else {
				# Safeguard
				if ($self->get_num_slots($taskname)) {
					croak "$taskname should not be a job array";
				}
				$qsub .= " -o $logdir/$taskname.out -e $logdir/$taskname.err";
				$task_cmd = "$qsub $srcdir/$taskname.sh\n";
				print $task_cmd;
				unless($arg->{dryrun}) {
					my $stdout = `$task_cmd`;
					if ($stdout) {
						push @{$jobids{$taskname}} => _parse_jobid($stdout);
					}
					else {
						croak "Cannot submit job $taskname";
					}
				}
				else {
					push @{$jobids{$taskname}} => $counter ++;
				}
			}
		}
	}
	return $self;
}

# helper function: parse job ID from qsub output
sub _parse_jobid {
  my ($stdout) = @_;
  my $jobid;
  if ($stdout =~ /^Your job (\d+)/) {
	$jobid = $1;
  }
  elsif ($stdout =~ /^Your job-array (\d+)\.(\d+)\-(\d+):(\d+)/) {
	$jobid = $1;
  }
  else {
	croak "Cannot determine job IDs from $stdout";
  }
  return $jobid;
}

# helper function: recursively identify depended jobs to run
sub _recur_lookup {
	my $self = shift @_;
	my $tasks2run = shift @_;
	my $alljobs = shift @_;

	# Go through each mandatory jobs (i.e., who depends on jobs that will be rerun)
	foreach my $job (@_) {
		my @deps = grep { defined $alljobs->{$_} } $self->{depends}->successors($job);
		if ($self->get_dep_array($job)) {
			foreach my $dep (@deps) {
				if (!defined $tasks2run->{$dep}) {
					$self->_recur_lookup($tasks2run, $alljobs, $dep);
				}
				else {
					# In case of deparray, all successors shuold be job arrays.
					if (ref $tasks2run->{$dep} ne 'Set::IntSpan') {
						croak "Not a Set::IntSpan object";
					}
				}
			}
			my $set = Set::IntSpan->new(map { $tasks2run->{$_}->run_list() } @deps);
			if (defined $tasks2run->{$job}) {
				if ($tasks2run->{$job} lt $set) {
					carp "Run slots for $job should be expanded";
					$tasks2run->{$job} = $set;
				}
				elsif ($tasks2run->{$job} gt $set) {
					carp "Run slots for $job may be incorrect (more than needed)";
				}
			}
			else {
				$tasks2run->{$job} = $set;
			}
		}
		else {
			unless (defined $tasks2run->{$job}) {
				if ($self->get_num_slots($job)) {
					my $nslots = $self->get_num_slots($job);
					$tasks2run->{$job} = Set::IntSpan->new([1..$nslots]);
				}
				else {
					$tasks2run->{$job} = 1;
				}
			}
		}
	}

	return $self;
}


=head1 AUTHOR

Xueya Zhou, C<< <xueyazhou at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-utils at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Utils>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Utils::Workflow


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

1; # End of Utils::Workflow

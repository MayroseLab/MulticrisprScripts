package multicrispr_CONSTS_and_Functions;

use strict;
use warnings;
use File::Slurp;
use File::Path;# qw(make_path);
use File::Basename;


use vars qw(@ISA @EXPORT);

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw(WriteToFile WriteToFileWithTimeStamp ReadFromFile CheckInputTypeValid);



######### Consts ###########
use constant SERVER_NAME 				=> "MultiCRISPR";
use constant RESULTS_DIR_ABSOLUTE_PATH 	=> "/bioseq/data/results/multicrispr/";
use constant RESULTS_LINK 				=> "results";
use constant LOG_FILE 					=> "output.log";
use constant ERROR_STATUS_LOG_FILE 		=> "error.OU";
#use constant OUTPUT_FILES_LIST 			=> "output_files.dat";
use constant OUTPUT_FILE 				=> "output.txt"; # Ofer
use constant LOG_DIR_ABSOLUTE_PATH 		=> "/bioseq/data/results/multicrispr/logs/";
use constant MAIN_PIPELINE_LOG 			=> "multicrispr.log";
#use constant PYTHON_MODULE_TO_LOAD 		=> 'anaconda3-4.0.0';#"python/anaconda_python-3.5"; ##'python/anaconda3-4.0.0';#'python/python-3.3.0';#'python/python-2.7.6';"python/anaconda-python3.5"; anaconda3-4.0.0
use constant PYTHON_MODULE_TO_LOAD 		=> 'python/python-anaconda3.6.5';#"python/anaconda_python-3.5"; ##'python/anaconda3-4.0.0';#'python/python-3.3.0';#'python/python-2.7.6';"python/anaconda-python3.5"; anaconda3-4.0.0
use constant PERL_MODULE_TO_LOAD 		=> 'perl/perl-5.28.1';	#'perl/perl-5.16.3'; (JeKyl)
use constant QSUB_JOB_NUM_FILE			=> "qsub_job_num.dat";
use constant RESULTS_PAGE_URL			=> "/results.html";
use constant RESULTS_PAGE_URL1			=> "/results1.html";
use constant FILENAME_PARAMS			=> "params.html";
use constant SERVER_LINK  => "http://multicrispr.tau.ac.il";
use constant RESULT_TREE  => "tree.newick";
use constant DAILY_TESTS_DIR => "/bioseq/bioSequence_scripts_and_constants/daily_tests/";




######### Functions ############
sub WriteToFileWithTimeStamp
{
	my ($file, $message) = @_;
	
	my $timestamp = localtime();
	
	&WriteToFile($file, "$timestamp: $message");
}

sub WriteToFile
{
	my ($file, $message, $shouldOverwrite) = @_;
	
	# creating file dir if doesnt exist
	my ($fileName, $fileDir) = fileparse($file);
	#make_path($fileDir);
	mkpath($fileDir);
	
	$message =~ s/^\s+//;
	
	if (defined $shouldOverwrite && $shouldOverwrite){
		write_file($file, "$message\n");
	}
	else {
		append_file($file, "$message\n");
	}	
}

sub ReadFromFile
{
	my ($file, $defaultValue) = @_;
	
	if (defined $defaultValue)
	{
		if (-e $file)
		{
			my $line = read_file($file);
			return $line;
		}
		else
		{
			# this is ugly. delete this after renaming all relevant files to dat.
			my ($name, $dir, $ext) = fileparse($file, ".dat");
		 	if ($ext eq ".dat")	{
		 		my $newFile = substr($file, 0, -4);
		 		
		 		if (-e $newFile) {
					my $line = read_file($newFile);
					return $line;
				}
		 	}
			
			return $defaultValue;
		}
	}
	else
	{
		my @lines;
		
		if (-e $file)
		{
			@lines = read_file($file);
		}
		else {
			# this is ugly. delete this after renaming all relevant files to dat.
			my ($name, $dir, $ext) = fileparse($file, ".dat");
		 	if ($ext eq ".dat")	{
		 		my $newFile = substr($file, 0, -4);
		 		
		 		if (-e $newFile) {
					@lines = read_file($newFile);
				}
		 	}
		}
		
		return @lines;
	}
}
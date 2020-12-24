#!/usr/bin/perl
#-------------------------------------------------------------------------#
=begin

Name:			send_commands_to_queue.pl

Last update:	August 2014

Description:	The script receives a commands list file, along with required
				arguments, and send the commands to queue. Works on jekyl and
				lecs2.

Usage:			Should be run from command line, while NOT connected (ssh) to
				any node. A typical running example should look like:
				perl /SVN/scripts/mayrose_lab/send_commands_to_queue.pl /project/cmds_to_queue /temp/ itaym
				* When passing arguments, avoid relative paths such as
				  those including ./ or ../

Arguments:		1. commands file - a file of commands to be sent to queue.
					Each command should be in a separate line, followed by
					a tab, and a job name. For example:
					python /project/my_script.py argument1 argument2	SAMPLE_JOB
					* REQUIRED
				2. temp directory - a location where sh, OU and ER files
					will be saved. Make sure you have writing permission!
					* REQUIRED
				3. queue name - the queue name to send commands to,
					for example: itaym, pupko etc...
					* REQUIRED
				4. time limit - time limit for the job running time, after
					which it will be automatically killed. Time should be
					given in hours!
					* OPTIONAL - default is two weeks.
				5. priority - you can assign lower priority to a job by
					assigning it a negative number (usually just use -1,-2,-3)
					* OPTIONAL - default is 0 (hight priority)

=cut
#-------------------------------------------------------------------------#

use strict;

# read command line args
my $cmds_file=shift; # commands file to run
my $tmp_dir=shift; # where the STDOUT,STDERR,Q_SHELL are created
my $queue=shift;	# queue name
my $time_limit=shift;
my $priority = shift;

# check args and assign defaults
die "Missing argument commands file!\n" unless defined $cmds_file;
die "Missing argument temp directory!\n" unless defined $tmp_dir;
die "Missing argument queue name!\n" unless defined $queue;
if (not defined $time_limit){
	$time_limit = 1209600;
}
else{
	$time_limit = $time_limit*60*60; # convert hours to seconds
}
if (not defined $priority){
	$priority = 0;
}

open (CMDS,$cmds_file) || die "Can't open cmds file: '$cmds_file' $!";
while (my $line=<CMDS>)
{
	chomp ($line);
	my ($cmd,$prefix_name)=split(/\t/,$line);
	chomp ($cmd);

	open (SH_SCRIPT,">",$tmp_dir."/".$prefix_name."\.sh");
	print SH_SCRIPT '#!/bin/bash',"\n";
	print SH_SCRIPT '#$ -N ',"$prefix_name\n";
	print SH_SCRIPT '#$ -S /bin/bash',"\n";
	print SH_SCRIPT '#$ -cwd',"\n";
	print SH_SCRIPT '#$ -e ',$tmp_dir,'/$JOB_NAME.$JOB_ID.ER',"\n";
	print SH_SCRIPT '#$ -o ',$tmp_dir,'/$JOB_NAME.$JOB_ID.OU',"\n";
	print SH_SCRIPT "$cmd\n";
	close(SH_SCRIPT);
	system ("chmod +x $tmp_dir".$prefix_name.".sh");

	system ("qsub -l $queue -l h_rt=$time_limit -p $priority $tmp_dir".$prefix_name.".sh");
}
close (CMDS);
#!/usr/bin/perl -w
#
# TODO:
# perlscript should redirect it's output, especially if running from at
#  (maybe sth like 'tee'?)

use strict;
use Getopt::Long;

use constant VERSION => '$Revision 0.0$';

my $debug=1;
my $mpirun;
my $llstat;
my $schedbgl_cmd;
my $runpath;
if ($debug) {
    $mpirun = "/bin/echo ";
    $schedbgl_cmd="/bin/echo 'Q   1131   $ENV{USER}  2006.04.06_18:00 2006.04.06_19:00 R11'";
    $llstat="echo 'R11 $ENV{USER} init'";
    $runpath = $ENV{HOME};
} else {
    $mpirun = "/bin/echo";
    $schedbgl_cmd="sched_bgl -l";
    $llstat="llstat";
    $runpath = "/home5/hch02/hch026/bglhome/b4.05_L32T64_k0.157025_mu0.003/";    
}
my $atShellScript="atScript.tmp.sh";
# FIXME this should first canonicalize the path...
my $atPerlScript=$0;
my $partition = 'R11 ';
my $mode = 'VN';
my $executable = "/home5/hch02/hch026/bglhome/bin/hmc_tm_xyzt";
my $logfile = "runlog";

# 
#my @startingtimes = (0, 6, 12);
#my @startingdays  = (29, 29, 30);



my ($verbose,$dryrun,$at)=(0,0);
my $getoptRet=GetOptions ( 'help|?' => \&usage,
                           'verbose|v+' => \$verbose,
                           'quiet'   => sub { $verbose = 0 },
                           'dryrun|n!' => \$dryrun,
                           'at' => \$at,
			   'executable' => \$executable,
			   'runpath' => \$runpath,
			   'logfile' => \$logfile,
			   'mode' => \$mode,
                           );
exit -1 unless ($getoptRet);

my ($resid)=@ARGV;
usage() unless (defined($resid));

if ($at) {
    submitAtJob($resid);
} else {
    run($resid)
}
    

#################### submit at script ##########################################
sub submitAtJob {
    my ($resid)=@_;

    my %reservationParameters=bglJobParameters($resid);
    unless (defined($reservationParameters{resid})) {
	printf(STDERR "Job with resID %s does not seem to exist (according to sched_bgl -l output)\n",$resid);
	exit -1;
    }
    if (!($ENV{USER} eq $reservationParameters{user})) {
	printf(STDERR "Wrong username (running under %s, submitting for %s)\n",$ENV{USER},$reservationParameters{user});
	exit -1;
    }
    
    open(ATSCRIPT, "> $atShellScript");
    print ATSCRIPT <<EOF;
#!/bin/sh
$atPerlScript --executable $executable --runpath $runpath --logfile $logfile --mode $mode $resid
EOF
    close ATSCRIPT;
    my $attime=bglTime2atTime($reservationParameters{start});
    if ($dryrun) {
	print 'Would submit at job for '.$attime."\n";
    } else {
	system "at -f $atShellScript $attime";
    }
    
    unlink $atShellScript unless ($debug);
}


############################## run #############################################
sub run {
    my ($resid)=@_;
    
    bglAvailable($ENV{USER},$resid,-1);

    my %reservationParameters=bglJobParameters($resid);
    my $partition=$reservationParameters{partition};
    
    my $command = "$mpirun -partition $partition -mode $mode -exe $executable -cwd $runpath > ".
	"$runpath/$logfile.$reservationParameters{resid}";
    
    waitForBglInitialize(%reservationParameters);
    
    print "Running: $command\n";
    my $errorstatus = system "$command" unless($dryrun);
    
    exit 0;
    ############## we might want to extend here one day ###########################
    
    sleep 60;
    my $status = `llstat | grep $partition`;
    if(  bglAvailable($ENV{USER},$resid) && waitForBglInitialize(%reservationParameters) && bglAvailable($ENV{USER},$resid) ) {
	for(my $j = 0; $j < 2; $j++) {
	    if ($errorstatus != 0) {
		print "Job finished with error, try again!\n";
		$errorstatus = system "$command";
	    }
	}
	if ($errorstatus != 0) {
	    print "Job failed to restart twice, Aborting\n";
	}
    }
}

###################################################################################

sub bglTime2atTime {
    my ($bglTime)=@_;
    # format of bgl time is
    # YYYY.MM.DD_HH:MM
    # format of at time must be
    # HH:MM MMDDYY
    my ($date,$time)=split '_', $bglTime;
    my $year=substr $date, 2,2;
    my $month=substr $date, 5,2;
    my $day=substr $date, 8,2;
    return (sprintf("%s %s.%s.%s",$time,$month,$day,$year));
}

sub bglJobParameters {
    my ($ResID)=@_;
    # output of sched_bgl -l is
    # Q   ResID       User              Start                End             Resource Prio
    open(SCHEDBGL, "$schedbgl_cmd|");
    my %reservationParameters;
    while (<SCHEDBGL>) {
	next unless (/.*\s+$ResID/);
	my (@line)=split;
	$reservationParameters{"resid"}=$line[1];
	$reservationParameters{"user"}=$line[2];
	$reservationParameters{"start"}=$line[3];
	$reservationParameters{"end"}=$line[4];
	$reservationParameters{"partition"}=$line[5]
    }
    close SCHEDBGL;
    return %reservationParameters;
}

sub bglAvailable {
    my ($user,$resid,$exit)=@_;
    my $sched_bgl=`$schedbgl_cmd`;
    if ($sched_bgl =~ /$resid\s+$user/) {
	return(1);
    }
    if (defined($exit) && ($exit != 0)) {
	printf STDERR "vaffanculo: Somebody stole our slot for %s, exiting.\n",$resid;
	exit $exit;
    }
    return;
}

sub waitForReservationTime {
    my ($starttime,$startday)=@_;
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    while( ($hour < $starttime) || ($mday != $startday)) {
	print "sleeping ...\n";
	sleep 60;
	($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    }
}

sub waitForBglInitialize {
    my (%jobPars)=@_;
    my $status = `$llstat | grep $jobPars{partition}`;
    while (($status !~ /$jobPars{user}/) || ($status !~ /init/) ) {
	printf "%s: Waiting for partition to be initialised...\n",$jobPars{resid};
	sleep 60;
	bglAvailable($jobPars{user},$jobPars{resid},-1);
	$status = `$llstat | grep $jobPars{partition}`;
    }
    return 1;
}

sub usage {
    use File::Basename;
    my $basename=basename $0;
    printf <<EndOfUsage ,VERSION;
    $basename [opts] {resID}
   Version: %s
   No extensive help yet. Write it!
    --help                    print this help message
    --verbose                 increase verbosity level (multiple times possible)
    --quiet                   set verbosity to a minimum
    --[no]dryrun, -n          enable (disable) dryrun mode
    --at                      start at script (otherwise start directly)
    --executable              full path to executable
    --runpath                 path in which directory to run
    --logfile                 logfile
    --mode                    mode of the blue gene

EndOfUsage
    exit(0);
}

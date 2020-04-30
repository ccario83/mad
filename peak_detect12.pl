#----------------------------( peak_detect12 )---------------------------------
#    Copyright 2007 by Clinton Cario.  All Rights Reserved.                                                                                                                                                                          
#--------------------------------------------------------------------------------
#

#! usr/bin/perl  
use Cwd;
use Bio::Trace::ABIF;
use GD::Graph::lines; 
use GD::Graph::bars; 
use Spreadsheet::WriteExcel;
#use strict;
use warnings;
select((select(STDOUT), $|=1)[0]);	#to autoflush the write buffer   

#==================
# Set configuration variables 
#=====================
my @conf_values = @{&read_conf_file()};
#for trace cleaning
my $low_peak_thresh = $conf_values[0]; 
my $std_low_peak_thresh = $conf_values[1];
#parameters for gaussian smooth
my $smooth_sigma = $conf_values[2];			
my $mask_size = $conf_values[3];	
#for average smooth
my $average_smooth_width = $conf_values[4];			
#for BP trace calibration
my $std_error_thresh = $conf_values[5];
my $num_before_fail = $conf_values[6];
my $avg_first_peak_loc = $conf_values[7];
my $first_50_offset = $conf_values[8];
my $second_50_offset = $conf_values[9];
#for allele discrimination
my $allele_peak_low_threshold = $conf_values[10];
my $allele_peak_min_percent = $conf_values[11];
#for frequency analysis and display 
my $create_images = $conf_values[12];
my @total_count = (0) x 1000; 

#==================
# User Input Variables
#=====================
my $directory;
my @file_list;
my $standard_dye = "ZZZ";
my $analyzed_dye = "ZZZ";
my $lower_expected = -1;
my $upper_expected = -1;
my $database_location;  
my $marker;
my $call_alleles = "ZZZ";

#==================
# Get Input and Load Information
#=====================
$marker  = &promptUser("\n\nEnter the analyzed marker (as appears in the database)");
$directory  = &promptUser("Enter the location of your *.fsa files");
	if (chdir $directory) { print "\nLoading File information..."; @file_list = <*.fsa>; } 
	else { die "\nDirectory and/or files not found!"; }
	#Open the first ABIF file for trace information
	my $abif = Bio::Trace::ABIF->new();
	$abif->open_abif($file_list[0]);
	if (!$abif->is_abif_format()) { print "\n$file_list[0] is invalid!"; exit; }
	my @dye;
	my %trace;
	for (1..4)
	{ 	
		$dye[$_] = $abif->dye_name($_);
		$trace{ $dye[$_] } = [ $abif->raw_data_for_channel($_) ];
		print "\nLoaded trace $_ for $dye[$_]"; 
	}
	$abif->close_abif();
print "\n\n";
while ($standard_dye ne $dye[1] and $standard_dye ne $dye[2] and $standard_dye ne $dye[3] and $standard_dye ne $dye[4])
	{ $standard_dye = &promptUser("Which dye is your MW standard marker? ", "ROX"); }
while ($analyzed_dye ne $dye[1] and $analyzed_dye ne $dye[2] and $analyzed_dye ne $dye[3] and $analyzed_dye ne $dye[4])
	{ $analyzed_dye = &promptUser("Which dye is your sample label? ", "HEX"); }
while ($lower_expected < 0 or $lower_expected > 500)
	{ $lower_expected = &promptUser("What is the lower bound of the expected allele range? ", 0); }
while ($upper_expected < 0 or $upper_expected > 500)
	{ $upper_expected = &promptUser("What is the upper bound of the expected allele range? ", 500); }
while ($call_alleles ne "Y" and $call_alleles ne "N")
{$call_alleles = &promptUser("Do you want allele calls to be automatically made (Y or N)? ", "Y");}

#==================
# Open and setup output files
#=====================
#error log
unless (open (ERROR, '>', cwd()."/".$marker."_errors.log"))						
	{die "Cannot create file: $!";}
print ERROR "\nError Log";
#excel sheet
my $row = 4;
my $workbook = Spreadsheet::WriteExcel->new(cwd()."/".$marker.".xls");
my $worksheet = $workbook->add_worksheet();
$worksheet->write(1, 2, "Marker:");
$worksheet->write(1, 3, "$marker");
$worksheet->write(3, 1, "Extraction number");
$worksheet->set_column('B:B', 17);
$worksheet->write(3, 2, "Filename");
$worksheet->set_column('C:C', 17);
$worksheet->write(3, 3, "Sample");
$worksheet->set_column('D:D', 7);
$worksheet->write(3, 4, "Warnings");
$worksheet->set_column('E:E', 9);
$worksheet->write(3, 5, "Alleles");
$worksheet->set_column('F:F', 8);
$worksheet->set_column('G:G', 8);
$worksheet->set_column('H:H', 8);
my $format;

#==============
# File Analysis
#=================
my $sample_num = 0;
my $current_file;
my @BP_trace;
foreach $current_file (@file_list)
{
	#LOAD AND CHECK: get the next file, check its validity and load the traces
	#-------------------------------------------------------------------------
	$sample_num++;
	my $abif = Bio::Trace::ABIF->new();
	$abif->open_abif($current_file);
	if (!$abif->is_abif_format())
		{ print "\n$current_file is invalid!"; next;	}
	my $sample_name = $abif->gene_scan_sample_name();
	$worksheet->write($sample_num+3, 2, $current_file);
	$worksheet->write($sample_num+3, 3, $sample_name);
	if ($sample_name eq  "a")
	{ 
		#print "\n$current_file is designated as a blank, skipping";
		$sample_num --;
		next; 
	}
	#print "\n_______________________________________________";
	print "\nAnalyzing: $current_file"; 
	#print "      Sample name: $sample_name";
	print ERROR "\n_______________________________________________";
	print ERROR "\n\n\n$current_file"; 
	print ERROR "\tSample: $sample_name";
	my @dye;
	my %trace;
	for (1..4)
	{ 	
		$dye[$_] = $abif->dye_name($_);
		$trace{ $dye[$_] } = [ $abif->raw_data_for_channel($_) ];
	}

	#CLEAN TRACES: Threshold the traces for low values and slopes, then apply average and gaussian filters
	#-----------------------------------------------------------------------------------------------------
	my @raw_std;
	my @raw_samp;
	my $std_baseline;
	foreach my $value (1..20)
	{
		$std_baseline += @{$trace{$standard_dye}}[$value];
	}
	#get std baseline and threshold, store raw data
	$std_baseline = int ($std_baseline/20);
	foreach my $value (@{$trace{$standard_dye}})
	{
		$value -= $std_baseline;
		if ($value < 0)
			{ $value = 0; }
		push @raw_std, $value;
		if ($value < $std_low_peak_thresh)
			{ $value = 0; }
	}
	#get sample baseline and threshold, store raw data
	my $samp_baseline;
	foreach my $value (1..20)
	{
		$samp_baseline += @{$trace{$analyzed_dye}}[$value];
	}
	$samp_baseline = int ($samp_baseline/20);
	foreach my $value (@{$trace{$analyzed_dye}})
	{
		$value -= $samp_baseline;
		if ($value < 0)
			{ $value = 0; }
		push @raw_samp, $value;
		if ($value < $low_peak_thresh)
			{ $value = 0; }
	}
	#smooth both traces
	for ($analyzed_dye, $standard_dye)
	{
		#Average smooth 
		$trace{$_} = &average_smooth($trace{$_}, $average_smooth_width);
		#Gaussian smooth the trace to remove any insignificant peaks
		$trace{$_} = &gauss_smooth($smooth_sigma, $mask_size, $trace{$_});
	}	
	my @s_trace = @{$trace{$standard_dye}};
	my @a_trace = @{$trace{$analyzed_dye}};
	my $length = (scalar @s_trace) -1;
	#get inital peak information
	my @std_peaks = @{&find_peaks(\@s_trace, undef, undef)};
	my @ana_peaks = @{&find_peaks(\@a_trace, 1500, 2000)};

	#CALIBRATE: Use the standard trace to find the BP trace conversion
	#-------------------------------------------------------------------------------------------
	my $trace_ref = &dynamic_calibrate_dbp($sample_num, \@ana_peaks, \@std_peaks, $length, $std_error_thresh, $num_before_fail, $avg_first_peak_loc, $first_50_offset, $second_50_offset);
	if (defined($trace_ref))
	{
		@BP_trace = @{$trace_ref}; 
		if ($create_images)
			{ &graph([ \@BP_trace, \@raw_std,\@raw_samp,\@s_trace, \@a_trace ], cwd()."/".$current_file.".png"); }	}
	else 
	{
		print ERROR "\nERROR: Analysis could not be completed.";
		if ($create_images)
			{ &graph([ [0..$length], \@raw_std,\@raw_samp,\@s_trace, \@a_trace ], cwd()."/".$current_file.".png"); }
		$worksheet->write($sample_num+3 , 4, "Failed");
		next;
	}
	
	#FIND PEAKS: Find the peaks in the analyzed trace using the user specified range
	#-------------------------------------------------
	for (0..$length)
	{
		#print "\n($BP_trace[$_],$lower_expected)";
		if ($BP_trace[$_] < $lower_expected)
		{ $from = $_; }
		if ($BP_trace[$_] > $upper_expected)
		{ $to = $_; last;}
	}
	@ana_peaks = sort down_by_height @{&find_peaks(\@a_trace, $from, $to)};
	
	#MAKE ALLELE CALLS: If specified, remove any peaks that are not alleles as decided by the configuration file
	#-------------------------------------------------
	if ($call_alleles eq "Y")
		{ @ana_peaks = @{&make_allele_calls(\@ana_peaks, $allele_peak_low_threshold, $allele_peak_min_percent)}; }

	#OUTPUT: Update error log and the excel output sheet
	#-------------------------------------------------
	print ERROR "\nAllele calls in requested range:"; 

	foreach my $num (0..4)
	{
		if (defined(@{$ana_peaks[$num]}[3]))
		{
			my $loc = @{$ana_peaks[$num]}[3];
			my $height = @{$ana_peaks[$num]}[2];
			#my $bp_conv = int ($BP_trace[$loc]+0.5);
			my $bp_conv = $BP_trace[$loc];
			#makes it a bit more accurate (traces shifted a bit by smooth)
			$bp_conv -= 0.2;
			my @standards = (50,75,100,125,150,200,250,300,350,400,450,475,500);
			foreach (@standards)
			{
				if ( abs(int($bp_conv)-$_) <= 1)
				{
					print ERROR "\n\tWARNING: NEXT PEAK MAY BE A STANDARD BLEED. Verify by hand.";
					$format = $workbook->add_format(bold => 1);
				}
				else 
				{
					$format = $workbook->add_format(bold => 0);
				}
			}
			$bp_conv = int($bp_conv*10)/10;
			print ERROR "\n\t$bp_conv BP, height of $height"; 
			$worksheet->write($sample_num+3, $num+5, $bp_conv, $format);
			$total_count[int($bp_conv+0.5)]++;
		}
	}
	#close the current file
	$abif->close_abif();
}

#==============
# Frequency Graph
#=================
&graph_freq([ [0..500], \@total_count ], cwd()."/".$marker."_frequencies.png");

#==============
# Cleanup and say goodbye
#=================
#close the error log and workbook
close ERROR;
 $workbook->close();
#print "\n_______________________________________________";
print "\n\nAnalysis Completed Successfully.";
print "\n3 files were generated in your sample directory: ";
print "\n\t1) $marker.xls - An excel file with allele information.";
print "\n\t2) $marker"."_errors.log - A log file with calibration information and any\n\t    encountered warnings.";
print "\n\t3) frequency.png - An image showing allele frequencies encountered.";
print "\nIndividual trace images were also created if specified by the config file.\n";
print "\n\n\t\t(Press Enter to close this window...)";
<STDIN>;

#==============
# Functions
#=================
	#find_peaks: the body of code that detects peaks
	#-------------------------------------------------
	#INPUT: a cleaned trace reference. OPTIONAL: a range specified by the start and stop x values
	#OUTPUT: a reference to an array of peaks, each containing [ peak start, peak stop, max peak height, center location]
sub find_peaks()
{
	#input
	my ($trace_ref, $from, $to) = @_;
	my @trace= @{$trace_ref};
	my $length = (scalar @trace) -1;
	my @peaks;
	#check range or initialize it
	if (!defined($from) or $from-150 <= 1)
		{ $from = 150; }
	if (!defined($to) or $to+150 >= $length)
		{ $to = $length-150; }
	#print "\n\nSearch for peaks from $from to $to.";
	#state variables
	my $state = "out";
	my $height = $trace[$from-1];
	my $last_height;
	my $pos;
	#peak information variables
	my $start;
	my $peak_height;
	my $center;
	
	#main loop
	for ($pos = $from-150; $pos < $to+150; $pos++)
	{
		#get the new height value, check it, and then round to the tenth place
		$last_height = $height;	
		$height = $trace[$pos];
		if (!defined($height))
			{ print ERROR "\n\tERROR: BAD TRACE VALUE at $pos. Analysis failed."; }
		$height = int($height*10)/10;
		#print "\n$state-$pos: $height";
			#if the state is outside of a peak
			if ($state eq "out")
			{
				if ($height > $last_height)
				{
					$state = "inc";
					$start = $pos;
					next;
				}
			}
			#if the state is inside a peak and increasing
			if ($state eq "inc")
			{
				$center = $pos;
				$peak_height = $height;
				if ($height == $last_height)
					{ $state = "flat"; next;}
				if ($height < $last_height)
					{ $state = "dec"; next;}
			}
			#if the state is inside a peak and decreasing
			if ($state eq "dec")
			{
				#no change in height, we have a peak, and are now outside a peak
				if ($height == $last_height)
				{					
						#print "\n\tPeak ended (too flat): [$start, $pos-10]";					
						push @peaks, [ $start, $pos-10, $peak_height, int($center) ];
						$start = $pos;
						$peak_height = 0;
						$center = 0;
						$state = "out";
						next;
				}
				#if the height has increased, we have a peak, but are still increasing in a new peak
				if ($height > $last_height)
				{
					#print "\n\tDEC=>INC Peak found: [$start, $pos]";					
					push @peaks, [ $start, $pos, $peak_height, int($center) ];
					$start = $pos;
					$peak_height = 0;
					$center = 0;
					$state = "inc";
					next;
				}
			}
			#if the state is inside a peak but no change in height has occured
			if ($state eq "flat")
			{
				#still now change in height, we are now outside a peak
				if ($height == $last_height)
				{
						#print "\n\tPeak ended (too flat)...not recorded\n";
						$start = 0;
						$peak_height = 0;
						$state = "out";
						next;
				}
				#if height has increased, enter the increasing state
				if ($height > $last_height)
				{
					$state = "inc";
					next;
				}
				#if height has decreased, enter the decrease state
				if ($height < $last_height)
					{ $state = "dec"; }
			}
	}		

#check to see if each peak's center is in range and store it to return
my @in_range;
foreach (@peaks)
{
	if (@{$_}[3] >= $from and @{$_}[3] <= $to)
		{ push @in_range, $_; }
}
return \@in_range;
}


	#gauss_smooth: a function to gaussian smooth a trace
	#-------------------------------------------------
	#INPUT: a sigma value, mask size, and reference to a trace
	#OUTPUT: a reference to a smoothed trace 
sub gauss_smooth()
{
	#input
	my ($sigma, $mask_size, $array_ref) = @_;
	my @trace = @{ $array_ref };
	my @filter;	
	my $scale = 0;
	my $half_mask = ($mask_size-1)/2;
	#build filter mask
	for (my $pos = 0; $pos <= $mask_size-1; $pos++)
	{
		my $x = $pos - $half_mask;
		my $weight = (1/($sigma*2.5066))*exp(-1*(($x)**2)/(2*($sigma**2)));
		$filter[$pos] = $weight;
		$scale += $weight;
	}
	$scale = 1/$scale;
	foreach (@filter)
		{ $_ *= $scale; }
	#apply the filter to each value; 
	for (my $pos = $half_mask; $pos <= (scalar @trace)-$half_mask-1; $pos++)
	{ 
		my @view = @trace[($pos-$half_mask)..($pos+$half_mask)];
		$trace[$pos] = &PMAS(\@filter, \@view);
	}
return \@trace;
}


	#PMAS:  Pairwise Multiply And Sum
	#-------------------------------------------------
	#INPUT: two array references representing arrays of the same size
	#OUTPUT: the sum of the pairwise multiplication of both arrays
sub PMAS()
{
	#input
	my ($ref1, $ref2) = @_;
	my @a1 = @{$ref1};
	my @a2 = @{$ref2};
	my $sum = 0;
	for (my $i = 0; $i <= (scalar @a1)-1; $i++)
		{ $sum += ($a1[$i]*$a2[$i]); }
return $sum;
}

	#graph: graph a trace file 
	#-------------------------------------------------
	#INPUT: a reference to an array containing references to arrays of x-values, and then to any number of y-values, and a filename
	#OUTPUT: a png file with the given filename
sub graph
{	
	#input
	my ($t_ref, $filename) = @_;
	my @data = @{$t_ref};
	#get the maximum height of the first set of y-values for the y-axis (using the first is arbitrary)
	my $max_height = 0;
	foreach (@{$data[1]})
	{
		if ($_ > $max_height)
		{
			$max_height = $_;
		}
	}
	#set the graph parameters
	my $graph = GD::Graph::lines->new(10000, 800);
	$graph->set
	( 
	     transparent       => 0,
	     x_label           => 'Time',
 	     y_label           => 'Intensity',
 	     title             => 'Traces',
 	     x_label_skip      => 100,
	     x_tick_offset     => 50,
	     y_max_value      => $max_height+1000,
             y_tick_number     => 20,
	     y_label_skip      => 2, 
	     dclrs 		=> [ qw(gray gray red blue black black) ]
	)or die $graph->error;
	#plot and create the file
	my $gd = $graph->plot(\@data) or die $graph->error;
	open(IMG, ">$filename") or die $!;
	binmode IMG;
	print IMG $gd->png;
	close IMG;
}
	
	
	#graph_freq: graph an allele frequency 
	#-------------------------------------------------
	#INPUT: a reference to an array containing references to arrays of x-values, and then to y-values
	#OUTPUT: a png file with the given filename
sub graph_freq
{
	my ($t_ref, $filename) = @_;
	my @data = @{$t_ref};
	#get the maximum height of the y-axis
	my $max_height = 0;
	foreach (@{$data[1]})
	{
		if ($_ > $max_height)
		{
			$max_height = $_;
		}
	}
	#set the graph parameters
	my $graph = GD::Graph::bars->new(10000, 800);
	$graph->set
	( 
	     transparent		=> 0,
	     x_label			=> 'Base Pair',
 	     y_label			=> 'Frequency',
 	     title			=> 'Allele BP Frequencies',
	     x_label_skip		=> 5,
	     y_tick_number	=> int($max_height),
	     y_max_value		=> int($max_height),
	     y_label_skip		=> 2,
	     dclrs 			=> [ qw(orange red green black) ],
	)or die $graph->error;
	#plot and create the file
	my $gd = $graph->plot(\@data) or die $graph->error;
	open(IMG, ">$filename") or die $!;
	binmode IMG;
	print IMG $gd->png;
	close IMG;
}


	#average_smooth: smooth a trace by averaging x values together
	#-------------------------------------------------
	#INPUT: a reference to a trace, and the number of points to average together
	#OUTPUT: a reference to a trace whose ever x points have been averaged and set to that average
sub average_smooth
{
	my ($t_ref, $size) = @_;
	my @trace = @{$t_ref};
	my @weights = (1) x $size;	
	
	my @smoothed_trace;
	for (my $pos = 1; $pos < (scalar @trace)-$size+1; $pos += $size)
	{ 
		my @slice = @trace[$pos..$pos+$size-1];
		my $val = &PMAS(\@weights, \@slice);
		$val = int($val/$size);
		for (1..$size)
			{ push @smoothed_trace, $val; }
	}
return \@smoothed_trace;
}


	#Sorting subroutine: down_by_height
	#-------------------------------------------------
	#INPUT: a array of peak references
	#OUTPUT: that array sorted descendingly by the peak heights 
sub down_by_height
	{ @{$b}[2] <=> @{$a}[2]; }


	#dynamic_calibrate_dbp: generate an array of BP values that corresponds to the x-values of the trace file
	#-------------------------------------------------
	#INPUT: The current sample number, a reference to the sample peaks, a reference to the standard peaks, the length of the trace, the standard error threshold, the number of above threshold
	#	     errors before failure, the average first peak location, the average location of the first of a 50bp peak pair, and the average location of the second of a 50bp pair 
	#OUTPUT: undef if the BP values cannot be generated, otherwise a reference to an array of corresponding BP values to the x-values
sub dynamic_calibrate_dbp
{
	#input and vars
	my @warnings;
	my ($sample_num, $sample_peaks_ref, $std_peaks_ref, $length, $std_error_thresh, $num_before_fail, $avg_first_peak_loc, $first_50_offset, $second_50_offset) = @_;
	my @ana_peaks = @{$sample_peaks_ref};
	my @standards = (0,50,75,100,125,150,200,250,300,350,400,450,500);
	my $last_dbpdx = 0.0765;
	my $plusminus = 163;
	my $zeroBP_peak;
	#get the first and largest sample peaks
	my $first_sample_peak = $ana_peaks[0][3]; 
	@ana_peaks = sort down_by_height @ana_peaks;
	my $largest_sample_peak = $ana_peaks[0][3];
	#get the first standard peak, and the closest sample peak to it
	my $first_std_peak = &find_nearest_peak($std_peaks_ref, $avg_first_peak_loc);
	my $closest_to_std =  &find_nearest_peak($sample_peaks_ref, $first_std_peak);

	#check if they are all defined, update the error log, and attempt to find a suitable location for the 0BP peak
	if (defined($first_std_peak) and defined($largest_sample_peak) and defined($closest_to_std))
	{
		print ERROR "\n\tSTD: $first_std_peak LARGE: $largest_sample_peak CLOSEST: $closest_to_std FIRST SAMP: $first_sample_peak";
		if (abs($first_std_peak - $largest_sample_peak) < 25)
			{ $zeroBP_peak = $first_std_peak; }
		elsif (abs($first_std_peak - $closest_to_std) < 10)
			{ $zeroBP_peak = $first_std_peak; }
		elsif ($first_std_peak < $largest_sample_peak)
			{ $zeroBP_peak = $first_std_peak; }
		else
		{
			push @warnings, "\nWarning: ODD STANDARD 0BP PEAK. The standard trace appears to be missing a 0BP calibration peak. Using the sample's 0BP.";
			if (abs($avg_first_peak_loc-$largest_sample_peak) <= abs($avg_first_peak_loc-$first_sample_peak))
				{ $zeroBP_peak = $largest_sample_peak; }
			else
				{ $zeroBP_peak = $first_sample_peak; }
		}
	}
	#otherwise, use the default value for the 0BP peak (and warn)
	else 
	{
		push @warnings, "\nWarning: NO 0BP CALIBRATION PEAK in sample or standard. Using the default value ($avg_first_peak_loc) for the first peak."; 
		$zeroBP_peak = $avg_first_peak_loc;
	}
	
	#get the actual location of the first 50BP peak
	my $fifty_peak_1 = &find_nearest_peak($std_peaks_ref, $zeroBP_peak+$first_50_offset);
	if (!defined($fifty_peak_1))
		{ $fifty_peak_1 = $zeroBP_peak+$first_50_offset; }
	#get the actual location of the second 50BP peak
	my $fifty_peak_2 = &find_nearest_peak($std_peaks_ref, $fifty_peak_1+$second_50_offset);
	if (!defined($fifty_peak_2))
		{ $fifty_peak_2 = $fifty_peak_1+$second_50_offset; }
	# if they are the same, attempt to find the next 50BP peak by searching gradualy searching farther from the first 50BP peak
	my $offset = 0;
	while ($fifty_peak_1 == $fifty_peak_2)
	{
		$offset += 50;
		my $now_look = $fifty_peak_1+$second_50_offset+$offset;
		$fifty_peak_2 = &find_nearest_peak($std_peaks_ref, $now_look);
		if (!defined($fifty_peak_2))
		{
			print ERROR "\nERROR: ODDLY SPACED STANDARDS. 50bp spaced peaks cannot be found. Analysis failed.";
			return undef;
		}
	}
	#if this fails, warn, but continue with the default value
	if (!defined($fifty_peak_1) or !defined($fifty_peak_2))
	{
		push @warnings, "\nWarning: ODDLY SPACED STANDARDS. No 50bp spaced peaks found. Continuing anyway.";
		$last_dbpdx = 0.0765;
	}
	#otherwise conpute the initial dbpdx
	else
	{ $last_dbpdx = 50/($fifty_peak_2-$fifty_peak_1); }
	#update the error log
	print ERROR "\n\t0BP Peak:\t\t$zeroBP_peak.\n\tstarting dBPdX:\t\t$last_dbpdx.\n\t50bps peaks:\t\t$fifty_peak_1   $fifty_peak_2";
	#if the dbpdx value is obsurd, set it back to the default
	if ($last_dbpdx > 0.09 or $last_dbpdx < 0.05)
	{ 
		push @warnings, "\nWarning: ODDLY SPACED STANDARDS. Peaks too close or far apart in expected range. Using the default value for dbpdx";
		$last_dbpdx = 0.0765;
	}
	
	#peak calibration
	my @bp_trace;
	my $last_loc = $zeroBP_peak;
	my $last_bp_std = 0;
	my $value = 0;
	my $skip_count = 0;
	print ERROR "\n\tStandard Peak Information:\n\t\tBP\tExpected\tActual\t\tError\tdBPdX\n\t\t--\t--------\t------\t\t-----\t-----";
	#set all points before the 0BP peak to 0BP
	for (1..$zeroBP_peak)
	{
		push @bp_trace, 0;
	}
	#calibrate for each standard
	foreach (@standards)
	{
		#the 0BP has already been calibrated, goto the next one
		if ($_==0) 
			{ next; }
		#get the expected x-cord for this standard peak
		my $expected = int ($last_loc + (($_-$last_bp_std)/$last_dbpdx));
		print ERROR "\n\t\t$_\t$expected";
		#and the closest actual peak
		my $best_peak_loc = &find_nearest_peak($std_peaks_ref, $expected);
		if (!defined($best_peak_loc))
		{
			push @warnings, "\nWarning: CANNOT FIND the $_ BP STANDARD. Skipping.";
			$best_peak_loc = 0;
		}
		#compute the error
		my $error = abs($expected - $best_peak_loc);
		print ERROR "\t\t$best_peak_loc\t\t$error";
		#if the best peak is far from expected, generate an error up to the threshold, and then fail
		if ($error > $std_error_thresh)
		{ 
			print ERROR "*"; 
			if ($_ == 0)
				{ print " (no zero peak in standard)"; $best_peak_loc = $expected; }
			if ($skip_count++ > $num_before_fail)
			{
				print ERROR "\nERROR: ODDLY SPACED STANDARDS. Did not find standard peaks where expected. Analysis failed.";
				return undef;
			}
			else
				{ next; }
		}
		#compute the dbpdx for this peak
		my $dx = $best_peak_loc - $last_loc;
		if ($dx == 0)
			{ next; }
		my $dbp = $_ - $last_bp_std;
		$dbpdx = $dbp/$dx;
		print ERROR "\t$dbpdx";
		#and update the trace
		for (1..$dx)
		{
			push @bp_trace, $value+=$dbpdx;
		}
		#use this peak for the next calibration
		$last_loc = $best_peak_loc;
		$last_bp_std = $_;
	}

	#update the error log with all warnings so far
	print ERROR @warnings;
	if (defined($warnings[0]))
	{
		my $num=0;
		foreach (@warnings)
			{ $num++; }
		$worksheet->write($sample_num+3 , 4, $num);
	}
	#make sure the bp_trace is defined
	if (!defined($bp_trace[0]))
		{ return undef; }
	#and fill all values > 500BP with 500BP if necessary
	for ((scalar @bp_trace)..$length)
		{ push @bp_trace, $standards[-1]; }
#and return
return \@bp_trace;
}


#----------------------------(  promptUser  )-----------------------------------
#    Copyright 1998 by DevDaily Interactive, Inc.  All Rights Reserved.                                                                                                                                                                          
#  PURPOSE:	Prompt the user for some type of input, and return the    
#		input back to the calling program.                                                                           
#  ARGS:	$promptString - what you want to prompt the user with     
#		$defaultValue - (optional) a default value for the prompt                                                                        
#--------------------------------------------------------------------------------
sub promptUser()
{
	local($promptString,$defaultValue) = @_;
	if (defined($defaultValue)) 
	{
		print $promptString, "[$defaultValue]: ";
	} 
	else 
		{ print $promptString, ": "; }
	$| = 1;			# force a flush after our print
	$_ = <STDIN>;		# get the input from STDIN (presumably the keyboard)
	chomp;
	if ($_ ne	"") 
	{
		return $_;		# return $_ if it has a value
	} 
	else 
		{ return $defaultValue; }
}


	#make_allele_calls: automatically make allele calls using information in the configuration file
	#-------------------------------------------------
	#INPUT: a reference to a peak array, the minimum threshold, and the minimum height percent as compared to the largest allele
	#OUTPUT: a reference to a peak array containing peaks with heights larger than the min threshold, and that are at least 60% of the height of the largest allele peak
sub make_allele_calls()
{
	my ($peaks_ref, $min_thresh, $min_percent) = @_;
	my @peaks = @{$peaks_ref};
	my @called_peaks;
	
	#get the largest peaks height
	my $height = @{$peaks[0]}[2];
	
	if (!defined($height))
		{ return $peaks_ref; }
	else
		{ push @called_peaks, $peaks[0]; }
	#add all peaks that meet the criteria
	for (1..2)
	{
		if (defined(@{$peaks[$_]}[2]))
		{
			my $height2 = @{$peaks[$_]}[2];
			if ( (($height2/$height) >= $min_percent) && $height2 >= $min_thresh)
				{ push @called_peaks, $peaks[$_]; }
		}
	}
	return \@called_peaks;
}
		
	
	#find_nearest_peak: returns the location of the nearest peak to specified value
	#-------------------------------------------------
	#INPUT: a reference to a peak array and a location x-value
	#OUTPUT: the location of the closest peak to the given location 
sub find_nearest_peak()
{
	my ($peaks_ref, $location) = @_;
	my @peaks = @{$peaks_ref};

	my $best_diff = 10000;
	my $best_peak_loc;
	foreach (@peaks)
	{
		my $difference = abs(@{$_}[3] - $location);
		if ($difference < $best_diff)
		{
			$best_diff = $difference;
			$best_peak_loc = @{$_}[3];
		}
	}
return $best_peak_loc;
}


	#read_conf_file: attemps to load all non # lines into an array, or loads default values instead
	#-------------------------------------------------
	#INPUT: a file called PeakDetect.conf containing 13 ordered configuration values and anynumber of lines starting with the '#'character
	#OUTPUT: a reference to an array of configuration values in order
sub read_conf_file()
{
	my @conf_values;
	unless (open (CONF, '+<', "PeakDetect.conf"))
	{
		print "\nCannot find file PeakDetect.conf. Using default values.";
		@conf_values = (300,50,5,31,2,200,2,1680,2650,650,400,0.6,0);
		#$low_peak_thresh = 300;
		#$std_low_peak_thresh = 50;
		#$smooth_sigma = 5;			
		#$mask_size = 31;	
		#$average_smooth_width = 2;			
		#$std_error_thresh = 200;
		#$num_before_fail = 2;
		#$avg_first_peak_loc = 1680;
		#$first_50_offset = 2650;
		#$second_50_offset = 650;
		#$allele_peak_low_threshold = 400;
		#$allele_peak_min_percent = .60;
		return \@conf_values;
	}

	while (<CONF>)
	{
		chomp($_); 
		if (!($_ =~ m/#/))
			{ push @conf_values, $_; }
	}
print "\nConfiguration file \"PeakDetect.conf\" successfully loaded.";
return \@conf_values;
}

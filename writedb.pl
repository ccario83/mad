use warnings;
use Cwd;
use Win32::OLE qw(in with);
use Win32::OLE::Const 'Microsoft Excel';
$Win32::OLE::Warn = 3;                                # die on errors...
select((select(STDOUT), $|=1)[0]);	#to autoflush the write buffer

my $location = cwd();
my $information_location = &promptUser("Where is the excel file containing the information to add? ");
if (!($information_location =~ m/:/))
	{ $information_location = $location."/".$information_location; }
print "\n$information_location\n";
my $database_location = &promptUser("Enter the location of Raul's database ", "Z:/McPheronLab/Anastrepha Thesis 2005-2008/Obliqua Veracruz phylogeography/Microsatellite Work/Data/Main Data/Obliqua microsat main data set.xls");
my $Drow = &promptUser("What row in the database contains marker lables? ", 6);
my $Dcol = &promptUser("What column in the database contains the extraction numbers?", 3);
my $Irow = &promptUser("What row in the information file contains the first extraction number?", 5); 
my $Icol = &promptUser("What column in the information file contains the extraction numbers?", 2);
my $Icol2 = &promptUser("What column in the information file contains the first allele?", 6);
my $Mrow = &promptUser("What row in the information file contains the marker?", 2);
my $Mcol = &promptUser("What column in the information file contains the marker?", 4);

my $Excel = Win32::OLE->GetActiveObject('Excel.Application') || Win32::OLE->new('Excel.Application', 'Quit'); 
my $Excel2 = Win32::OLE->GetActiveObject('Excel.Application') || Win32::OLE->new('Excel.Application', 'Quit'); 
#Open the information file
my $Info = $Excel->Workbooks->Open($information_location); 
my $InfoSheet = $Info->Worksheets(1);
#open the database
my $Db = $Excel2->Workbooks->Open($database_location); 
my $DbSheet = $Db->Worksheets(1);

#get the marker value
my $marker = $InfoSheet->Cells($Mrow,$Mcol)->{'Value'};
if (!defined($marker))
{
	die "\nNo marker fount at $Mrow, $Mcol!";
}

#get the marker column for database writing
my $Dbm_col = -1;
for (1..100)
{
	next unless defined $DbSheet->Cells($Drow,$_)->{'Value'};
	my $value = $DbSheet->Cells($Drow,$_)->{'Value'};
	#print "\n$value";
	if ($value eq $marker)
	{
		$Dbm_col = $_;
		#print "\nmatch for $marker at $Dbm_col, $_";
		last;
	}
}
if ($Dbm_col == -1)
{
	die "\nMarker $marker not found in row $Drow of database";
}


#get the extraction number from info file and corresponding row in database
my $extraction = "Z";
until (!defined($extraction))
{
	#get the extraction number
	$extraction = $InfoSheet->Cells($Irow,$Icol)->{'Value'};
	if (!defined($extraction))
		{ exit; }
		
	print "\nWriting information for extraction: $extraction";	
	
	#get extraction row from database
	my $Dbe_row = -1;
	for (1..800)
	{
		next unless defined $DbSheet->Cells($_,$Dcol)->{'Value'};
		my $value = $DbSheet->Cells($_,$Dcol)->{'Value'};
		if ($value eq $extraction)
		{
			$Dbe_row = $_;
			#print "\nmatch for $extraction at $_, $Dbe_row";
			last;
		}
	}
	if ($Dbe_row == -1)
	{
		die "\nExtraction $extraction not found in database";
	}
	
	#for each allele, write it to the database at Dbe_row, 
	for ($Icol2..$Icol2+2)
	{
		my $allele = $InfoSheet->Cells($Irow,$_)->{'Value'};
		if (!defined($allele))
			{ next; }
		print "\n\tWriting allele: $allele to cell: $Dbe_row,$Dbm_col";
		$DbSheet->Cells($Dbe_row,$Dbm_col+$_-$Icol2)->{'Value'} = $allele;
	}
$Irow++;
}
$Info->close();
$Db->close();


	
#$homeDir  = &promptUser("Enter the home directory ", "/home/$username");

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
	$| = 1;               # force a flush after our print
	$_ = <STDIN>;         # get the input from STDIN (presumably the keyboard)
	chomp;
	if ($_ ne	"") 
	{
		return $_;    # return $_ if it has a value
	} 
	else 
		{ return $defaultValue; }
}



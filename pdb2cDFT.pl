#! /usr/bin/perl

#                        ( ~!~ ~!~ )
#  .-----------------oOOo----(_)----oOOo----------------------------.
#  |                                                                |
#  |                      A Perl Script for:                        |
#  |                                                                |
#  | ..--** Converting pdb files into G09 input files        **--.. |
#  | ..--** for the calculation of Conceptual DFT parameters **--.. |
#  |                                                                |
#  | ..--** Sister script of cDFTExtract.pl                  **--.. |
#  |                                                                |
#  | Version 0.1 - First Version                                    |
#  | Author David Hose                                              |
#  | (c) July 2014                                                  |
#  |                      .oooO                                     |
#  |                      (   )   Oooo.                             |
#  .-----------------------\ (----(   )-----------------------------.
#                           \_)    ) /
#                                 (_/

# VARIABLES:
	# User varialbles:
		my	$Memory = "200MW";	# Memory allocation.
		my	$Processors = 4;	# Maximum number of processor to be used in parallel excution.
		my	$G09Ext = "com";	# The extension to be used for the G09 files.
	# General Variables:
		my	@FileList;			# Holds the list of pdb files to be processed.
		my	@Temp;				# An array that holds the contents of the PDB file current line.
	# Charge/Population flags and methods:
		@Population = (0,0,0);	# An array which holds which population methods are to be used (switches).
		@PopulationMethod = ("npa", "mk", "hirshfeld");	# An array which holds the names of the population methods available.
		$PopString = "pop=(";	# Opening section of the population command.

#	MAIN CODE STARTS HERE:
	Hello();	# Say hello.
	Parser(@ARGV);	# Get the command line input (STDIN).
	# Construct the population cammand depending upon the options selected.
		$cnt = 0;	# A positional counter for the Population Arrays.
		foreach $i (@Population) {
			if ($i == 1) {
				$PopString = $PopString . $PopulationMethod[$cnt] . ",";
			}
			$cnt++; # Increment positional counter.
		}
		$PopString = (substr $PopString, 0, length($PopString) - 1) . ")";	# Trims off the last unwanted comma and closes the bracket.
	# If no population methods have been selected, then default to NPA and Hirshfeld charges.
		if ($PopString eq "pop=)") {$PopString = "pop=(npa,hirshfeld)"};

	print "Processing files...\n";	# Tell the user that the files are being processed.
	# BEGIN of FileList Loop.
	foreach my $PDBFile (@FileList) {
		print "\t$PDBFile ... Creating G09 file.\n";	# Tell the user which file is being processed.
		# Create the new G09 file.
			$G09File = $PDBFile;	# The new filename is based upon the PDB filename.
			$G09File =~ s/.pdb$//i;	# Crude strip off the extension. (Need to anchor this to the end.)
			open(OUTPUT, ">$G09File.$G09Ext") or die "Can't create the G09 file: $!";

		# PART 1:
			# Write the Link0 information:
				print OUTPUT "%mem=$Memory\n";				# Specify the amount of memory.
				print OUTPUT "%nproc=$Processors\n";		# Specify the number of processors.
				print OUTPUT "%chk=$G09File\n";				# Specify the Checkpoint file.
				print OUTPUT "#p opt b3lyp/6-31+G(d) nosymm $PopString\n\n";
				print OUTPUT "$G09File: Parent Compound Optimisation [B3LYP/6-31+G(d)]\n\n";		# Title of calculation including filename.
				print OUTPUT "0 1\n";		# Charge and multiplicity for the neutral closed shell parent.
			# Open the current pdb file and extract the atomic coordinates.
				open(INPUT, "<$PDBFile") or die "Can't open the input PDB file: $!";
				while(<INPUT>)
				{	# BEGIN Reading the PDB coordinates and Writing them to the G09 .com file.
					chomp;
					@Temp = split(/\s+/);
					#print OUTPUT "$Temp[2]\t$Temp[5]\t$Temp[6]\t$Temp[7]\n" if ($Temp[0] eq "HETATM");	# Writes element and coordinate data to the G09 file.
					#print OUTPUT "$Temp[2]\t$Temp[5]\t$Temp[6]\t$Temp[7]\n" if ($Temp[0] eq "HETATM" && $Temp[5] ne "0001");	# Writes element and coordinate data to the G09 file.
					print OUTPUT "$Temp[2]\t$Temp[5]\t$Temp[6]\t$Temp[7]\n" if ($Temp[0] eq "HETATM" && $Temp[4] eq "0001");	# Writes element and coordinate data to the G09 file.
				}	# END Reading and Writing the PDB coordinates to the G09 file.
			# Close the current PDB file.
				close INPUT;
		# PART 2:
			# Link1.
				print OUTPUT "\n--link1--\n";
			# Write the Link0 information:
				print OUTPUT "%mem=$Memory\n";	# Specify the amount of memory.
				print OUTPUT "%nproc=$Processors\n";	# Specify the number of processors.
				print OUTPUT "%chk=$G09File\n";	# Specify the Checkpoint file.
				print OUTPUT "#p rob3lyp/6-31+G(d) scf=tight nosymm guess=read geom=checkpoint $PopString\n\n";	# The calculation routecard.
				print OUTPUT "$G09File: Radical Anion Single Point.\n\n";	# Title of calculation including filename.
				print OUTPUT "-1 2\n";	# Charge and multiplicity for the Open-Shell Radical Anion.
		# PART 3:
			# Link1.
				print OUTPUT "\n--link1--\n";
				# Write the Link0 information:
					print OUTPUT "%mem=$Memory\n";	# Specify the amount of memory.
					print OUTPUT "%nproc=$Processors\n";	# Specify the number of processors.
					print OUTPUT "%chk=$G09File\n";	# Specify the Checkpoint file.
					print OUTPUT "#p rob3lyp/6-31+G(d) scf=tight nosymm guess=read geom=checkpoint $PopString\n\n";	# The calculation routecard.
					print OUTPUT "$G09File: Radical Cation Single Point.\n\n";	# Title of calculation including filename.
					print OUTPUT "1 2\n";	# Charge and multiplicity for the Open-Shell Radical Cation.
		# All parts completed.  Close OUTPUT.
			close OUTPUT;
	}	# END of FileList Loop.
	GoodBye();
	exit;	# Catch all END
#	MAIN CODE ENDS:

# SUBROUTINES.
sub Hello {
print <<ENDHELLO;
			                      ( ~!~ ~!~ )
			.-----------------oOOo----(_)----oOOo----------------------------.
			|                                                                |
			|                      A Perl Script for:                        |
			|                                                                |
			| ..--** Converting pdb files into G09 input files        **--.. |
			| ..--** for the calculation of Conceptual DFT parameters **--.. |
			|                                                                |
			| ..--** Sister script of cDFTExtract.pl                  **--.. |
			|                                                                |
			| Version 0.1 - First Version                                    |
			| Author David Hose                                              |
			| (c) July 2014                                                  |
			|                      .oooO                                     |
			|                      (   )   Oooo.                             |
			.-----------------------\\ (----(   )----------------------------.
			                         \\_)    ) /
			                               (_/


ENDHELLO
}

sub GoodBye {
		print <<ENDGOODBYE;

		Script completed. Goodbye.

ENDGOODBYE
}

sub Usage {
print <<ENDUSAGE;

			             ..--** NMR_pdb2G09 **--..

			This script reads pdb files and sets up a Conceptual DFT calculation which consists of:
				* Geometry optimisation of the Neutral (Parent) species.
				* Single point energy of the Radical Anion.
				* Single point energy of the Radical Cation.

			NMR_pdb2G09.pl [--help] [-npa] [-mk] [-hirshfeld] [-all] File01.pdb File02.pdb

			Optional Switches:

				--help     : This help page.
				-npa       : Requests Natural Population Analysis.
				-mk        : Requests Merz-Singh-Kollman Population Analysis.
				-hirshfeld : Requests Hirshfeld Population Analysis.
				-all       : Requests Natural, Merz-Singh-Kollman and Hirshfeld Population Analysis.

				Multiple .pdb files can be GLOB'd in by using *.pdb.

			Version 0.2	First Version.
			(c) David Hose, July 2014.

ENDUSAGE
exit;
}

sub ParsingError {
	# Prints out error messages from the Parser subroutine.
		print "\nERROR:\t$_[0]\n\n";
		Usage();
}

sub Parser {
	# Parses the Command Line Arguments and sets the variables up for the main routine.
	my @Input = @_;
	# Are there any arguments in the command line? If not, abort with help message.
		Usage() if (scalar(@Input) == 0);
	# Check each input argument.
		foreach $Argument (@Input) {
			# Test for the Help switch (--help).
				Usage() if ($Argument =~ /^--help$/i);
			# Test for which population method is to be run.
				$Population[0] = 1 if ($Argument =~ /^-npa$/i);
				$Population[1] = 1 if ($Argument =~ /^-mk$/i);
				$Population[2] = 1 if ($Argument =~ /^-hirshfeld$/i);
				@Population = (1,1,1) if ($Argument =~ /^-all$/i);
				@Population = (1,0,1) if ($Argument =~ /^-recommended$/i);
			# Test for input file (must be .pdb file).
				push (@FileList, $Argument) if ($Argument =~ /\.pdb$/i);
		}
	# Have any files been parsed?
		ParsingError("No PDB files have been inputed.") if (scalar(@FileList) == 0);
}
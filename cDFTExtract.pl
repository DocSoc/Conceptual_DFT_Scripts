#! /usr/bin/perl

# This script extracts the cDFT data.
# Sister script to pdb2cDFT.pl script.

# VARIABLES:
	# User defined:
		my $SpinThreshold = 10; # Percentage
	
	
	# 
		my @FILELIST; # Holds the names of all of the files to be processed.
		my $FILE; # Holds the name of the current file being processed.
		my $NormTerms = 0; # Holds the number of times "Normal termination" is found in the current file.
		my @Energies = (0,0,0); # Holds the neutral, radical anion and radical anion SCF energies.



		my $au2eV = 27.21138505; # Hartrees (au) to eV convertion factor.
		my $SpinThres = ($SpinThreshold / 100) * 0.7500; # Calculates the acceptable limit of spin contamination for the radical calculations.
			# For a molecular multiplicity of 2, S**2 will be 0.7500.
		my $SpinContaminationCnt = 0; # Counters the number of spin contamination violations in the current file.
		my $SpinContaminationErrorCnt = 0; # Counters the number of spin contamination violations in all files.
		my $Section = 0; # Holds the section number: 0 for Neutral, 1 for Radical Anion, 2 for Radical Cation.

		my @MullikenAtomNum = undef;
		my @MullikenElement = undef;
		my @MullikenNeutral = undef;
		my @MullikenAnion = undef;
		my @MullikenCation = undef;
		
		my @HirshfeldAtomNum = undef;
		my @HirshfeldElement = undef;
		my @HirshfeldNeutral = undef;
		my @HirshfeldAnion = undef;
		my @HirshfeldCation = undef;
		
		my @ESPAtomNum = undef;
		my @ESPElement = undef;
		my @ESPNeutral = undef;
		my @ESPAnion = undef;
		my @ESPCation = undef;


# MAIN CODE STARTS HERE #
Hello();




# parse STDIN for:
	# Filenames (GLOB)
	# Switches that are set.
		# Atoms.
			# "all": All atoms.
			# "heavies": All atoms except hydrogen.
			# "Default": Reactive centre atoms; C, N, S, Si, P,... (exclude halogens, oxygen,...)
	@FileList = ("Pyridine.out"); # Temp.


print "Processing files...\n";	# Tell the user that the files are being processed.
foreach $File (@FileList) {
	print "\t$File\t";	# Tell the user which file is being processed.
	open(INPUT, "<$File") or die "Can't open the G09 file: $!";

	# START WHILE loop to read each line of the current file being processed.
	LINE: while (<INPUT>) {
		$Line = $_;
		chomp($Line);
		# Find "SCF Done" energies and store the last three energies in @Energies array.
			# If there are no errors, this will contain the neutral, radical anion and radical anion SCF energies.
			if($Line =~ /SCF Done/) {
				my @Temp = split(/[\s]+/, $Line); # Splits the line by spaces.
				shift(@Energies); # Discard the unwanted energy value (from the LHS).
				push @Energies, $Temp[5]; # Add the new energy value (to the RHS).
				next LINE;
			}

		# Count the number of "Normal termination" in the file.  There should be 3. If less than this, issue an error message.
			if($Line =~ /Normal termination/) {
				print "Found Normal Termination\n";
				$NormTerm++; # Increment the Normal termination counter (also indicates which section of the output is being processed).
				next LINE;
			}

		# Check for spin contamination of the radicals.  If outside of a specified range, issue a warning (see below).
			if($Line =~ /S\*\*2 before annihilation/) {
				my @Temp = split(/[\s]+/, $Line); # Splits the line by spaces.
				if(abs($Temp[6] - 0.75) > $SpinThres) {
					$SpinContaminationCnt++; # Increment the spin contamination counter (file).
					$SpinContaminationErrorCnt++; # Increment the spin contamination counter (files).
				}
				next LINE;
			}

		# Check if the geometry has been optimised (to ensure that the initial charge calculations are ignored).
			if($Line =~ /^\s+-- Stationary point found/) {
				$Stationary = 1; # Set the Stationary Point Found Flag.
				next LINE;
			}
			
		# Extract charges
		 	# Only extract charges for the required elements.
			# Charge for Neutral, Radical Anion and Radical Cation species.
			# Set Flag that Charges are available for charge method for further calculation and output.

			# Find and extract the Mulliken charges.
			{
				# Find the start of the Mulliken charge section.
					if ($Line =~ /^\sMulliken atomic charges/ && $Stationary == 1) {
						$MullikenFlag = 1; # In the Mulliken charges section of the output.
						$MullikenCharges = 1; # Flag that Mulliken charges are available for futher calculations.
						next LINE;
					}
				# In the Mulliken charge section, extract lines of the correct form.
					if ($MullikenFlag == 1 && $Line =~ /^\s+\d+\s+\p{IsL}+(\s+|\s+-)\d+.\d+$/) {
						my @Temp = split(/\s+/, $Line);
						# Test to see if element corresponds to list of elements to be captured.
							if($Temp[2] =~ /\p{IsL}+/) {
								# Push data to the correct subarray according to which section of the calculation results the script is bring processing.
									if($NormTerm == 0) {
										push @MullikenAtomNum, $Temp[1];	# This holds the structural atom number.
										push @MullikenElement, $Temp[2];	# This holds the atom element.
										push @MullikenNeutral, $Temp[3];	# This holds the NEUTRAL SPECIES charge.
									}
									if($NormTerm == 1) {push @MullikenAnion, $Temp[3];}		# This holds the ANIONIC SPECIES charge.
									if($NormTerm == 2) {push @MullikenCation, $Temp[3];}	# This holds the CATIONIC SPECIES charge.
							}
						next LINE;
					}
				# At the end of the Mulliken charge section, reset the flag.
					if ($MullikenFlag == 1 && $Line =~ /^\sSum of Mulliken atomic charges/) {
						$MullikenFlag = 0;
						next LINE;
					}
			}

			# Find and extract the Hirshfeld charges.
			{
				# Find the start of the Hirshfeld charges section.
					if ($Line =~ /^\sHirshfeld populations/ && $Stationary == 1) {
						$HirshfeldFlag = 1; # In the Hirshfeld charges section of the output.
						$HirshfeldCharges = 1; # Flag that Hirshfeld charges are available for futher calculations.
						next LINE;
					}
				# In the Hirshfeld charge section, extract lines of the correct form.
					if ($HirshfeldFlag == 1 && $Line =~ /^\s+\d+\s+\p{IsL}+(\s+|\s+-)\d+.\d+/) {
							my @Temp = split(/\s+/, $Line);
							if($Temp[2] =~ /\p{IsL}+/) {
								# Determine the correct Hirschfeld charge from CA and CB (count alpha and beta electrons) columns.
									$Electrons = $Temp[4] + $Temp[5]; # Total electrons associated with atom.
									$Charge = AtomicNumber($Temp[2]) - $Electrons; # Determine the Charge.
								# Push data to the correct subarray according to which section of the calculation results the script is processing.
									if($NormTerm == 0) {
										push @HirshfeldAtomNum, $Temp[1];	# This holds the structural atom number.
										push @HirshfeldElement, $Temp[2];	# This holds the atom element.
										push @HirshfeldNeutral, $Charge;	# This holds the Neutral species charge.
									}
									if($NormTerm == 1) {push @HirshfeldAnion, $Charge;}		# This holds the Anionic species charge.
									if($NormTerm == 2) {push @HirshfeldCation, $Charge;}	# This holds the Cationic species charge.
							}
						next LINE;
					}
				# At the end of the Hirshfeld charge section, reset the flag.
					if ($HirshfeldFlag == 1 && $Line =~ /^\sHirshfeld spin densities/) {
						$HirshfeldFlag = 0;
						next LINE;
					}
			}

			# Find and extract the ESP charges.
			{
				# Find the start of the ESP charge section.
					if ($Line =~ /^\s+Electrostatic Properties Using The SCF Density/ && $Stationary == 1) {
						$ESPFlag = 1; # In the ESP charges section of the output.
						$ESPCharges = 1; # Flag that ESP charges are available for futher calculations.
						next LINE;
					}
				# In the ESP charge section, extract lines of the correct form.
					if ($ESPFlag == 1 && $Line =~ /^\s+\d+\s+\p{IsL}+(\s+|\s+-)\d+.\d+$/) {
						my @Temp = split(/\s+/, $Line);
						# Test to see if element corresponds to list of elements to be captured.
							if($Temp[2] =~ /\p{IsL}+/) {
								# Push data to the correct subarray according to which section of the calculation results the script is bring processing.
									if($NormTerm == 0) {
										push @ESPAtomNum, $Temp[1];	# This holds the structural atom number.
										push @ESPElement, $Temp[2];	# This holds the atom element.
										push @ESPNeutral, $Temp[3];	# This holds the NEUTRAL SPECIES charge.
									}
									if($NormTerm == 1) {push @ESPAnion, $Temp[3];}		# This holds the ANIONIC SPECIES charge.
									if($NormTerm == 2) {push @ESPCation, $Temp[3];}	# This holds the CATIONIC SPECIES charge.
							}
						next LINE;
					}
				# At the end of the ESP charge section, reset the flag.
					if ($ESPFlag == 1 && $Line =~ /^\sCharges from ESP fit with hydrogens summed into heavy atoms/) {
						$ESPFlag = 0;
						next LINE;
					}
			}


			# Find and extract the Natural charges.
			{
				# Find the start of the Natural charge section.
					if ($Line =~ /^\sSummary of Natural Population Analysis/ && $NaturalFlag == 0 && $Stationary == 1) {
						print "FOUND NATURAL CHARGES SECTION.\n";
						$NaturalFlag = 1;
						$NaturalCharges = 1;
						next LINE;
					}
				# In the Natural charge section, extract lines of the correct form.
					if ($NaturalFlag == 1 && $Line =~ /^\s+\p{IsL}+\s+\d+(\s+|\s+-)\d+.\d+/) {
							print "NATURAL LINE: $Line\n";
							@Temp = split(/\s+/, $Line);
							# Do something with the data collectted for this line.  Data structure?
							push @NaturalValue, $Temp[3];		# This holds the atom Natural Charge.
							next LINE;
					}
				# At the end of the Natural charge section (Combined Alpha and Beta), reset the flag.
					if ($NaturalFlag == 1 && $Line =~ /^\s+\* Total \*/) {
						$NaturalFlag = 2;
						print "NATURAL SECTION END\n";
						next LINE;
					}
				# Ensure that Alpha and Beta spin charges are ignored within this Link.
					if ($Line =~ /^\s+NBO analysis skipped by request/i) {
					print "END OF NBO.\n";
					$NaturalFlag = 0;
					next LINE;
					}
			}

	}	# END WHILE loop to read each line of the current file being processed.
	close INPUT; # Close the INPUT File.
	# Issue ERROR message if there are an incorrect number of Normal terminations (G09 claculation didn't succesfully complete).
		if($NormTerm != 3) {
			print "\tERROR: G09 calculation isn't complete. Aborting this file.\n";
			$ErrorCnt++;
			$TerminationErrorCnt++;
			goto CLEANUP; # Abort calculations. Jump to the clean up operations and then process next file.
		}
	# Issue spin contamination warning.
		print "\tWARNING: Spin Contamination Detected." if($SpinContaminationCnt > 1);
	
	# Open File to recieve the cDFT results.
		$FileName = $File; # The new filename is based upon the original filename.
		$FileName =~ s/.out$//i; # Crude strip off the extension.
		$ResultsFile = "cDFT_Results_" . $FileName .".txt";
		open(OUTPUT, ">$ResultsFile") or die "Can't create the results file: $!";
		print "\tSaving to $ResultsFile\n"; # Inform the user of the filename of the results.
	# Write data to output file.
		print OUTPUT "Conceptual DFT Results\n\n";
		print OUTPUT "File:\t$File\n\n";
		# Print energies and Global Indices:
			print OUTPUT "E Neutral (au)\tE Anion (au)\tE Cation (au)\tIP (eV)\tEA (eV)\tSoftness (eV)\tElectrophilicity (eV)\n";
			print OUTPUT "$Energies[0]\t$Energies[1]\t$Energies[2]\t";
			foreach my $x (@Energies) {$x = $x * $au2eV} # Convert energies from au to eV.
			$IP = $Energies[2] - $Energies[0];
			$EA = $Energies[1] - $Energies[0];
			$Hardness = ($IP - $EA) / 2;
			$Softness = $Hardness**-1;
			$Electrophilicity = (($Energies[2] - $Energies[0])+($Energies[0] - $Energies[1])) / (4 * ($Energies[2] + $Energies[1] - 2*$Energies[0]));
			print OUTPUT "$IP\t$EA\t$Softness\t$Electrophilicity\n\n";

		# If Mulliken Charges available then print Fukui Table.
			if($MullikenCharges == 1) {
				print OUTPUT "Mulliken Based Fukui and Related Parameters\n\n";
				print OUTPUT "Atom\tElement\tq Neutral\t q Anion\tq Cation\tf+\tf-\tf·\tS+\tS-\tS·\tw+\tw-\tw·\n"; # Table header.
				# Remove first, blank, element from the Mulliken arrays.
					shift @MullikenAtomNum;
					shift @MullikenElement;
					shift @MullikenNeutral;
					shift @MullikenAnion;
					shift @MullikenCation;
				# Determine the number of Atoms to be processed.
					$Length = scalar(@MullikenAtomNum);
				# Loop through the atoms, calculating and printing the Fukui related results.
					for (my $Index = 0; $Index < $Length; $Index++) {
						# Labels and Charges.
							print OUTPUT "$MullikenAtomNum[$Index]\t$MullikenElement[$Index]\t$MullikenNeutral[$Index]\t$MullikenAnion[$Index]\t$MullikenCation[$Index]\t";
						# Fukui Functions.
							my $FukuiPlus = $MullikenAnion[$Index] - $MullikenNeutral[$Index];
							my $FukuiMinus = $MullikenNeutral[$Index] - $MullikenCation[$Index];
							my $FukuiRad = $MullikenAnion[$Index] - $MullikenCation[$Index];
							print OUTPUT "$FukuiPlus\t$FukuiMinus\t$FukuiRad\t";
						# Softness Functions.
							my $SoftnessPlus = $FukuiPlus * $Softness;
							my $SoftnessMinus = $FukuiMinus * $Softness;
							my $SoftnessRad = $FukuiRad * $Softness;
							print OUTPUT "$SoftnessPlus\t$SoftnessMinus\t$SoftnessRad\t";
						# Electrophilicity Functions.
							my $ElectroPlus = $FukuiPlus * $Electrophilicity;
							my $ElectroMinus = $FukuiMinus * $Electrophilicity;
							my $ElectroRad = $FukuiRad * $Electrophilicity;
							print OUTPUT "$ElectroPlus\t$ElectroMinus\t$ElectroRad\n";
					}
					print OUTPUT "\n";
			}
		
		# If Hirschfeld Charge avaiable then print Fukui Table.
			if($HirshfeldCharges == 1) {
				print OUTPUT "Hirshfeld Based Fukui and Related Parameters\n\n";
				print OUTPUT "Atom\tElement\tq Neutral\t q Anion\tq Cation\tf+\tf-\tf·\tS+\tS-\tS·\tw+\tw-\tw·\n"; # Table header.
				# Remove first, blank, element from the Mulliken arrays.
					shift @HirshfeldAtomNum;
					shift @HirshfeldElement;
					shift @HirshfeldNeutral;
					shift @HirshfeldAnion;
					shift @HirshfeldCation;
				# Determine the number of Atoms to be processed.
					$Length = scalar(@HirshfeldAtomNum);
				# Loop through the atoms, calculating and printing the Fukui related results.
					for (my $Index = 0; $Index < $Length; $Index++) {
						# Labels and Charges.
							print OUTPUT "$HirshfeldAtomNum[$Index]\t$HirshfeldElement[$Index]\t$HirshfeldNeutral[$Index]\t$HirshfeldAnion[$Index]\t$HirshfeldCation[$Index]\t";
						# Fukui Functions.
							my $FukuiPlus = $HirshfeldAnion[$Index] - $HirshfeldNeutral[$Index];
							my $FukuiMinus = $HirshfeldNeutral[$Index] - $HirshfeldCation[$Index];
							my $FukuiRad = $HirshfeldAnion[$Index] - $HirshfeldCation[$Index];
							print OUTPUT "$FukuiPlus\t$FukuiMinus\t$FukuiRad\t";
						# Softness Functions.
							my $SoftnessPlus = $FukuiPlus * $Softness;
							my $SoftnessMinus = $FukuiMinus * $Softness;
							my $SoftnessRad = $FukuiRad * $Softness;
							print OUTPUT "$SoftnessPlus\t$SoftnessMinus\t$SoftnessRad\t";
						# Electrophilicity Functions.
							my $ElectroPlus = $FukuiPlus * $Electrophilicity;
							my $ElectroMinus = $FukuiMinus * $Electrophilicity;
							my $ElectroRad = $FukuiRad * $Electrophilicity;
							print OUTPUT "$ElectroPlus\t$ElectroMinus\t$ElectroRad\n";
					}
				print OUTPUT "\n";
			}

		# If ESP Charge avaiable then print Fukui Table.
			if($ESPCharges == 1) {
				print OUTPUT "Electrostatic Potential (MK) Based Fukui and Related Parameters\n\n";
				print OUTPUT "Atom\tElement\tq Neutral\t q Anion\tq Cation\tf+\tf-\tf·\tS+\tS-\tS·\tw+\tw-\tw·\n"; # Table header.
				# Remove first, blank, element from the ESP arrays.
					shift @ESPAtomNum;
					shift @ESPElement;
					shift @ESPNeutral;
					shift @ESPAnion;
					shift @ESPCation;
				# Determine the number of Atoms to be processed.
					$Length = scalar(@ESPAtomNum);
				# Loop through the atoms, calculating and printing the Fukui related results.
					for (my $Index = 0; $Index < $Length; $Index++) {
						# Labels and Charges.
							print OUTPUT "$ESPAtomNum[$Index]\t$ESPElement[$Index]\t$ESPNeutral[$Index]\t$ESPAnion[$Index]\t$ESPCation[$Index]\t";
						# Fukui Functions.
							my $FukuiPlus = $ESPAnion[$Index] - $ESPNeutral[$Index];
							my $FukuiMinus = $ESPNeutral[$Index] - $ESPCation[$Index];
							my $FukuiRad = $ESPAnion[$Index] - $ESPCation[$Index];
							print OUTPUT "$FukuiPlus\t$FukuiMinus\t$FukuiRad\t";
						# Softness Functions.
							my $SoftnessPlus = $FukuiPlus * $Softness;
							my $SoftnessMinus = $FukuiMinus * $Softness;
							my $SoftnessRad = $FukuiRad * $Softness;
							print OUTPUT "$SoftnessPlus\t$SoftnessMinus\t$SoftnessRad\t";
						# Electrophilicity Functions.
							my $ElectroPlus = $FukuiPlus * $Electrophilicity;
							my $ElectroMinus = $FukuiMinus * $Electrophilicity;
							my $ElectroRad = $FukuiRad * $Electrophilicity;
							print OUTPUT "$ElectroPlus\t$ElectroMinus\t$ElectroRad\n";
					}
				print OUTPUT "\n";
			}

	# Based upon which charges have been found in the output file.
		# Calculate the Condensed to Atom Fukui functions (f-, f+, f°) for each charge system.
		# Calculate the corresponding Softness Indices (S-, S+, S°) for each charge system.
		# Calculate the electrophilicity Indices (w+, w-, w°) for each charge system.

	# Print results to an output table (one results file per input file). If error free.
		# (cDFTResults_FILENAME.txt)
		# OUTPUT:
			# Filename.
			# Energies (au)
			# IP, EA, Hardness, Softness, Electophilicity index (eV)
			# Mulliken Charge and Index table if data available.
			# NPA Charge and Index table if data available.
			# Hirshfeld Charge and Index table if data available.
			# MK Charge and Index table if data available.
			# Issue any errors/warnings about the data.

	# Issue warnings wrt Spin Contamination:
		print "\tSpin $SpinContaminationCnt\n";
		print "\tSpin Contamination Detected. Treat results with caution.\n" if($SpinContaminationCnt == 1);
		print "\tTwo Spin Contaminations Detected. Treat results with caution.\n" if($SpinContaminationCnt == 2);

		close OUTPUT;


	CLEANUP:	# Reset all variables in preparation for the next file.
		@Energies = (0,0,0);	# Reset the energies.
		$NormTerm = 0;	# Reset the Normal terminations counter.
		$SpinContaminationCnt = 0; # Reset the Spin Contamination counter.
		$Stationary = 0; # Reset the Stationary Point Found Flag.




	} # END FOREACH loop to process each file in turn.

# Any tidy up required.

	GoodBye();
	exit; # Catch All exit.

# MAIN CODE ENDS HERE #


# SUBROUTINES:

sub Hello
{
print <<ENDHELLO;

		                       ( ~!~ ~!~ )
		 .-----------------oOOo----(_)----oOOo---------------------------.
		 |                                                               |
		 |                     A Perl Script for:                        |
		 |                                                               |
		 |  ..--** The extraction of key parameters from G09     **--..  |
		 |  ..--** output files for Conceptual DFT calculations. **--..  |
		 |                                                               |
		 | Version 0.1 - First Release                                   |
		 | Author David Hose                                             |
		 | (c) July 2014                                                 |
		 |                      .oooO                                    |
		 |                      (   )   Oooo.                            |
		 .-----------------------\\ (----(   )---------------------------.
		                          \\_)    ) /
		                                (_/

ENDHELLO
}

sub GoodBye
{
print <<ENDGOODBYE;

Script completed. Goodbye.

ENDGOODBYE
}


# Parser


# Parser errors

# Usage


sub AtomicNumber {
# This function returns the atomic number of the element passed to it.
	# Define atomic numbers in a hash.
		my %AtomicZ =	('H', '1', 'He', '2', 'Li', '3', 'Be', '4', 'B', '5', 
						'C', '6', 'N', '7', 'O', '8', 'F', '9', 'Ne', '10', 
						'Na', '11', 'Mg', '12', 'Al', '13', 'Si', '14', 'P', '15', 
						'S', '16', 'Cl', '17', 'Ar', '18', 'K', '19', 'Ca', '20', 
						'Sc', '21', 'Ti', '22', 'V', '23', 'Cr', '24', 'Mn', '25', 
						'Fe', '26', 'Co', '27', 'Ni', '28', 'Cu', '29', 'Zn', '30', 
						'Ga', '31', 'Ge', '32', 'As', '33', 'Se', '34', 'Br', '35', 
						'Kr', '36', 'Rb', '37', 'Sr', '38', 'Y', '39', 'Zr', '40', 
						'Nb', '41', 'Mo', '42', 'Tc', '43', 'Ru', '44', 'Rh', '45', 
						'Pd', '46', 'Ag', '47', 'Cd', '48', 'In', '49', 'Sn', '50', 
						'Sb', '51', 'Te', '52', 'I', '53', 'Xe', '54', 'Cs', '55', 
						'Ba', '56', 'La', '57', 'Ce', '58', 'Pr', '59', 'Nd', '60', 
						'Pm', '61', 'Sm', '62', 'Eu', '63', 'Gd', '64', 'Tb', '65', 
						'Dy', '66', 'Ho', '67', 'Er', '68', 'Tm', '69', 'Yb', '70', 
						'Lu', '71', 'Hf', '72', 'Ta', '73', 'W', '74', 'Re', '75', 
						'Os', '76', 'Ir', '77', 'Pt', '78', 'Au', '79', 'Hg', '80', 
						'Tl', '81', 'Pb', '82', 'Bi', '83', 'Po', '84', 'At', '85', 
						'Rn', '86', 'Fr', '87', 'Ra', '88', 'Ac', '89', 'Th', '90', 
						'Pa', '91', 'U', '92', 'Np', '93', 'Pu', '94', 'Am', '95', 
						'Cm', '96', 'Bk', '97', 'Cf', '98', 'Es', '99', 'Fm', '100', 
						'Md', '101', 'No', '102', 'Lr', '103', 'Rf', '104', 'Db', '105', 
						'Sg', '106', 'Bh', '107', 'Hs', '108', 'Mt', '109', 'Ds', '110', 
						'Rg', '111', 'Cn', '112', 'Uut', '113', 'Fl', '114', 'Uup', '115', 
						'Lv', '116', 'Uus', '117', 'Uuo', '118');
	my $Element = @_[0];
	$Element =~ s/([A-Z][a-z]?)/\$AtomicZ{$1}/g;
	my $Charge = eval($Element);
	return($Charge);
}


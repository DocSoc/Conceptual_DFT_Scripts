#! /usr/bin/perl




# This script extracts the cDFT data.
# Sister script to pdb2cDFT.pl script.


# Say "Hi" to the user.


# parse STDIN for:
	# Filenames (GLOB)
	# Switches that are set.
		# Atoms.
			# "all": All atoms.
			# "heavies": All atoms except hydrogen.
			# "Default": Reactive centre atoms; C, N, S, Si, P,... (exclude halogens, oxygen,...)

# START FOREACH loop to process each file in turn.


	# START WHILE loop to read each line of the current file being processed.
	
		# Find "SCF Done" energies.
			# Store the last three energies in @Energies array.  If there are no errors, this will contain the neutral, radical anion and radical anion energies.
		# Count the number of "Normal termination" in the file.  There should be 3. If less than this, issue an error message.
		# Check the spin contamination for the radicals.  If outside of a specified range, issue a warning.
		# Keep track of which section (Neutral, Radical Anion or Radial Cation) of the calculation outout the script in operating in.
		# Extract charges (only extract charges for the required elements):
			# Extract Mulliken charges:
				# For Neutral species.
				# For Radical Anion species.
				# For Radical Cation species.
				# Set Flag that Mulliken Charges are available for calculation.
			# Extract NPA charges:
				# For Neutral species.
				# For Radical Anion species.
				# For Radical Cation species.
				# Set Flag that NPA Charges are available for calculation.
			# Extract Hirshfeld charges:
				# For Neutral species.
				# For Radical Anion species.
				# For Radical Cation species.
				# Set Flag that Hirshfeld Charges are available for calculation.
			# Extract MK (ESP) charges:
				# For Neutral species.
				# For Radical Anion species.
				# For Radical Cation species.
				# Set Flag that MK Charges are available for calculation.

	# END WHILE loop to read each line of the current file being processed.

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

	# Reset all variables in preparation for the next file.






# END FOREACH loop to process each file in turn.

# Any tidy up required.

# Say "Bye" to the user.



# SUBROUTINES:

# Parser


# Parser errors

# Usage


# 


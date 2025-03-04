function initialize_drawing(tag )
% initialize_drawing( tag )
%
%  Master function for starting RiboDraw, based on output of
%  Rosetta rna_motif run on a PDB file.
%
% INPUT: 
%  tag = name of input PDB file, possibly including path. 
%            If tag is "RNA.pdb", this function expects filenames
%            with the names
%
%                RNA.pdb.fasta
%                RNA.pdb.base_pairs.txt
%                RNA.pdb.stacks.txt
%                RNA.pdb.other_contacts.txt
%                RNA.pdb.stems.txt
%                RNA.pdb.ligands.txt    
%                RNA.pdb.motifs.txt    
%
% (C) R. Das, Stanford University, 2017
% (C) Hao Sun, GuangZhou National Labortory, 2025

%% 1. Read Input Files
% -------------------------
% Retrieve the sequence, residue numbers, chain identifiers, segment IDs, 
% and non-standard residue information from the .fasta file.
disp("1. Starting to read input files");
[sequence, resnum, chains, segid, non_standard_residues] = get_sequence([tag, '.fasta']);

% Read base pairs information, including noncanonical pairs.
base_pairs = read_base_pairs([tag, '.base_pairs.txt']);

% Check whether the indices in base_pairs match the indices from the FASTA file.
check_index(base_pairs, resnum, chains)

% Read base stacks information, including noncanonical stacks.
base_stacks = read_base_stacks([tag, '.stacks.txt']);

% Read other contact information (may include non-standard intermolecular contacts).
other_contacts = read_other_contacts([tag, '.other_contacts.txt']);

% Read the helix information (stems) for the RNA secondary structure.
stems = read_stems([tag, '.stems.txt']);

% Read ligands information, such as small molecules or ions.
ligands = read_ligands([tag, '.ligands.txt']);

% Read motifs information.
motifs = read_motifs([tag, '.motifs.txt']);

% If base pair information is empty but stems exist, generate base pairs from stems.
if length(base_pairs) == 0 & length(stems) > 0
    disp("Base pairs are empty but stems exist. Generating base pairs from stems.");
    base_pairs = get_base_pairs_from_stems(stems);
end
disp("1. Finished reading input files");

%% 2. Set Up Drawing Environment
% -------------------------
% Clear the current figure and set the current axes (gca) to occupy the full window.
disp("2. Setting up the drawing environment");
clf; 
set(gca, 'Position', [0 0 1 1]);

% Hold the current plot to allow multiple drawings.
hold on;

% Initialize a zero vector the same length as the sequence (possibly used for later calculations or markings).
t = zeros(1, length(sequence));

% Set the axis limits.
axis([0 200 0 200]);

% Save the default plotting settings to the current axes' application data for later use.
setappdata(gca, 'plot_settings', default_plot_settings);
disp("2. Drawing environment setup complete");

%% 3. Initialize Structural Objects and Data
% -------------------------
% Set the default positions for the stem (helix) objects, e.g., determining each helix's center.
disp("3. Initializing structural objects and data: Processing helix objects");
disp("3. Setting default positions for helix objects");
stems = set_default_stem_positions(stems);

% Initialize each residue object based on the sequence, residue numbers, chain identifiers, and segment IDs,
% and assign the corresponding helix_tag to each residue.
disp("3. Initializing each residue object based on sequence, residue numbers, chain IDs, and segment IDs, and assigning helix_tag");
stems = setup_residues(stems, sequence, resnum, chains, segid);

% Set up the pairing between stems (e.g., preparing for drawing connecting lines later).
disp("3. Setting up pairing between helix objects");
setup_stem_partner(stems);

%% 4. Initialize Linkers
% -------------------------
% Set up arrow linkers, typically used to indicate directionality or connectivity, 
% based on residue numbers, chain identifiers, and segment IDs.
disp("4. Initializing linkers");
disp("4. Setting up arrow linkers based on residue numbers, chain IDs, and segment IDs");
setup_arrow_linkers(resnum, chains, segid);

% Set up linkers for base pairs to display the connections between base pairs.
disp("4. Setting up base pair linkers to display connections between base pairs");
setup_base_pair_linkers(base_pairs);

% Set up linkers for base stacks to display the interactions in base stacking.
disp("4. Setting up base stack linkers to display interactions in base stacking");
setup_base_stack_linkers(base_stacks);

% Set up linkers for other contacts to display noncanonical contacts.
disp("4. Setting up linkers for other contacts to display noncanonical contacts");
setup_other_contact_linkers(other_contacts);

% Attempt to process non-standard residue names, possibly normalizing or converting names.
disp("4. Attempting to process non-standard residue names (normalization/conversion)");
try_non_standard_names(sequence, resnum, chains, segid, non_standard_residues);

%% 5. Initialize Other Structures (Motifs, Ligands, etc.)
% -------------------------
% The following two lines (regarding coaxial stacking) are commented out; 
% enable if needed:
% coaxial_stacks = get_coaxial_stacks(base_pairs, base_stacks, stems);
% setup_coaxial_stacks(coaxial_stacks);
disp("5. Initializing other structures (motifs, ligands, etc.)");
% Set up ligand information to prepare for drawing ligand molecules.
disp("5. Setting up ligand information for drawing ligand molecules");
setup_ligands(ligands);

% Set up motifs information to prepare for drawing structural motifs.
disp("5. Setting up motifs information for drawing structural motifs");
setup_motifs(motifs);

%% 6. Draw Structures and Final Settings
% -------------------------
% Draw the RNA helix structure based on the stem information.
disp("6. Drawing the RNA helix structure based on stem information (using stems)");
draw_helices(stems);

% Initialize zoom functionality to allow users to zoom in and out of the drawing.
disp("6. Initializing zoom functionality to allow user zooming");
setup_zoom();

% Set the properties of the axes to enhance the overall appearance.
disp("6. Setting axes properties for improved appearance");
set_nice_axes();

% Update the graphics size to adapt to the current window or screen settings.
disp("6. Updating graphics size to adapt to the current window/screen");
update_graphics_size();
update_graphics_size();

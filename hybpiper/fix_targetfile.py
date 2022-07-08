#!/usr/bin/env python

"""
########################################################################################################################

Takes a nucleotide or protein target file in fasta format as input.

If a nucleotide target file is provided:

    - For each sequence, identifies all forwards frames that can be translated without stop codons. If only one
      such frame exists, this sequence is writen to a new target file with any trailing 3' nucleotides removed (i.e.
      it writes a sequence that is a multiple of three, containing full codons only).

    - If more than one possible translation frame exists, an attempt to identify the correct frame is made by:
        * comparing each possible translation to a reference protein (if provided), or
        * comparing each possible translation to the longest protein translation from all other representative
          sequences for the given gene (if such sequences are present in the target file and have only one possible
          translation frame)

    - If all forwards frames contain at least one unexpected stop codon, these sequences will be written to a second
      fasta file.

    - If a sequence has more than one possible translation frame and the correct frame can't be determined,
      sequences for all possible translation frames will be written to a third fasta file.

    - If a sequence has more than one possible translation frame and all frames are greater than a maximum distance
      from the reference protein (if present), sequences for all possible translation frames will be written to a
      fourth fasta file. The sensitivity of this comparison can be adjusted with the parameter "--maximum_distance",
      and is useful to filter out sequences with frameshifts that do NOT introduce stop codons. 0.0 means identical
      sequences, 1.0 means completely different sequences. Default is 0.5.

    - OPTIONAL filter by length.

    - OPTIONAL remove seqs with low complexity regions.


If a protein target file is provided:

    - OPTIONAL filter by length.

    - OPTIONAL remove seqs with low complexity regions.

########################################################################################################################
"""

import logging
import sys
import argparse
import os
import socket
import re
from collections import defaultdict
import datetime
import copy
from concurrent.futures.process import ProcessPoolExecutor
from multiprocessing import Manager
from concurrent.futures import wait
import glob
import textwrap
import io
import shutil

# Import non-standard-library modules:
unsuccessful_imports = []
try:
    from Bio import SeqIO, SeqRecord, AlignIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align.Applications import MafftCommandline
    from Bio.Phylo.TreeConstruction import DistanceCalculator
except ImportError:
    unsuccessful_imports.append('Bio')
try:
    import progressbar
except ImportError:
    unsuccessful_imports.append('progressbar')

if unsuccessful_imports:
    package_list = '\n'.join(unsuccessful_imports)
    sys.exit(f'The required Python packages are not found:\n\n{package_list}\n\nAre they installed for the Python '
             f'installation used to run this script?')

# f-strings will produce a 'SyntaxError: invalid syntax' error if not supported by Python version:
f'Must be using Python 3.6 or higher.'

if sys.version_info[0:2] < (3, 6):
    sys.exit(f'Must be using Python 3.6 or higher. You are using version {sys.version_info[0]}.{sys.version_info[1]}.')

########################################################################################################################
########################################################################################################################
# Get current working directory and host name:

cwd = os.getcwd()
host = socket.gethostname()

# Set widget format for progressbar:
widgets = [' ' * 11,
           progressbar.Timer(),
           progressbar.Bar(),
           progressbar.ETA()]


# Configure logger:
def setup_logger(name, log_file, console_level=logging.INFO, file_level=logging.DEBUG,
                 logger_object_level=logging.DEBUG):
    """
    Function to create a logger instance.

    By default, logs level DEBUG and above to file.
    By default, logs level INFO and above to stderr and file.

    :param string name: name for the logger instance
    :param string log_file: filename for log file
    :param string console_level: level for logging to console
    :param string file_level: level for logging to file
    :param string logger_object_level: level for logger object
    :return: a logger object
    """

    # Get date and time string for log filename:
    date_and_time = datetime.datetime.now().strftime("%Y-%m-%d-%H_%M_%S")

    # Log to file:
    file_handler = logging.FileHandler(f'{log_file}_{date_and_time}.log', mode='w')
    file_handler.setLevel(file_level)
    file_format = logging.Formatter('%(asctime)s - %(filename)s - %(name)s - %(funcName)s - %(levelname)s - %('
                                    'message)s')
    file_handler.setFormatter(file_format)

    # Log to Terminal (stdout):
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(console_level)
    console_format = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_format)

    # Setup logger:
    logger_object = logging.getLogger(name)
    logger_object.setLevel(logger_object_level)  # Default level is 'WARNING'

    # Add handlers to the logger
    logger_object.addHandler(console_handler)
    logger_object.addHandler(file_handler)

    return logger_object


# Create logger(s):
logger = setup_logger(__name__, 'fix_targetfile')

########################################################################################################################
########################################################################################################################
# Define functions:


def createfolder(directory):
    """
    Attempts to create a directory named after the name provided, and provides an error message on failure
    """
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        logger.info(f'Error: Creating directory: {directory}')


def file_exists_and_not_empty(file_name):
    """
    Check if file exists and is not empty by confirming that its size is not 0 bytes
    """

    return os.path.isfile(file_name) and not os.path.getsize(file_name) == 0


def done_callback(future_returned):
    """
    Callback function for ProcessPoolExecutor futures; gets called when a future is cancelled or 'done'.
    """
    if future_returned.cancelled():
        logger.info(f'{future_returned}: cancelled')
        return
    elif future_returned.done():
        error = future_returned.exception()
        if error:
            logger.info(f'{future_returned}: error returned: {error}')
            return future_returned.result()
        else:
            result = future_returned.result()

    return result


def align(fasta_file, algorithm, output_folder, counter, lock, num_files_to_process, threads_per_concurrent_alignment):
    """
    Uses mafft to align a fasta file of sequences, using the algorithm and number of threads provided. Returns filename
    of the alignment produced.

    :param fasta_file:
    :param algorithm:
    :param output_folder:
    :param counter:
    :param lock:
    :param num_files_to_process:
    :param threads_per_concurrent_alignment:
    :return:
    """

    createfolder(output_folder)
    fasta_file_basename = os.path.basename(fasta_file)
    expected_alignment_file = f'{output_folder}/{re.sub(".fasta", ".aln.fasta", fasta_file_basename)}'

    try:
        assert file_exists_and_not_empty(expected_alignment_file)
        logger.debug(f'Alignment exists for {fasta_file_basename}, skipping...')
        with lock:
            counter.value += 1
        return os.path.basename(expected_alignment_file)
    except AssertionError:
        mafft_cline = (MafftCommandline(algorithm, thread=threads_per_concurrent_alignment, input=fasta_file))
        stdout, stderr = mafft_cline()
        with open(expected_alignment_file, 'w') as alignment_file:
            alignment_file.write(stdout)
            with lock:
                counter.value += 1
        return os.path.basename(expected_alignment_file)
    finally:
        sys.stderr.write(f'\r{"[INFO]:":10} Finished generating alignment for file {fasta_file_basename}, '
                         f' {counter.value}/{num_files_to_process}')


def align_extractions_multiprocessing(unaligned_folder, output_folder, concurrent_alignments,
                                      threads_per_concurrent_alignment):
    """
    Generate alignments via function <align> using multiprocessing.

    :param unaligned_folder:
    :param output_folder:
    :param concurrent_alignments:
    :param threads_per_concurrent_alignment:
    :return:
    """

    createfolder(output_folder)
    alignments = [file for file in sorted(glob.glob(f'{unaligned_folder}/*.fasta'))]

    with ProcessPoolExecutor(max_workers=concurrent_alignments) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(align,
                                      alignment,
                                      'linsi',
                                      output_folder,
                                      counter,
                                      lock,
                                      num_files_to_process=len(alignments),
                                      threads_per_concurrent_alignment=threads_per_concurrent_alignment)
                          for alignment in alignments]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")

        alignment_list = [alignments for alignment in glob.glob(f'{output_folder}/*.aln.fasta') if
                          file_exists_and_not_empty(alignment)]

        logger.info(f'\n{"[INFO]:":10} {len(alignment_list)} alignments generated from {len(future_results)} '
                    f'alignment files...\n')

        if len(alignments) != len(alignment_list):
            sys.exit(f'Only {len(alignment_list)} alignments were generated from {len(alignments)} genes, check '
                     f'for errors!')


def pad_seq(sequence):
    """
    Pads a sequence Seq object to a multiple of 3 with 'N'.

    :param Bio.SeqRecord.SeqRecord sequence: sequence to pad
    :return: Bio.SeqRecord.SeqRecord sequence padded with Ns if required, int for number of Ns added to 3' end
    """

    remainder = len(sequence.seq) % 3
    if remainder == 0:
        return sequence, 0
    else:
        sequence.seq = sequence.seq + Seq('N' * (3 - remainder))
        return sequence, (3 - remainder)


def choose_best_match_translation(sequence_list, gene_to_inframe_seqs, ref_is_protein=False, maximum_distance=0.5):
    """
    Compares translations from multiple frames to a reference protein sequence via a distance matrix, and selects the
    translation frame that is most similar.

    :param list sequence_list: list of SeqRecord protein translation seqs to test against reference
    :param Bio.SeqRecord.SeqRecord/list gene_to_inframe_seqs: list of SeqRecord protein translation seqs for good
    gene references
    :param bool ref_is_protein: if True, the reference sequence is amino acids from external fasta file
    :param float maximum_distance: maximum distance allowed for any frame to be selected
    :return:
    """

    seq_to_keep = None
    distance = 1.0  # i.e. all positions different
    temp_fasta_file = f'temp.fasta'
    for seq_to_test in sequence_list:
        seq_to_test_copy = copy.deepcopy(seq_to_test)
        seq_record = SeqRecord(Seq(seq_to_test_copy.seq.translate()),
                               id=f'{seq_to_test.id}',
                               name=f'{seq_to_test.name}',
                               description=f'')
        if ref_is_protein:
            ref_protein = gene_to_inframe_seqs
            temp_seqs_to_write = [ref_protein, seq_record]
        else:
            longest_reference = max(gene_to_inframe_seqs, key=len)
            ref_protein = copy.deepcopy(longest_reference)
            ref_protein.seq = ref_protein.seq.translate()
            logger.debug(f'{"[INFO]:":10} Longest reference is: {ref_protein.name}')
            temp_seqs_to_write = [ref_protein, seq_record]

        with open(temp_fasta_file, 'w') as temp_handle:
            SeqIO.write(temp_seqs_to_write, temp_handle, 'fasta')

        mafft_cline = (MafftCommandline('linsi', thread=1, input=temp_fasta_file))
        stdout, stderr = mafft_cline()
        alignment = AlignIO.read(io.StringIO(stdout), 'fasta')
        for sequence in alignment:  # Distance matrix requires capital letters
            sequence.seq = sequence.seq.upper()

        # Create a distance matrix
        skip_letters = set(letter for sequence in alignment for letter in sequence.seq if letter not in
                           ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W',
                            'Y', 'V'])
        my_calculator = DistanceCalculator('blosum62', skip_letters=''.join(skip_letters))
        trimmed_dm = my_calculator.get_distance(alignment)
        distance_values = trimmed_dm[seq_to_test.name]
        if distance_values[0] < distance:
            distance = distance_values[0]
            seq_to_keep = seq_to_test

    logger.debug(f'{"[INFO]:":10} Keeping translation {seq_to_keep.name} with distance {distance}')

    os.remove(temp_fasta_file)

    if distance > maximum_distance:
        logger.debug(f'{"[INFO]:":10} Distance between selected sequence {seq_to_keep.name} and protein reference'
                     f' {ref_protein.name} is {distance}. Maximum distance is {maximum_distance}. No "best" '
                     f'translation frame found!')
        return None

    return seq_to_keep


def get_inframe_sequence(target_fasta_file, no_terminal_stop_codons, reference_protein_file, maximum_distance,
                         filter_by_length_percentage, allow_gene_removal, keep_low_complexity_sequences,
                         low_complexity_seq_names):
    """
    Recovers sequences with no stop codons and full codon triplets only, if present, and writes them to a fasta file.
    Sequences without a full-length open reading frame are written to a second fasta file. Sequences with more than
    one full-length open reading frame are tested against a reference protein (if provided) or other known protein
    translations of the same gene (if present); if the correct frame can't be determined, such sequences are written
    to a third file.

    :param str target_fasta_file: path to the target fasta file.
    :param bool no_terminal_stop_codons: if True, do not allow a translated amino-acid sequence to have a stop codon at
    the C-terminus.
    :param None or str reference_protein_file: path to a fasta file containing reference proteins
    :param float maximum_distance: maximum distance allowed for any frame to be selected
    :param float filter_by_length_percentage: only include sequences longer than this % of the longest gene sequence
    :param bool allow_gene_removal: if True, allow all seqs for a given gene to be removed
    :param bool keep_low_complexity_sequences: if True, keep sequences that contain regions of low-complexity
    :param list low_complexity_seq_names: list of sequence names for seqs with regions of low-complexity
    :return collections.defaultdict gene_to_inframe_seq_dictionary: a dictionary of gene to SeqRecord objects for all
    fixed sequences
    """

    file, ext = os.path.splitext(os.path.basename(target_fasta_file))
    fixed_file_output = f'{file}_inframe{ext}'
    fixed_file_output_with_stop_codons = f'{file}_stop_codons_all_frames{ext}'
    fixed_file_output_with_undetermined_frame = f'{file}_undetermined_frames{ext}'
    fixed_file_output_with_multiple_frames_above_maximum_distance = f'{file}_above_maximum_distance_frames{ext}'
    fixed_file_output_low_complexity_sequences= f'{file}_low_complexity_regions{ext}'

    # If a file of reference proteins is provided, create a gene_name:protein dictionary:
    if reference_protein_file:
        with open(reference_protein_file) as reference_handle:
            seqs = SeqIO.parse(reference_handle, 'fasta')
            reference_protein_dict = dict()
            for seq in seqs:
                gene_id = seq.name.split('-')[-1]
                reference_protein_dict[gene_id] = seq
        logger.info(f'{"[INFO]:":10} A fasta file of reference protein sequences has been provided, containing a '
                    f'sequence for {len(reference_protein_dict)} genes.')

    # Check that seqs in target_fasta_file can be translated in one of the forwards frames:
    sequences = list(SeqIO.parse(target_fasta_file, 'fasta'))
    original_gene_names = set([sequence.name.split('-')[-1] for sequence in sequences])
    logger.info(f'{"[INFO]:":10} Input fasta file contains at least one sequence for {len(original_gene_names)} genes')

    sequences_with_stop_codons_all_frames = []
    sequences_with_multiple_possible_frames_above_similarity_threshold = []
    sequences_with_multiple_possible_frames_above_similarity_threshold_count = 0
    gene_to_multiframe_seq_dictionary = defaultdict(lambda: defaultdict(list))  # nucleotide when no external prot ref
    gene_to_inframe_seq_dictionary = defaultdict(list)  # nucleotide dict for prot alignments

    logger.info(f'{"[INFO]:":10} Identifying sequences with a single forwards open reading frame...')
    for sequence in progressbar.progressbar(sequences,
                                            max_value=len(sequences),
                                            min_poll_interval=10,
                                            widgets=widgets):

        gene_id = sequence.name.split('-')[-1]
        taxon_id = '-'.join(sequence.name.split('-')[:-1])
        frames_without_stop_codons_seqs = []

        for frame_start in [0, 1, 2]:
            sequence_to_test = copy.deepcopy(sequence)
            sequence_to_test.seq = sequence_to_test.seq[frame_start:]
            padded_sequence, ns_added = pad_seq(sequence_to_test)
            padded_sequence_translated = padded_sequence.seq.translate()
            num_stop_codons = padded_sequence_translated.count('*')

            if num_stop_codons == 0 or \
                    (num_stop_codons == 1 and
                     re.search('[*]', str(padded_sequence_translated)[-1]) and not
                     no_terminal_stop_codons):

                logger.debug(f'{"[INFO]:":10} Translated sequence {padded_sequence.name} does not contain any '
                             f'unexpected stop codons in frame {frame_start + 1}')

                if ns_added == 0:
                    seq_record = SeqRecord(Seq(padded_sequence.seq),
                                           id=f'{sequence.id}_frame_{frame_start + 1}',
                                           name=f'{sequence.name}_frame_{frame_start + 1}',
                                           description=f'')
                    frames_without_stop_codons_seqs.append(seq_record)
                else:
                    seq_record = SeqRecord(Seq(padded_sequence.seq[:-3]),
                                           id=f'{sequence.id}_frame_{frame_start + 1}',
                                           name=f'{sequence.name}_frame_{frame_start + 1}',
                                           description=f'')
                    frames_without_stop_codons_seqs.append(seq_record)

            elif num_stop_codons:

                logger.debug(f'{"[INFO]:":10} Translated sequence {padded_sequence.name} contains at least one '
                             f'unexpected stop codon in frame {frame_start + 1}, moving on...')

                if len(frames_without_stop_codons_seqs) == 0 and frame_start in [2]:
                    sequences_with_stop_codons_all_frames.append(sequence)

        # Check if more than one forward frame could be translated without stop codons:
        if len(frames_without_stop_codons_seqs) == 1:
            seq_correct_frame = frames_without_stop_codons_seqs[0]
            seq_correct_frame.name = '_'.join(seq_correct_frame.name.split('_')[:-2])
            seq_correct_frame.id = '_'.join(seq_correct_frame.id.split('_')[:-2])
            gene_to_inframe_seq_dictionary[gene_id].append(seq_correct_frame)
        elif len(frames_without_stop_codons_seqs) > 1:
            if reference_protein_file:
                try:
                    ref_protein = reference_protein_dict[gene_id]
                    logger.debug(f'{"[INFO]:":10} Found reference protein sequence for gene {gene_id} sequence '
                                 f'{sequence.id}')
                    best_match_seq = choose_best_match_translation(frames_without_stop_codons_seqs,
                                                                   ref_protein,
                                                                   ref_is_protein=True,
                                                                   maximum_distance=maximum_distance)
                    if best_match_seq:
                        best_match_seq.name = '_'.join(best_match_seq.name.split('_')[:-2])
                        best_match_seq.id = '_'.join(best_match_seq.id.split('_')[:-2])
                        gene_to_inframe_seq_dictionary[gene_id].append(best_match_seq)
                    else:
                        logger.debug(f'{"[INFO]:":10} Multiple possible translation frames are present for {gene_id} '
                                     f'sequence {sequence.id}, but none are less than the maximum allowed distance of'
                                     f' {maximum_distance}.')
                        sequences_with_multiple_possible_frames_above_similarity_threshold.extend(
                            frames_without_stop_codons_seqs)
                        sequences_with_multiple_possible_frames_above_similarity_threshold_count += 1

                except KeyError:
                    logger.debug(f'{"[INFO]:":10} Sequence {sequence.name} can be translated in more than one '
                                 f'forwards frame, and no reference protein gene was provided.')
                    gene_to_multiframe_seq_dictionary[gene_id][taxon_id].extend(frames_without_stop_codons_seqs)
            else:
                gene_to_multiframe_seq_dictionary[gene_id][taxon_id].extend(frames_without_stop_codons_seqs)

    # For sequences that can be translated in more than one forwards frame without stop codons, check if there are good
    # reference sequences for that gene, compare frame translations against the longest good reference, and select
    # the most similar frame translation:
    seqs_with_undetermined_frame = []
    seqs_with_undetermined_frame_count = 0

    logger.info(f'{"[INFO]:":10} Attempting to identify the correct reading frame for sequences with more than one '
                f'possible forward translation...')
    for gene_id, taxon_dict in progressbar.progressbar(gene_to_multiframe_seq_dictionary.items(),
                                                       max_value=len(gene_to_multiframe_seq_dictionary),
                                                       min_poll_interval=10,
                                                       widgets=widgets):
        for taxon, sequence_list in taxon_dict.items():
            if len(gene_to_inframe_seq_dictionary[gene_id]) > 0:  # i.e. at least one good ref sequence for this gene
                best_match_seq = choose_best_match_translation(sequence_list,
                                                               gene_to_inframe_seq_dictionary[gene_id],
                                                               ref_is_protein=False,
                                                               maximum_distance=maximum_distance)
                if best_match_seq:
                    best_match_seq.name = '_'.join(best_match_seq.name.split('_')[:-2])
                    best_match_seq.id = '_'.join(best_match_seq.id.split('_')[:-2])
                    gene_to_inframe_seq_dictionary[gene_id].append(best_match_seq)
                else:
                    logger.debug(f'{"[INFO]:":10} Multiple possible translation frames are present for {gene_id} '
                                 f'sequence {sequence_list[0].id.split("_")[:-2]}, but none are less than the maximum '
                                 f'distance of {maximum_distance}.')
                    sequences_with_multiple_possible_frames_above_similarity_threshold.extend(
                        sequence_list)
                    sequences_with_multiple_possible_frames_above_similarity_threshold_count += 1
            else:
                logger.debug(f'{"[INFO]:":10} The correct translation frame can not be determined for sequence '
                             f'{taxon}-{gene_id}. Sequences corresponding to possible open reading frames will be '
                             f'written to file.')
                seqs_with_undetermined_frame.extend(sequence_list)
                seqs_with_undetermined_frame_count += 1

    # Log some info messages:
    logger.info(f'{"[INFO]:":10} Number of sequences with at least one stop codon in all forwards frames:'
                f' {len(sequences_with_stop_codons_all_frames)}')

    logger.info(f'{"[INFO]:":10} Number of sequences with multiple possible reading frames but no reference; correct '
                f'frame undetermined: {seqs_with_undetermined_frame_count}')

    logger.info(f'{"[INFO]:":10} Number of sequences with all possible reading frames are greater than minimum '
                f'distance from reference: {sequences_with_multiple_possible_frames_above_similarity_threshold_count}')

    inframe_seqs_to_write = [seq for gene, sequences in gene_to_inframe_seq_dictionary.items() for seq in sequences]

    assert (len(inframe_seqs_to_write) + len(sequences_with_stop_codons_all_frames)) + seqs_with_undetermined_frame_count + \
           sequences_with_multiple_possible_frames_above_similarity_threshold_count == len(sequences)

    # Check if frame_correction has removed all sequences for a given gene:
    inframe_seqs_names = set([seq.name.split('-')[-1] for gene, sequences in gene_to_inframe_seq_dictionary.items() for
                             seq in sequences])

    logger.info(f'{"[INFO]:":10} Frame-corrected target file contains at least one sequence for'
                f' {len(inframe_seqs_names)} genes')

    genes_removed = set([name for name in original_gene_names if name not in inframe_seqs_names])
    if genes_removed:
        fill = textwrap.fill(f'{"[WARNING]:":10} After correcting reading frames and filtering out undetermined '
                             f'sequences, the following genes have ALL representative sequences removed:',
                             width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
        logger.warning(fill)
        logger.warning('')

        for gene in genes_removed:
            logger.warning(f'{" " * 10} {gene}')
        logger.warning('')

        if not allow_gene_removal:
            fill = textwrap.fill(f'{"[WARNING]:":10} By default, removing all representative sequences for a gene is '
                                 f'not allowed. If you would like to enable this, please supply the '
                                 f'--allow_gene_removal flag. Exiting now...)',
                                 width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
            logger.warning(fill)
            sys.exit(0)

    # Filter by length if more than one representative sequence for a gene:
    fill = textwrap.fill(f'{"[INFO]:":10} If there is more than one representative sequence for a given gene, '
                         f'sequences less than {filter_by_length_percentage} the length of the longest representative '
                         f'will be removed.',
                         width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
    logger.info(fill)

    logger.info(f'{"[INFO]:":10} Number of in-frame sequences before filtering by minimum length percentage:'
                f' {len(inframe_seqs_to_write)}')

    gene_to_inframe_seq_dictionary_filtered = defaultdict(list)
    for gene, sequences in gene_to_inframe_seq_dictionary.items():
        if len(sequences) > 1:
            max_seq_length = len(max(sequences, key=len))
            for seq in sequences:
                if len(seq) > (max_seq_length * filter_by_length_percentage):
                    gene_to_inframe_seq_dictionary_filtered[gene].append(seq)
        elif len(sequences) == 1:
            gene_to_inframe_seq_dictionary_filtered[gene].append(sequences[0])
        else:
            pass

    filtered_seqs = [seq for gene, sequences in gene_to_inframe_seq_dictionary_filtered.items() for seq in sequences]
    filtered_seqs_names = [seq.name.split('-')[-1] for gene, sequences in
                           gene_to_inframe_seq_dictionary_filtered.items() for seq in sequences]

    logger.info(f'{"[INFO]:":10} Number of in-frame sequences after filtering by minimum length percentage:'
                f' {len(filtered_seqs)}')

    # Remove sequences with low complexity regions:
    low_complexity_seqs = []
    if keep_low_complexity_sequences and low_complexity_seq_names:
        logger.info(f'{"[INFO]:":10} Keeping sequences with regions of low-complexity... ')
    elif low_complexity_seq_names:
        fill = textwrap.fill(f'{"[INFO]:":10} Removing sequences with regions of low-complexity',
                             width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
        logger.info(fill)
        print(f'low_complexity_seq_names is: {low_complexity_seq_names}')

        logger.info(f'{"[INFO]:":10} Number of in-frame sequences before removing sequences with '
                    f'low-complexity regions: {len(filtered_seqs)}')

        gene_to_inframe_seq_dictionary_filtered_complexity = defaultdict(list)
        for gene, sequences in gene_to_inframe_seq_dictionary_filtered.items():
            for seq in sequences:
                if seq.name not in low_complexity_seq_names:
                    gene_to_inframe_seq_dictionary_filtered_complexity[gene].append(seq)
                else:
                    low_complexity_seqs.append(seq)

        filtered_seqs_complexity = [seq for gene, sequences in
                                    gene_to_inframe_seq_dictionary_filtered_complexity.items()
                                    for seq in sequences]
        filtered_seqs_complexity_names = [seq.name.split('-')[-1] for gene, sequences in
                                          gene_to_inframe_seq_dictionary_filtered_complexity.items()
                                          for seq in sequences]

        logger.info(f'{"[INFO]:":10} Number of in-frame sequences after removing sequences with '
                    f'low-complexity regions: {len(filtered_seqs_complexity)}')

        genes_removed = set([name for name in filtered_seqs_names if name not in filtered_seqs_complexity_names])
        if genes_removed:
            fill = textwrap.fill(f'{"[WARNING]:":10} After removing sequences with low-complexity regions, '
                                 f'the following genes have ALL representative sequences removed:',
                                 width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
            logger.warning(fill)
            logger.warning('')

            for gene in genes_removed:
                logger.warning(f'{" " * 10} {gene}')
            logger.warning('')

            if not allow_gene_removal:
                fill = textwrap.fill(f'{"[WARNING]:":10} By default, removing all representative sequences for a gene '
                                     f'is =not allowed. If you would like to enable this, please supply the '
                                     f'--allow_gene_removal flag. Exiting now...)',
                                     width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
                logger.warning(fill)
                sys.exit(0)

    # Write output fasta files:
    if low_complexity_seqs:
        with open(fixed_file_output, 'w') as fixed_handle:
            SeqIO.write(filtered_seqs_complexity, fixed_handle, 'fasta')
        with open(fixed_file_output_low_complexity_sequences, 'w') as low_complexity_handle:
            SeqIO.write(low_complexity_seqs, low_complexity_handle, 'fasta')
    else:
        with open(fixed_file_output, 'w') as fixed_handle:
            SeqIO.write(filtered_seqs, fixed_handle, 'fasta')

    with open(fixed_file_output_with_stop_codons, 'w') as stops_handle:
        SeqIO.write(sequences_with_stop_codons_all_frames, stops_handle, 'fasta')

    with open(fixed_file_output_with_undetermined_frame, 'w') as undetermined_handle:
        SeqIO.write(seqs_with_undetermined_frame, undetermined_handle, 'fasta')

    with open(fixed_file_output_with_multiple_frames_above_maximum_distance, 'w') as maximum_distance_handle:
        SeqIO.write(sequences_with_multiple_possible_frames_above_similarity_threshold, maximum_distance_handle,
                    'fasta')

    if low_complexity_seqs:
        return gene_to_inframe_seq_dictionary_filtered_complexity
    else:
        return gene_to_inframe_seq_dictionary_filtered


def inframe_seq_alignments(gene_to_inframe_seq_dictionary, concurrent_alignments, threads_per_concurrent_alignment):
    """
    Generate per-gene alignments of protein translations of each 'fixed' inframe sequence. Useful for identifying
    problems in the 'fixed' sequences (e.g. when framshifts are present but do not introduce stop codons).

    :param collections.defaultdict gene_to_inframe_seq_dictionary: a dictionary of gene to SeqRecord objects for all
    fixed sequences.
    :param int concurrent_alignments: number of alignments to run concurrently.
    :param int threads_per_concurrent_alignment: Number of threads to run each concurrent alignment with.
    :return:
    """

    output_folder_unaligned = f'01_gene_fixed_seqs_unaligned'
    output_folder_aligned = f'02_gene_fixed_seqs_aligned'
    createfolder(output_folder_unaligned)
    createfolder(output_folder_aligned)

    for gene, seqrecord_list in gene_to_inframe_seq_dictionary.items():
        translated_seqs = []
        for seq in seqrecord_list:
            seq_translated = SeqRecord(Seq(seq.seq.translate()), id=seq.id, name=seq.name, description=seq.description)
            translated_seqs.append(seq_translated)
        with open(f'{output_folder_unaligned}/{gene}_unaligned.fasta', 'w') as unaligned_handle:
            SeqIO.write(translated_seqs, unaligned_handle, 'fasta')

    logger.info(f'{"[INFO]:":10} Running alignments for {len(gene_to_inframe_seq_dictionary)} genes...')

    align_extractions_multiprocessing(output_folder_unaligned,
                                      output_folder_aligned,
                                      concurrent_alignments,
                                      threads_per_concurrent_alignment)


def parse_control_file(control_file):
    """
    Checks that the control file provided exists and can be parsed correctly.

    :param str control_file: path to the control file output of `hybpiper check_targetfile`
    """

    logger.info(f'{"[INFO]:":10} Reading the provided control file "{control_file}"...')

    if not os.path.isfile(control_file) or not file_exists_and_not_empty(control_file):
        logger.error(f'{"[ERROR]:":10} The control file "{control_file}" can not be found or is empty!')
        sys.exit(1)

    expected_parameters = ['TARGETFILE_TYPE', 'TRANSLATE_TARGET_FILE', 'NO_TERMINAL_STOP_CODONS',
                           'SLIDING_WINDOW_SIZE', 'COMPLEXITY_MINIMUM_THRESHOLD', 'ALLOW_GENE_REMOVAL',
                           'LOW_COMPLEXITY_SEQUENCES']

    control_dict = dict()
    with open(control_file, 'r') as control_handle:
        lines = [line for line in control_handle.readlines() if line.rstrip()]
        for line in lines:
            entry = line.split('\t')[0].rstrip()
            if entry not in expected_parameters:
                logger.error(f'{"[ERROR]:":10} The control file "{control_file}" contains unexpected entry'
                             f' {entry}, please check!')
                sys.exit(1)
            else:
                try:
                    entry_value = line.split('\t')[1].rstrip()
                    if len(entry_value) == 0:
                        sys.exit(f'{"[ERROR]:":10} The control file entry {entry} does not have a corresponding value, '
                                 f'please check!')
                    if entry in ['LOW_COMPLEXITY_SEQUENCES']:
                        control_dict[entry] = [seq_name.rstrip() for seq_name in line.split('\t')[1:]]
                    else:
                        control_dict[entry] = entry_value
                except IndexError:
                    sys.exit(f'{"[ERROR]:":10} The control file entry {entry} does not have a corresponding value, '
                             f'please check!')

    logger.info(f'{"[INFO]:":10} Control file parsed!')

    print(control_dict)

    return control_dict


def check_reference_protein_file(reference_protein_file):
    """

    """

    print('TO DO - check ref protein file!')


########################################################################################################################
########################################################################################################################

def standalone():
    """
    Used when this module is run as a stand-alone script. Parses command line arguments and runs function main().:
    """

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    group_1 = parser.add_mutually_exclusive_group(required=True)
    group_1.add_argument('--targetfile_dna', '-t_dna',
                         dest='targetfile_dna',
                         default=False,
                         help='FASTA file containing DNA target sequences for each gene. If there are multiple '
                              'targets for a gene, the id must be of the form: >Taxon-geneName')
    group_1.add_argument('--targetfile_aa', '-t_aa',
                         dest='targetfile_aa',
                         default=False,
                         help='FASTA file containing amino-acid target sequences for each gene. If there are multiple '
                              'targets for a gene, the id must be of the form: >Taxon-geneName')
    parser.add_argument('control_file',
                        help='The *.ctl file output by the command "hybpiper check_targetfile".')
    parser.add_argument('--no_terminal_stop_codons', action='store_true', default=False,
                        help='When testing for open reading frames, do not allow a translated frame to have a single '
                             'stop codon at the C-terminus of the translated protein sequence. Default is False. If '
                             'supplied, this parameter will override the setting in the *.ctl file.')
    parser.add_argument('--allow_gene_removal', action='store_true', default=False,
                        help='Allow frame-correction and filtering steps to remove all representative sequences for a '
                             'given gene. Default is False; script will exit with an information message instead. If '
                             'supplied, this parameter will override the setting in the *.ctl file.')
    parser.add_argument('--reference_protein_file', default=None,
                        help='If a given sequence can be translated in more than one forward frame without stop '
                             'codons, choose the translation that best matches the corresponding reference protein '
                             'provided in this fasta file. Sequence ids must be of the form: >Taxon-geneName')
    parser.add_argument('--maximum_distance', default=0.5, type=float, metavar='FLOAT',
                        help='When comparing possible translation frames to a reference protein, the maximum distance '
                             'allowed between the translated frame and the reference sequence for any possible '
                             'translation frame to be selected. Useful to filter out sequences with frameshifts that '
                             'do NOT introduce stop codons. 0.0 means identical sequences, 1.0 means completely '
                             'different sequences. Default is 0.5')
    parser.add_argument('--filter_by_length_percentage', default=0.0, type=float, metavar='FLOAT',
                        help='If more than one sequence is present for a given gene, only include sequences longer '
                             'than this percentage of the longest gene sequence. Default is 0.0 (all sequences '
                             'retained)')
    parser.add_argument('--keep_low_complexity_sequences', action='store_true', default=False,
                        help='Keep sequences that contain regions of low-complexity, as identified by the command '
                             '"hybpiper check_targetfile". Default is to remove these sequences.')
    parser.add_argument('--alignments', action='store_true', default=False,
                        help='Create per-gene alignments for in-frame translated protein sequences.')
    parser.add_argument('--concurrent_alignments', default=1, type=int, metavar='INTEGER',
                        help='Number of alignments to run concurrently. Default is 1.')
    parser.add_argument('--threads_per_concurrent_alignment', default=1, type=int, metavar='INTEGER',
                        help='Number of threads to run each concurrent alignment with. Default is 1.')

    results = parser.parse_args()
    return results


########################################################################################################################
########################################################################################################################
# Run script:

def main(args):
    """
    Entry point for the assemble.py module.

    :param argparse.Namespace args:
    """

    # Check targefile matches control file tupe or exit

    logger.info(f'{"[INFO]:":10} Script was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.info(f'{fill}\n')

    # Check for any external executables:
    logger.info(f'{"[INFO]:":10} Checking for external dependencies:')
    executables = ['mafft']

    for exe in executables:
        exe_loc = shutil.which(exe)
        if exe_loc:
            logger.info(f'{"[INFO]:":10} {exe} found at {exe_loc}')
        else:
            logger.info(f'{"[ERROR]:":10}{exe:20} not found in your $PATH!')
            sys.exit(1)

    # Parse the control file
    control_dict = parse_control_file(args.control_file)

    # Set target file type and path, check it exists and isn't empty, and check type matches control file:
    if args.targetfile_dna:
        targetfile = os.path.abspath(args.targetfile_dna)
        targetfile_type = 'DNA'
    elif args.targetfile_aa:
        targetfile = os.path.abspath(args.targetfile_aa)
        targetfile_type = 'protein'

    expected_targetfile_type = control_dict['TARGETFILE_TYPE']
    if targetfile_type != expected_targetfile_type:
        fill = textwrap.fill(f'{"[ERROR]:":10} The target file type you have provided ({targetfile_type}) does not '
                             f'match the type in your control file ({expected_targetfile_type}). Please run '
                             f'"hybpiper fix_targetfile" with the same targetfile you used as input to '
                             f'"hybpiper check_targetfile"!',
                             width=90, subsequent_indent=' ' * 11)
        logger.error(fill)
        sys.exit(1)

    if os.path.isfile(targetfile) and not os.path.getsize(targetfile) == 0:
        logger.debug(f'Input target file {os.path.basename(targetfile)} exists and is not empty, proceeding...')
    else:
        sys.exit(f'Input target file {os.path.basename(targetfile)} does not exist or is empty!')

    logger.debug(f'The target file {os.path.basename(targetfile)} has been provided, containing {targetfile_type} '
                 f'sequences')

    # Check the reference protein file, if provided:
    check_reference_protein_file(args.reference_protein_file)

    # Get parameters from command line if provided, control file if not:
    no_terminal_stop_codons = args.no_terminal_stop_codons if \
        args.no_terminal_stop_codons else True if control_dict['NO_TERMINAL_STOP_CODONS'] == 'True' else False

    allow_gene_removal = args.allow_gene_removal if \
        args.allow_gene_removal else True if control_dict['ALLOW_GENE_REMOVAL'] == 'True' else False

    low_complexity_seq_names = control_dict['LOW_COMPLEXITY_SEQUENCES'] if control_dict['LOW_COMPLEXITY_SEQUENCES'][0] \
                                                                           != 'None' else ''

    logger.info(f'{"[INFO]:":10} Terminal stop codon not allowed: {no_terminal_stop_codons}')
    logger.info(f'{"[INFO]:":10} Allow removal of all representative sequences for a given gene:'
                f' {allow_gene_removal}')
    logger.info(f'{"[INFO]:":10} Number of sequences with low-complexity regions: {len(low_complexity_seq_names)}')

    if targetfile_type == 'DNA':
        gene_to_inframe_seq_dictionary = get_inframe_sequence(targetfile,
                                                              no_terminal_stop_codons,
                                                              args.reference_protein_file,
                                                              args.maximum_distance,
                                                              args.filter_by_length_percentage,
                                                              allow_gene_removal,
                                                              args.keep_low_complexity_sequences,
                                                              low_complexity_seq_names)

        if args.alignments:
            inframe_seq_alignments(gene_to_inframe_seq_dictionary,
                                   args.concurrent_alignments,
                                   args.threads_per_concurrent_alignment)


########################################################################################################################
########################################################################################################################

if __name__ == '__main__':
    if not len(sys.argv) >= 1:
        print(__doc__)
        sys.exit()
    sys.exit(standalone())

########################################################################################################################
########################################################################################################################



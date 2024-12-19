import os
import pysam
import pytest
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import emu


@pytest.mark.parametrize(
    "output_dir, output_basename, input_file, expected_output_basepath",
    [
        (
            "results",  # output_dir
            "",  # output_basename
            "some_file.fastq",  # input_file
            "results/some_file.fastq",  # expected_output_basepath
        ),
        (
            "results",  # similar as above ...
            "",
            ["some_file.fastq", "another_file.fastq"],
            "results/some_file-another_file",
        ),
        (
            "results",
            "my-output-basename",
            "some_file.fastq",
            "results/my-output-basename",
        ),
    ],
)
def test_create_output_basepath(
    output_dir, output_basename, input_file, expected_output_basepath
):
    got_output_basepath = emu.create_output_basepath(
        output_dir, output_basename, input_file
    )
    assert got_output_basepath == expected_output_basepath


@pytest.mark.parametrize(
    "in_file_list, database, preset, num_threads, max_sec_alignments, minibatch_size, optional_args, output_file, expected_mm2_command",
    [
        (
            ["foo1.fastq", "foo2.fastq"],
            "emu_database",
            "map-ont",
            "4",
            "5",
            "500M",
            "",
            "out.sam",
            "minimap2 -ax map-ont -t 4 -N 5 -p .9 -K 500M emu_database/species_taxid.fasta foo1.fastq foo2.fastq -o out.sam",
        ),
        (
            ["foo1.fastq", "foo2.fastq"],
            "custom_db",
            "map-ont",
            "10",
            "1",
            "100M",
            True,
            "myoutput.sam",
            "minimap2 -ax map-ont -t 10 -N 1 -p .9 -u f -K 100M custom_db/species_taxid.fasta foo1.fastq foo2.fastq -o myoutput.sam",
        ),
    ],
)
def test_create_minimap2_command(
    in_file_list,
    database,
    preset,
    num_threads,
    max_sec_alignments,
    minibatch_size,
    optional_args,
    output_file,
    expected_mm2_command,
):
    got_mm2_command = emu.create_minimap2_command(
        in_file_list,
        database,
        preset,
        num_threads,
        max_sec_alignments,
        minibatch_size,
        optional_args,
        output_file,
    )
    assert expected_mm2_command == got_mm2_command


@pytest.mark.parametrize(
    "file_path, expected_bool",
    [
        ("myfile.sam", True),
        ("myfile.txt", False),
        ("myfile.sam.py", False),
    ],
)
def test_is_sam_file(file_path, expected_bool):
    is_sam = emu.is_sam_file(file_path)
    assert is_sam == expected_bool


@pytest.mark.parametrize(
    "file_path, expected_bool",
    [
        ("myfile.fastq", True),
        ("myfile.fq", True),
        ("myfile.fastq.txt", False),
    ],
)
def test_is_fastq_file(file_path, expected_bool):
    is_fastq = emu.is_fastq_file(file_path)
    assert is_fastq == expected_bool


def test_generate_alignments():
    sam_file = "results/somesamfile.sam"

    got_alignment_file = emu.generate_alignments(
        in_file_list=[sam_file],
        out_basename="mybasepath",
        database="emu_database",
        mm2_preset="map-ont",
        num_threads="4",
        max_sec_alignments="4",
        minibatch_size="500M",
        mm2_fwd_only=False,
    )

    assert got_alignment_file == sam_file


def test_get_cigar_op_log_probabilities():
    pass


@pytest.mark.parametrize(
    "cigar_string, expected_alignment_length",
    [
        ("2M1000X", 2),
        ("2M3I2M", 7),
        ("2M4D2M", 8),
        ("10M1000X10M", 20),
        ("2S4M2S", 8),
    ],
)
def test_get_align_len(cigar_string, expected_alignment_length):
    mock_alignment = pysam.AlignedSegment()
    mock_alignment.cigarstring = cigar_string

    got_alignment_length = emu.get_align_len(mock_alignment)

    assert got_alignment_length == expected_alignment_length


@pytest.mark.parametrize(
    "cigar_string, nm_tag, expected_align_stats",
    [
        ("1M2I3D4S5X", 2 + 3 + 4 + 5, [2, 3, 4, 9]),
        ("3S1I1M2D", 1 + 2 + 3, [1, 2, 3, 3]),
    ],
)
def test_get_align_stats(cigar_string, nm_tag, expected_align_stats):
    mock_alignment = pysam.AlignedSegment()
    mock_alignment.cigarstring = cigar_string
    mock_alignment.set_tag("NM", nm_tag)

    got_align_stats = emu.get_align_stats(mock_alignment)

    assert got_align_stats == expected_align_stats

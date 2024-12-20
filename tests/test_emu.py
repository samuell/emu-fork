import os
import pysam
import pytest
import sys
import tempfile

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


@pytest.mark.parametrize(
    "cigar_string, expected_alignment_length",
    [
        ("2M1000X", 2),
        ("2M3I2M", 7),
        ("2M4D2M", 8),
        ("10M1000X10M", 20),
        ("2S4M2S", 8),
        ("18S1154M", 1172),
        ("198M1D974M", 1173),
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
def test_get_misalign_type_counts(cigar_string, nm_tag, expected_align_stats):
    mock_alignment = pysam.AlignedSegment()
    mock_alignment.cigarstring = cigar_string
    mock_alignment.set_tag("NM", nm_tag)

    got_align_stats = emu.get_misalign_type_counts(mock_alignment)

    assert got_align_stats == expected_align_stats


def test_get_misalign_type_log_probs():
    with tempfile.NamedTemporaryFile(suffix=".sam") as samfile:
        samfile.write(bytearray(SAM_EXAMPLE_CONTENT, "utf-8"))
        samfile.flush()

        (
            misalign_type_log_probs,
            zero_locs,
            longest_alignments,
        ) = emu.get_misalign_type_log_probs(samfile.name)

    assert len(misalign_type_log_probs) + len(zero_locs) == 4

    for logprob in misalign_type_log_probs:
        assert -1 < logprob < 0

    assert "Sphingobacterium_puteal_r1" in longest_alignments
    assert "Mycobacterium_saskatchewanense_r1" in longest_alignments
    assert "Streptococcus_sobrinus_r1" in longest_alignments

    assert longest_alignments["Sphingobacterium_puteal_r1"] == 1248
    assert longest_alignments["Mycobacterium_saskatchewanense_r1"] == 1172
    assert longest_alignments["Streptococcus_sobrinus_r1"] == 1360


SAM_EXAMPLE_CONTENT = """@SQ	SN:2420510:emu_db:1	LN:1451
@SQ	SN:1933220:emu_db:46868	LN:1537
@SQ	SN:141349:emu_db:16952	LN:1452
@SQ	SN:185642:emu_db:18546	LN:1533
@SQ	SN:1310:emu_db:16914	LN:1538
@SQ	SN:1310:emu_db:9243	LN:1459
@SQ	SN:1310:emu_db:29451	LN:1560
@PG	ID:minimap2	PN:minimap2	VN:2.24-r1122	CL:minimap2 -ax map-ont -t 3 -N 50 -p .9 -K 500000000 -o ./results/full_length_emu_alignments.sam emu_database/species_taxid.fasta example/full_length.fa
Sphingobacterium_puteal_r1	0	2420510:emu_db:1	84	60	1248M	*	0	0	*	*	NM:i:0	ms:i:2496	AS:i:2496	nn:i:0	tp:A:P	cm:i:210	s1:i:1191	s2:i:1073	de:f:0	rl:i:211
Sphingobacterium_puteal_r1	256	1933220:emu_db:46868	115	0	1248M	*	0	0	*	*	NM:i:16	ms:i:1664	AS:i:2400	nn:i:0	tp:A:S	cm:i:173	s1:i:1073	de:f:0.0128	rl:i:211
Mycobacterium_saskatchewanense_r1	0	141349:emu_db:16952	1	0	18S1154M	*	0	0	*	*	NM:i:16	ms:i:1476	AS:i:2212	nn:i:0	tp:A:S	cm:i:136	s1:i:927	de:f:0.0139	rl:i:277
Mycobacterium_saskatchewanense_r1	256	185642:emu_db:18546	6	0	198M1D974M	*	0	0	*	*	NM:i:18	ms:i:1408	AS:i:2236	nn:i:0	tp:A:S	cm:i:142	s1:i:948	de:f:0.0153	rl:i:277
Streptococcus_sobrinus_r1	0	1310:emu_db:16914	81	0	1360M	*	0	0	*	*	NM:i:0	ms:i:2720	AS:i:2720	nn:i:0	tp:A:P	cm:i:208	s1:i:1179	s2:i:1179	de:f:0	rl:i:373
Streptococcus_sobrinus_r1	256	1310:emu_db:9243	66	0	1360M	*	0	0	*	*	NM:i:0	ms:i:2720	AS:i:2720	nn:i:0	tp:A:S	cm:i:208	s1:i:1179	de:f:0	rl:i:373
Streptococcus_sobrinus_r1	272	1310:emu_db:29451	113	0	1360M	*	0	0	*	*	NM:i:3	ms:i:2564	AS:i:2702	nn:i:0	tp:A:S	cm:i:202	s1:i:1162	de:f:0.0022	rl:i:373
"""
